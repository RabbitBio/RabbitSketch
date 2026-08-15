#include "api/SketchBuilder.h"
#include "api/SketchConfig.h"
#include "io/FastxReader.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>
#include <zlib.h>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Exception, typename Function>
void expectThrows(Function function, const std::string& message) {
    try {
        function();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error(message);
}

struct TemporaryFiles {
    std::string directory;
    TemporaryFiles() {
        char pattern[] = "/tmp/rabbitsketch-fastx-test.XXXXXX";
        char* path = ::mkdtemp(pattern);
        if (path == nullptr) throw std::runtime_error("mkdtemp failed");
        directory = path;
    }
    ~TemporaryFiles() {
        for (const char* name : {"input.fa", "input.fa.gz", "quality.fq",
                                 "abundance.fa", "broken.fq"})
            ::unlink((directory + "/" + name).c_str());
        ::rmdir(directory.c_str());
    }
};

void writeText(const std::string& path, const std::string& content) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output << content;
    if (!output) throw std::runtime_error("failed to write test FASTX");
}

void writeGzip(const std::string& path, const std::string& content) {
    gzFile output = gzopen(path.c_str(), "wb");
    if (output == nullptr ||
        gzwrite(output, content.data(), static_cast<unsigned>(content.size())) !=
            static_cast<int>(content.size())) {
        if (output != nullptr) gzclose(output);
        throw std::runtime_error("failed to write gzip test FASTX");
    }
    gzclose(output);
}

Sketch::API::SketchConfig fastConfig(uint16_t k = 5) {
    Sketch::API::SketchConfig config;
    config.algorithm = Sketch::API::Algorithm::FastKMV;
    config.sampling = Sketch::API::SamplingMode::BottomK;
    config.kmer_size = k;
    config.resolution = 64;
    return config;
}

void testReaderAndMultiSketchBuilder() {
    TemporaryFiles files;
    const std::string fasta =
        ">alpha first record\nACGTACGT\nACGT\n"
        ">beta\nTTTTNNNN\nacgt\n";
    const std::string fasta_path = files.directory + "/input.fa";
    const std::string gzip_path = files.directory + "/input.fa.gz";
    writeText(fasta_path, fasta);
    writeGzip(gzip_path, fasta);

    Sketch::IO::FastxReader reader(gzip_path);
    Sketch::IO::FastxRecord record;
    check(reader.next(record) && record.name == "alpha" &&
              record.sequence == "ACGTACGTACGT" && !record.hasQuality(),
          "gzip FASTA reader lost multiline sequence/name state");
    check(reader.next(record) && record.name == "beta" &&
              !reader.next(record) &&
              reader.format() == Sketch::IO::FastxFormat::Fasta &&
              reader.stats().records == 2,
          "FASTX reader EOF/format/stats are inconsistent");

    Sketch::API::SketchConfig fast = fastConfig();
    Sketch::API::SketchConfig frac = fast;
    frac.algorithm = Sketch::API::Algorithm::FracMinHash;
    frac.sampling = Sketch::API::SamplingMode::Scaled;
    frac.scaled = 2;
    frac.aggregation = Sketch::API::AggregationMode::OneSketchPerRecord;
    const auto plain = Sketch::API::buildFastx(fasta_path, {fast, frac});
    const auto compressed = Sketch::API::buildFastx(gzip_path, {fast, frac});
    check(plain.size() == 3 && compressed.size() == 3,
          "mixed aggregation did not emit one file sketch plus two record sketches");
    check(plain[0].stats.records == 2 && plain[0].stats.input_bases == 24 &&
              plain[0].stats.sealed,
          "file-level FASTX BuildStats are incomplete");
    check(plain[0].sketch.query(compressed[0].sketch).jaccard == 1.0,
          "plain and gzip FASTX builds differ");
    check(plain[1].record_name == "alpha" &&
              plain[2].record_name == "beta",
          "record-level build provenance was not preserved");
}

void testQualityAmbiguityHpcAndAbundance() {
    TemporaryFiles files;
    const std::string fastq_path = files.directory + "/quality.fq";
    writeText(fastq_path,
              "@quality\nACGTACGTAC\n+\nIIII!IIIII\n"
              "@ambiguous\nAAANAAAAAA\n+\nIIIIIIIIII\n");
    Sketch::API::SketchConfig quality = fastConfig(3);
    quality.minimum_base_quality = 20;
    const auto quality_result = Sketch::API::buildFastx(fastq_path, {quality});
    check(quality_result.size() == 1 &&
              quality_result[0].stats.candidate_kmers == 16 &&
              quality_result[0].stats.accepted_kmers == 10 &&
              quality_result[0].stats.skipped_low_quality_kmers == 3 &&
              quality_result[0].stats.skipped_ambiguous_kmers == 3,
          "FASTQ quality/ambiguous k-mer accounting is incorrect");

    Sketch::API::SketchConfig hpc = fastConfig(3);
    hpc.homopolymer_compressed = true;
    Sketch::IO::FastxRecord hpc_record;
    hpc_record.name = "hpc";
    hpc_record.sequence = "AAAAACCC";
    Sketch::API::MultiSketchBuilder hpc_builder({hpc});
    hpc_builder.update(hpc_record);
    const auto hpc_result = hpc_builder.finish("hpc");
    check(hpc_result[0].stats.input_bases == 8 &&
              hpc_result[0].stats.candidate_kmers == 0,
          "homopolymer compression was not applied before k-mer accounting");

    const std::string abundance_path = files.directory + "/abundance.fa";
    writeText(abundance_path, ">abundance\nAAAAAA\n");
    Sketch::API::SketchConfig abundance;
    abundance.algorithm = Sketch::API::Algorithm::ProbMinHash;
    abundance.sampling = Sketch::API::SamplingMode::AlgorithmNative;
    abundance.kmer_size = 3;
    abundance.resolution = 64;
    abundance.track_abundance = true;
    abundance.weight_semantics = Sketch::Rank::WeightSemantics::Frequency;
    abundance.min_abundance = 2;
    const auto abundance_result = Sketch::API::buildFastx(
        abundance_path, {abundance});
    check(abundance_result[0].stats.accepted_kmers == 4 &&
              abundance_result[0].stats.distinct_kmers == 1 &&
              abundance_result[0].stats.accepted_distinct_kmers == 1 &&
              abundance_result[0].sketch.metadata().weight_semantics ==
                  Sketch::Rank::WeightSemantics::Frequency,
          "frequency-aware ProbMinHash build did not aggregate canonical counts");

    Sketch::API::SketchConfig reject = quality;
    reject.ambiguous_policy = Sketch::API::AmbiguousPolicy::RejectRecord;
    expectThrows<std::invalid_argument>(
        [&]() { (void)Sketch::API::buildFastx(fastq_path, {reject}); },
        "RejectRecord policy accepted an ambiguous FASTQ record");

    const std::string broken_path = files.directory + "/broken.fq";
    writeText(broken_path, "@broken\nACGT\n+\nII\n");
    expectThrows<std::invalid_argument>(
        [&]() {
            Sketch::IO::FastxReader broken(broken_path);
            Sketch::IO::FastxRecord value;
            (void)broken.next(value);
        },
        "FASTX reader accepted a truncated FASTQ quality string");
}

void testAllAlgorithmFanout() {
    const std::string sequence =
        "ACGTACGTTGCAACGTACGTTGCAACGTACGTTGCAACGTACGTTGCA";
    std::vector<Sketch::API::SketchConfig> configs;

    Sketch::API::SketchConfig legacy = fastConfig(5);
    legacy.algorithm = Sketch::API::Algorithm::LegacyMinHash;
    configs.push_back(legacy);
    configs.push_back(fastConfig(5));

    Sketch::API::SketchConfig frac = fastConfig(5);
    frac.algorithm = Sketch::API::Algorithm::FracMinHash;
    frac.sampling = Sketch::API::SamplingMode::Scaled;
    frac.scaled = 2;
    configs.push_back(frac);

    Sketch::API::SketchConfig hll;
    hll.algorithm = Sketch::API::Algorithm::HyperLogLog;
    hll.sampling = Sketch::API::SamplingMode::RegisterPrecision;
    hll.kmer_size = 5;
    hll.resolution = 6;
    configs.push_back(hll);

    Sketch::API::SketchConfig set = hll;
    set.algorithm = Sketch::API::Algorithm::SetSketch;
    configs.push_back(set);

    Sketch::API::SketchConfig kssd;
    kssd.algorithm = Sketch::API::Algorithm::KSSD;
    kssd.sampling = Sketch::API::SamplingMode::KssdReduction;
    kssd.kmer_size = 20;
    kssd.resolution = 0;
    kssd.kssd_half_subk = 6;
    configs.push_back(kssd);

    Sketch::API::SketchConfig prob;
    prob.algorithm = Sketch::API::Algorithm::ProbMinHash;
    prob.sampling = Sketch::API::SamplingMode::AlgorithmNative;
    prob.kmer_size = 5;
    prob.resolution = 64;
    configs.push_back(prob);

    Sketch::API::SketchConfig bindash = prob;
    bindash.algorithm = Sketch::API::Algorithm::BinDash;
    configs.push_back(bindash);

    Sketch::API::SketchConfig order = prob;
    order.algorithm = Sketch::API::Algorithm::OrderMinHash;
    order.aggregation = Sketch::API::AggregationMode::OneSketchPerRecord;
    order.ambiguous_policy = Sketch::API::AmbiguousPolicy::RejectRecord;
    order.order_l = 2;
    order.order_m = 32;
    configs.push_back(order);

    Sketch::IO::FastxRecord record;
    record.name = "fanout";
    record.sequence = sequence;
    Sketch::API::MultiSketchBuilder builder(configs);
    builder.update(record);
    const auto built = builder.finish("fanout", "memory", record.name);
    check(built.size() == configs.size(),
          "single-pass fanout omitted an algorithm lane");
    for (const auto& item : built) {
        Sketch::Query::Request request;
        if (item.config.algorithm == Sketch::API::Algorithm::OrderMinHash)
            request.metric = Sketch::Query::Metric::OrderSimilarity;
        const Sketch::Query::Result identity = item.sketch.query(
            item.sketch, request);
        check(identity.estimate_available && identity.value == 1.0 &&
                  item.stats.records == 1 && item.stats.sealed &&
                  item.sketch.algorithm() == item.config.algorithm,
              std::string("fanout identity query failed for ") +
                  Sketch::API::algorithmName(item.config.algorithm));
    }
}

} // namespace

int main() {
    try {
        testReaderAndMultiSketchBuilder();
        testQualityAmbiguityHpcAndAbundance();
        testAllAlgorithmFanout();
        std::cout << "fastx_builder_test: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "fastx_builder_test: FAIL: " << error.what() << '\n';
        return 1;
    }
}
