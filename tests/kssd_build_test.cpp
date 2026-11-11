#include "KssdBuild.h"
#include "KssdIndex.h"
#include "io/KssdFastx.h"
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <unistd.h>

using namespace Sketch;
namespace KB = Sketch::KssdBuild;

static void check(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

// Independent, non-rolling oracle: explicitly encode each complete window,
// construct its reverse complement, and read the original permutation.
static std::vector<uint64_t> oracle(const std::vector<std::string>& records,
                                    const kssd_parameter_t& p) {
    std::vector<uint64_t> result;
    for (const auto& sequence : records) {
        for (size_t i = 0; i + p.kmer_size <= sequence.size(); ++i) {
            uint64_t forward = 0, reverse = 0;
            bool valid = true;
            for (unsigned j = 0; j < p.kmer_size; ++j) {
                unsigned code;
                switch (sequence[i + j]) {
                case 'A': case 'a': code = 0; break;
                case 'C': case 'c': code = 1; break;
                case 'G': case 'g': code = 2; break;
                case 'T': case 't': code = 3; break;
                default: valid = false; code = 0;
                }
                forward = (forward << 2) | code;
                reverse |= uint64_t(3 - code) << (2 * j);
            }
            if (!valid) continue;
            const uint64_t canonical = std::min(forward, reverse);
            const uint32_t dimension = (canonical & p.domask) >> (p.half_outctx_len * 2);
            const int rank = p.shuffled_dim[dimension];
            if (rank < p.dim_start || rank >= p.dim_end) continue;
            result.push_back((((canonical & p.undomask0) |
                ((canonical & p.undomask1) << (p.kmer_size * 2 - p.half_outctx_len * 4)))
                >> (p.drlevel * 4)) | uint64_t(rank - p.dim_start));
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

static std::vector<uint64_t> values(const Kssd& sketch) {
    if (sketch.get_half_k() - sketch.get_drlevel() > 8) return sketch.getHashes64();
    const auto& v = sketch.getHashes();
    return {v.begin(), v.end()};
}

static std::string randomDNA(size_t length) {
    std::mt19937 random(431);
    std::string result(length, 'A');
    for (char& base : result) base = "ACGT"[random() & 3];
    return result;
}

static void checkStream(const std::string& sequence,
                        std::shared_ptr<const kssd_parameter_t> p) {
    const auto expected = oracle({sequence}, *p);
    for (auto lookup : {KB::Lookup::Bitmap, KB::Lookup::Map, KB::Lookup::Direct}) {
        for (size_t chunk : {size_t(1), size_t(31), size_t(64), size_t(4093), sequence.size()}) {
            KB::Stream stream(*p, lookup);
            KB::Accumulator<uint64_t> accumulator;
            accumulator.reset();
            for (size_t i = 0; i < sequence.size(); i += chunk)
                stream.append(sequence.data() + i, std::min(chunk, sequence.size() - i),
                    [&](uint64_t hash) { accumulator.insert(hash); });
            std::vector<uint64_t> observed;
            accumulator.finish(observed);
            check(observed == expected, "stream/oracle mismatch");
        }
    }
    const std::vector<std::string> records = {sequence.substr(0, 16003),
        sequence.substr(16003, 32009), sequence.substr(48012)};
    Kssd sketch(p);
    std::vector<std::string> prefix;
    for (const auto& record : records) {
        sketch.update(record.data(), record.size());
        prefix.push_back(record);
        check(values(sketch) == oracle(prefix, *p), "query then update lost tuples");
    }
    const auto before = values(sketch);
    Kssd loaded(p);
    if (p->half_k - p->drlevel > 8) loaded.loadHashes64(before);
    else loaded.loadHashes(std::vector<uint32_t>(before.begin(), before.end()));
    loaded.update(sequence.data(), sequence.size());
    prefix.push_back(sequence);
    check(values(loaded) == oracle(prefix, *p), "load then update lost tuples");
    check(values(sketch.project(p->drlevel)) == before, "same-resolution projection mismatch");
    const auto lite = sketch.toLite();
    check(lite.hashList.size() + lite.hashList64.size() == before.size(), "toLite missed pending data");
    sketch.finalize(); sketch.finalize();
    check(values(sketch) == before && sketch.hashSet.empty() && sketch.hashSet64.empty(),
          "finalize changed hashes or retained build set");
    bool threw = false;
    try { sketch.update("ACGT", 4); } catch (const std::logic_error&) { threw = true; }
    check(threw, "sealed sketch accepted update");
    Kssd legacy(p);
    legacy.update("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    check(legacy.hashList.size() + legacy.hashList64.size() == values(legacy).size(),
          "legacy update left public lists stale");
}

static void checkFiles(const std::string& sequence, const kssd_parameter_t& p) {
    char path[] = "/tmp/rabbitsketch-kssd-XXXXXX";
    check(mkdtemp(path) != nullptr, "temporary directory failed");
    const std::string directory(path);
    const std::vector<std::string> records = {sequence.substr(0, 19003), "AC", sequence.substr(19003)};
    const auto expected = oracle(records, p);
    std::string fasta, fastq;
    for (size_t r = 0; r < records.size(); ++r) {
        fasta += ">record" + std::to_string(r) + "\r\n";
        for (size_t i = 0; i < records[r].size(); i += 79)
            fasta += records[r].substr(i, 79) + "\r\n";
        fastq += "@record" + std::to_string(r) + "\n" + records[r] + "\n+\n" +
            std::string(records[r].size(), 'I') + "\n";
    }
    fasta.resize(fasta.size() - 2); // final line without a newline
    for (int format = 0; format < 4; ++format) {
        const auto file = directory + "/" + std::to_string(format);
        const std::string& text = format & 1 ? fastq : fasta;
        if (format & 2) {
            gzFile output = gzopen(file.c_str(), "wb");
            check(output && gzwrite(output, text.data(), text.size()) == int(text.size()), "gzip write failed");
            check(gzclose(output) == Z_OK, "gzip close failed");
        } else { std::ofstream output(file, std::ios::binary); output << text; }
        KB::Stream stream(p);
        KB::Accumulator<uint64_t> accumulator; accumulator.reset();
        KB::FastxReader reader;
        check(reader.read(file, stream, [&](uint64_t hash) { accumulator.insert(hash); }), "FASTX read failed");
        std::vector<uint64_t> observed; accumulator.finish(observed);
        check(observed == expected, "FASTX changed record or line boundaries");
    }
    for (int format = 0; format < 4; ++format)
        std::remove((directory + "/" + std::to_string(format)).c_str());
    rmdir(directory.c_str());
}

template<typename Key>
static void checkIndex(unsigned bits, std::vector<std::vector<Key>> keys) {
    std::map<uint64_t, std::vector<uint32_t>> expected;
    for (size_t i = 0; i < keys.size(); ++i)
        for (Key hash : keys[i]) expected[hash].push_back(i);
    for (int threads : {1, 4, 8}) {
        const auto index = KB::buildIndex<Key>(keys.size(), threads, bits,
            [&](int i) -> const std::vector<Key>& { return keys[i]; });
        size_t total = 0;
        for (const auto& pair : expected) {
            const uint64_t hash = pair.first;
            size_t start = 0, count = 0;
            if (index.offsets.empty()) {
                auto found = index.sparse.find(hash);
                if (found != index.sparse.end()) { start = found->second.off; count = found->second.cnt; }
            } else { start = index.offsets[hash]; count = index.offsets[hash + 1] - start; }
            if (pair.second.size() == 1) check(count == 0, "singleton retained");
            else {
                check(count == pair.second.size(), "posting count mismatch");
                check(std::equal(pair.second.begin(), pair.second.end(), index.postings.begin() + start),
                      "posting IDs mismatch");
                total += count;
            }
        }
        check(total == index.postings.size(), "extra postings");
    }
}

int main() {
    try {
        auto sequence = randomDNA(80000);
        for (unsigned byte = 0; byte < 256; ++byte) sequence[100 + byte * 233] = static_cast<char>(byte);
        auto p = getKssdParameters();
        check(p->shuffled_map.size() == 4096 && p->shuffled_map.capacity() < 16384,
              "sparse rank map still overallocated");
        check(p->admitted_dimensions.size() * sizeof(uint64_t) == (1 << 21), "bitmap size mismatch");
        checkStream(sequence, p);
        checkStream(sequence, getKssdParameters(10, 6, 2));
        checkStream(sequence, getKssdParameters(10, 4, 1));
        checkStream(sequence, getKssdParameters(16, 3, 0));
        auto dense = getKssdParameters(10, 6, 2);
        Kssd left(dense), right(p);
        left.update(sequence.data(), sequence.size()); right.update(sequence.data(), sequence.size());
        check(values(left.project(3)) == values(right), "sparser projection changed sampled set");
        check(left.jaccard(&right) == 1.0, "mixed-resolution query mismatch");
        auto projected = left.project(3);
        const auto extra = randomDNA(120001);
        projected.update(extra.data(), extra.size()); right.update(extra.data(), extra.size());
        check(values(projected) == values(right), "project then update mismatch");
        std::vector<std::shared_ptr<const kssd_parameter_t>> cached(8);
        #pragma omp parallel for num_threads(4)
        for (int i = 0; i < 8; ++i) cached[i] = getKssdParameters();
        for (const auto& copy : cached) check(copy.get() == p.get(), "parameters not shared across builders");
        auto file_sequence = randomDNA(80000);
        for (size_t i = 31; i < file_sequence.size(); i += 127) file_sequence[i] = "SDWsdwN"[i % 7];
        checkFiles(file_sequence, *p);
        KB::Accumulator<uint32_t> repetitive; repetitive.reset();
        for (int i = 0; i < 1000000; ++i) repetitive.insert(i % 37);
        std::vector<uint32_t> unique; repetitive.finish(unique);
        check(repetitive.usingHashSet() && unique.size() == 37, "duplicate spill failed");
        checkIndex<uint32_t>(12, {{0, 1, 17, 4095}, {}, {0, 17}, {1, 17, 33}});
        checkIndex<uint32_t>(32, {{0, UINT32_MAX}, {UINT32_MAX}, {17}});
        checkIndex<uint64_t>(64, {{0, UINT64_MAX}, {UINT64_MAX}, {17}});
        checkIndex<uint32_t>(28, {{}, {}});
        std::cout << "KSSD build, cache, streaming, projection, FASTX and CSR checks passed\n";
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
