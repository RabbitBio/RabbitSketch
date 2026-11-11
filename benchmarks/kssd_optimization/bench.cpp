#include "Sketch.h"
#include "api/SketchBuilder.h"
#include "io/FastxReader.h"
#include "api/OpenMPCompat.h"
#if __has_include("KssdBuild.h")
#define OPTIMIZED_KSSD 1
#include "KssdBuild.h"
#include "io/KssdFastx.h"
#endif
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

static double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "usage: kssd_bench list mode threads repeats fragment dump warmups\n";
        return 2;
    }
    const std::string mode = argv[2];
    const int threads = std::stoi(argv[3]), repeats = std::stoi(argv[4]);
    const size_t fragment = std::stoull(argv[5]);
    const int warmups = std::stoi(argv[7]);
    std::ifstream listing(argv[1]);
    if (!listing) throw std::runtime_error("cannot open input list");
    std::vector<std::string> paths;
    for (std::string path; std::getline(listing, path);) if (!path.empty()) paths.push_back(path);
    std::vector<std::vector<std::string>> records(paths.size());
    uint64_t bases = 0, record_count = 0;
    for (size_t i = 0; i < paths.size(); ++i) {
        Sketch::IO::FastxReader reader(paths[i]); Sketch::IO::FastxRecord record;
        while (reader.next(record)) {
            bases += record.sequence.size();
            if (fragment) {
                for (size_t pos = 0; pos < record.sequence.size(); pos += fragment)
                    records[i].push_back(record.sequence.substr(pos, fragment));
            } else records[i].push_back(std::move(record.sequence));
        }
        record_count += records[i].size();
    }
    const double setup_start = now();
    std::shared_ptr<const Sketch::kssd_parameter_t> parameters;
    if (mode != "builder") {
#ifdef OPTIMIZED_KSSD
        parameters = Sketch::getKssdParameters();
#else
        parameters = std::make_shared<Sketch::kssd_parameter_t>();
#endif
    }
    const double setup = now() - setup_start;
    std::cout << std::setprecision(12);
    for (int repeat = -warmups; repeat < repeats; ++repeat) {
        std::vector<std::vector<uint32_t>> hashes(paths.size());
        std::vector<std::vector<Sketch::API::BuildResult>> built(paths.size());
        const double start = now();
        #pragma omp parallel for num_threads(threads) schedule(dynamic)
        for (size_t i = 0; i < paths.size(); ++i) {
            if (mode == "builder") {
                Sketch::API::SketchConfig config;
                config.algorithm = Sketch::API::Algorithm::KSSD;
                config.sampling = Sketch::API::SamplingMode::KssdReduction;
                config.kmer_size = 20; config.kssd_half_subk = 6; config.resolution = 3;
                Sketch::API::MultiSketchBuilder builder({config});
                for (const auto& sequence : records[i]) {
                    Sketch::IO::FastxRecord record;
                    record.sequence = sequence;
                    builder.update(record);
                }
                built[i] = builder.finish(std::to_string(i));
            } else if (mode == "api") {
                Sketch::Kssd sketch(parameters);
                for (const auto& sequence : records[i]) {
#ifdef OPTIMIZED_KSSD
                    sketch.update(sequence.data(), sequence.size());
#else
                    sketch.update(sequence.c_str());
#endif
                }
                sketch.finalize(); hashes[i] = sketch.getHashes();
            }
#ifdef OPTIMIZED_KSSD
            else {
                const auto lookup = mode == "map" ? Sketch::KssdBuild::Lookup::Map
                    : mode == "direct" ? Sketch::KssdBuild::Lookup::Direct
                    : Sketch::KssdBuild::Lookup::Bitmap;
                Sketch::KssdBuild::Stream stream(*parameters, lookup);
                Sketch::KssdBuild::Accumulator<uint32_t> accumulator(mode != "bitmap-set" && mode != "file-set");
                accumulator.reset();
                auto emit = [&](uint64_t hash) { accumulator.insert(hash); };
                if (mode == "file" || mode == "file-set") {
                    Sketch::KssdBuild::FastxReader reader;
                    if (!reader.read(paths[i], stream, emit)) throw std::runtime_error("FASTX read failed");
                } else {
                    for (const auto& sequence : records[i]) {
                        stream.reset(); stream.append(sequence.data(), sequence.size(), emit);
                    }
                }
                accumulator.finish(hashes[i]);
            }
#endif
        }
        const double seconds = now() - start;
        uint64_t total = 0;
        std::ofstream dump(argv[6], std::ios::binary | std::ios::trunc);
        if (!dump) throw std::runtime_error("cannot open hash dump");
        if (mode == "builder") {
            for (const auto& left : built)
                for (const auto& right : built) {
                    const auto result = left.at(0).sketch.query(right.at(0).sketch);
                    if (!result.estimate_available) throw std::runtime_error("query unavailable");
                    dump.write(reinterpret_cast<const char*>(&result.value), sizeof(double));
                }
        } else {
            for (const auto& row : hashes) {
                const uint64_t size = row.size(); total += size;
                dump.write(reinterpret_cast<const char*>(&size), sizeof(size));
                dump.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(uint32_t));
            }
        }
        std::cout << "{\"mode\":\"" << mode << "\",\"repeat\":" << repeat
                  << ",\"threads\":" << threads << ",\"files\":" << paths.size()
                  << ",\"records\":" << record_count << ",\"bases\":" << bases
                  << ",\"hashes\":" << total << ",\"setup_seconds\":" << setup
                  << ",\"seconds\":" << seconds << "}" << std::endl;
    }
}
