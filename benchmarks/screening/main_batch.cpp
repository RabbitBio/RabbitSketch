// Synthetic batch adapter for the existing examples/main.cpp -t entry point.
// As in sub_command.cpp's list mode, one OpenMP process handles a sample list.
#include "sub_command.h"
#include <omp.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

int screeningMeasureOne(int argc, char** argv);
void screeningWarmKssd(int k, int half_subk, int drlevel);
int run_regimes_manifest(const std::string& path, int threads);

int run_screening_manifest(const Args& args) {
    try {
        if (args.threads < 1) throw std::invalid_argument("threads must be positive");
        std::ifstream input(args.screeningManifest);
        if (!input) throw std::runtime_error("cannot read screening manifest");
        std::string marker;
        std::getline(input, marker);
        if (marker == "#regimes-v1") return run_regimes_manifest(args.screeningManifest, args.threads);
        input.clear(); input.seekg(0);
        std::vector<std::vector<std::string>> rows;
        std::set<std::string> identifiers;
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream stream(line);
            std::vector<std::string> row;
            std::string field;
            while (std::getline(stream, field, '\t')) row.push_back(field);
            if ((row.size() != 12 && row.size() != 14) || row[0].size() != 24 ||
                row[0].find_first_not_of("0123456789abcdef") != std::string::npos)
                throw std::invalid_argument("invalid screening manifest row");
            if (!identifiers.insert(row[0]).second)
                throw std::invalid_argument("duplicate job in screening manifest");
            if (!rows.empty()) {
                if (row.size()!=rows.front().size()) throw std::invalid_argument("mixed input modes");
                for (size_t index : {1u, 2u, 3u, 8u, 9u, 10u, 11u})
                    if (row[index] != rows.front()[index])
                        throw std::invalid_argument("a batch must use one algorithm/k/native configuration");
            }
            rows.push_back(std::move(row));
        }
        if (rows.empty()) throw std::invalid_argument("empty screening manifest");
        const auto setup_start = std::chrono::steady_clock::now();
        if (rows[0][1] == "kssd")
            screeningWarmKssd(std::stoi(rows[0][3]), std::stoi(rows[0][10]), std::stoi(rows[0][2]));
        const double setup = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - setup_start).count();
        omp_set_dynamic(0);
        int failed = 0;
        const char* keys[] = {"--job-id", "--algorithm", "--parameter", "--k",
            "--n-a", "--n-b", "--intersection", "--data-seed", "--hash-seed",
            "--query-min-ms", "--half-subk", "--bbits", "--fasta-a", "--fasta-b"};
        #pragma omp parallel num_threads(args.threads) reduction(+:failed)
        {
            #pragma omp single
            std::cout << "{\"batch_event\":\"started\",\"threads\":" << omp_get_num_threads()
                      << ",\"pairs\":" << rows.size() << ",\"shared_setup_seconds\":" << setup
                      << "}\n" << std::flush;
            #pragma omp for schedule(dynamic, 1)
            for (size_t index = 0; index < rows.size(); ++index) {
                std::vector<std::string> values{"rabbitsketch-screening"};
                for (size_t field = 0; field < rows[index].size(); ++field) {
                    values.emplace_back(keys[field]); values.push_back(rows[index][field]);
                }
                values.emplace_back("--batch-threads"); values.push_back(std::to_string(omp_get_num_threads()));
                values.emplace_back("--worker-thread"); values.push_back(std::to_string(omp_get_thread_num()));
                std::vector<char*> argv;
                for (auto& value : values) argv.push_back(value.data());
                failed += screeningMeasureOne(static_cast<int>(argv.size()), argv.data()) != 0;
            }
        }
        return failed ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "screening batch: " << error.what() << '\n';
        return 2;
    }
}
