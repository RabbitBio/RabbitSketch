#include "BinDashIndex.h"
#include "InvertedIndex.h"
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unistd.h>

namespace BI = Sketch::BinDashIndex;
namespace {
constexpr double infinity = std::numeric_limits<double>::infinity();
void check(bool okay, const char* message) {
    if (!okay) throw std::runtime_error(message);
}

Sketch::BinDash packed(const std::vector<uint64_t>& values, uint32_t bits) {
    const uint32_t groups = values.size() / 64;
    std::vector<uint64_t> words(size_t(groups) * bits);
    for (size_t i = 0; i < values.size(); ++i)
        for (uint32_t b = 0; b < bits; ++b)
            if ((values[i] >> b) & 1)
                words[(i / 64) * bits + b] |= UINT64_C(1) << (i % 64);
    return Sketch::BinDash::fromPacked(groups, 21, bits, 42, 0, words);
}

template<class Key>
void checkIndex(const std::vector<std::vector<uint64_t>>& values, uint32_t bits, int threads) {
    const int n = values.size();
    const uint32_t bins = values.front().size();
    std::vector<Sketch::BinDash> sketches;
    std::vector<std::vector<Key>> keys;
    std::vector<std::string> names;
    for (int i = 0; i < n; ++i) {
        sketches.push_back(packed(values[i], bits));
        for (uint32_t bin = 0; bin < bins; ++bin)
            check(sketches.back().getBinValue(bin) == values[i][bin], "bin extraction truncated bits");
        keys.push_back(BI::keys<Key>(sketches.back()));
        names.push_back(std::to_string(i));
    }
    auto index = Sketch::buildCSRIndexFromKeys<Key>(keys, n, threads, BI::Bucket{});
    std::vector<double> cutoffs = {0.0, 0.02, 0.05, 0.2, 1.0, infinity};
    for (int i = 0; i < n; ++i) {
        std::vector<uint32_t> counts(n);
        for (const auto& key : keys[i]) {
            auto found = index.postIdx.find(key);
            if (found == index.postIdx.end()) continue;
            const auto range = found->second;
            for (size_t p = range.off; p < range.off + range.cnt; ++p)
                ++counts[index.csrPosts[p]];
        }
        for (int j = i + 1; j < n; ++j) {
            uint32_t oracle = 0;
            for (uint32_t bin = 0; bin < bins; ++bin)
                oracle += values[i][bin] == values[j][bin];
            check(counts[j] == oracle, "posting count differs from position-wise oracle");
            check(BI::jaccard(counts[j], bins, bits) == sketches[i].jaccard(sketches[j]),
                  "index Jaccard differs from native packed comparison");
            const double d = sketches[i].distance(sketches[j]);
            check(BI::distance(counts[j], bins, bits, 21) == d, "index distance changed");
            if (std::isfinite(d)) {
                cutoffs.push_back(d);
                cutoffs.push_back(std::nextafter(d, infinity));
                if (d > 0.0) cutoffs.push_back(std::nextafter(d, 0.0));
            }
        }
    }
    std::sort(cutoffs.begin(), cutoffs.end());
    cutoffs.erase(std::unique(cutoffs.begin(), cutoffs.end()), cutoffs.end());
    const std::string path = "/tmp/rabbitsketch_bindash_index_test_" + std::to_string(getpid()) + ".tsv";
    const std::vector<int> sizes(n, bins);
    auto jaccard = [=](int count, int, int) { return BI::jaccard(count, bins, bits); };
    for (double cutoff : cutoffs) {
        const int minimum = static_cast<int>(BI::minimumMatches(bins, bits, 21, cutoff));
        auto minCommon = [=](int) { return minimum; };
        Sketch::computeDistances<Key, decltype(jaccard), decltype(minCommon), uint32_t>(
            index, keys, sizes, names, n, 21, cutoff, jaccard, minCommon, path, threads);
        std::vector<std::string> actual, expected;
        std::ifstream input(path);
        for (std::string line; std::getline(input, line);) actual.push_back(line);
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
            const double d = sketches[i].distance(sketches[j]);
            if (d < cutoff) {
                char row[100];
                std::snprintf(row, sizeof(row), "%d\t%d\t%.6f", i, j, d);
                expected.emplace_back(row);
            }
        }
        std::sort(actual.begin(), actual.end());
        std::sort(expected.begin(), expected.end());
        check(actual == expected, "indexed output differs from native distance cutoff oracle");
    }
    std::remove(path.c_str());
}

void cases(uint32_t bins, uint32_t bits) {
    const uint64_t mask = bits == 64 ? UINT64_MAX : (UINT64_C(1) << bits) - 1;
    std::mt19937_64 random(516);
    std::vector<std::vector<uint64_t>> values(6, std::vector<uint64_t>(bins));
    for (uint32_t i = 0; i < bins; ++i) {
        values[0][i] = random() & mask;
        values[1][i] = values[0][i];
        values[2][i] = i % 3 ? values[0][i] : random() & mask;
        values[4][i] = mask; // Native packed representation of an empty sketch.
        values[5][i] = mask;
    }
    values[3] = values[0];
    std::rotate(values[3].begin(), values[3].begin() + 1, values[3].end());
    const unsigned width = BI::keyBits(bins, bits);
    for (int threads : {1, 4}) {
        if (width <= 32) checkIndex<uint32_t>(values, bits, threads);
        else if (width <= 64) checkIndex<uint64_t>(values, bits, threads);
        else checkIndex<BI::WideKey>(values, bits, threads);
    }
}
} // namespace

int main() {
    try {
        for (uint32_t bits : {1u, 4u, 8u, 16u, 17u, 32u, 58u, 59u, 64u}) cases(64, bits);
        cases(1024, 16);
        cases(65536, 16); // Regression: default uint16_t posting counters overflow.
        std::cout << "bindash_index_test: PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "bindash_index_test: FAIL: " << e.what() << '\n';
        return 1;
    }
}
