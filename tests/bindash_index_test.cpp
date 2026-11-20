#include "BinDashIndex.h"
#include "BinDashInvertedIndex.h"
#include "BinDashCompressedIndex.h"
#include "InvertedIndex.h"
#include "SortedInvertedIndex.h"
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

template<class Key, class Value = uint64_t, bool Compressed = false>
void checkIndex(const std::vector<std::vector<uint64_t>>& values, uint32_t bits, int threads,
                BI::LookupMode mode = BI::LookupMode::Automatic) {
    const int n = values.size();
    const uint32_t bins = values.front().size();
    std::vector<Sketch::BinDash> sketches;
    std::vector<std::vector<Key>> keys;
    std::vector<std::vector<Value>> compact(n);
    std::vector<BI::ValueKeyView<Key, Value>> views(n);
    std::vector<std::string> names;
    for (int i = 0; i < n; ++i) {
        sketches.push_back(packed(values[i], bits));
        for (uint32_t bin = 0; bin < bins; ++bin)
            check(sketches.back().getBinValue(bin) == values[i][bin], "bin extraction truncated bits");
        keys.push_back(BI::keys<Key>(sketches.back()));
        compact[i].assign(values[i].begin(), values[i].end());
        views[i] = {compact[i].data(), bins, bits};
        for (uint32_t bin = 0; bin < bins; ++bin)
            check(views[i][bin] == keys.back()[bin], "compact key lost position or signature bits");
        names.push_back(std::to_string(i));
    }
    auto index = [&]() {
        if constexpr (sizeof(Value) <= 2) {
            BI::DirectIndexBuildOptions options;
            options.lookup_mode = mode;
            if constexpr (Compressed) return BI::buildCompressedDirectIndex(views, n, threads, options);
            else return BI::buildDirectIndex(views, n, threads, options);
        } else {
            return BI::buildIndex(views, n, threads);
        }
    }();
    std::vector<double> cutoffs = {0.0, 0.02, 0.05, 0.2, 1.0, infinity};
    for (int i = 0; i < n; ++i) {
        std::vector<uint32_t> counts(n);
        for (const auto& key : keys[i]) {
            size_t offset;
            uint32_t count;
            if (!index.findPostingRange(key, offset, count)) continue;
            index.forEachPosting(offset, count, [&](uint32_t sid) { ++counts[sid]; });
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
            index, views, sizes, names, n, 21, cutoff, jaccard, minCommon, path, threads, 64);
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
        if (width <= 32 && bits <= 8) checkIndex<uint32_t, uint8_t, true>(values, bits, threads);
        else if (width <= 32 && bits <= 16) checkIndex<uint32_t, uint16_t, true>(values, bits, threads);
        if (width <= 32 && bits <= 8) checkIndex<uint32_t, uint8_t>(values, bits, threads);
        else if (width <= 32 && bits <= 16) checkIndex<uint32_t, uint16_t>(values, bits, threads);
        else if (width <= 32) checkIndex<uint32_t, uint32_t>(values, bits, threads);
        else if (width <= 64 && bits <= 32) checkIndex<uint64_t, uint32_t>(values, bits, threads);
        else if (width <= 64) checkIndex<uint64_t>(values, bits, threads);
        else checkIndex<BI::WideKey>(values, bits, threads);
    }
    if (bins == 64 && bits == 16) {
        checkIndex<uint32_t, uint16_t>(values, bits, 4, BI::LookupMode::SortedValues);
        checkIndex<uint32_t, uint16_t>(values, bits, 4, BI::LookupMode::Bitmap);
    }
}

// Exercise the large-database histogram path. This independent oracle checks
// every retained and singleton value, its exact IDs, and bitmap word boundaries.
void histogramCases() {
    constexpr uint32_t bins = 64, domain = 65536, n = 100000;
    std::vector<std::vector<uint16_t>> values(n, std::vector<uint16_t>(bins));
    std::vector<BI::ValueKeyView<uint32_t, uint16_t>> views(n);
    for (uint32_t sid = 0; sid < n; ++sid) {
        for (uint32_t b = 0; b < bins; ++b) values[sid][b] = (sid + b * 123) % domain;
        views[sid] = {values[sid].data(), bins, 16};
    }
    BI::DirectIndexBuildOptions options;
    options.workspace_bytes = 2 * size_t(domain) * sizeof(uint32_t);
    auto index = BI::buildDirectIndex(views, n, 4, options);
    auto compressed_options = options;
    compressed_options.workspace_bytes = 1u << 20;
    auto compressed = BI::buildCompressedDirectIndex(views, n, 4, compressed_options);
    check(compressed.workspace_bound_bytes <= compressed_options.workspace_bytes, "compressed workspace limit exceeded");
    check(compressed.totalPostings == index.totalPostings, "compressed logical posting count changed");
    check(index.workspace_bound_bytes <= options.workspace_bytes, "histogram workspace limit exceeded");
    check(index.totalPostings == size_t(2) * (n - domain) * bins, "histogram posting count changed");
    for (uint32_t b = 0; b < bins; ++b) {
        check(!index.bins[b].present.empty(), "dense value population did not use compact bitmap lookup");
        for (uint32_t value = 0; value < domain; ++value) {
            const uint32_t first = (value + domain - b * 123) % domain;
            size_t off = 0;
            uint32_t count = 0;
            const bool found = index.findPostingRange((b << 16) | value, off, count);
            const bool repeated = first + domain < n;
            check(found == repeated, "histogram lost a repeated value or retained a singleton");
            if (found) {
                check(count == 2, "histogram range size changed");
                check(index.csrPosts[off] == first && index.csrPosts[off + 1] == first + domain,
                      "histogram posting IDs differ from independent oracle");
            }
            size_t compressed_off;
            uint32_t encoded_bytes;
            check(compressed.findPostingRange((b << 16) | value, compressed_off, encoded_bytes) == repeated,
                  "compressed index changed singleton membership");
            if (repeated) {
                std::vector<uint32_t> actual;
                compressed.forEachPosting(compressed_off, encoded_bytes, [&](uint32_t sid) { actual.push_back(sid); });
                check(actual == std::vector<uint32_t>({first, first + domain}), "compressed IDs differ from oracle");
            }
        }
    }
    // Large N with very few distinct values should keep sorted lookup arrays.
    constexpr uint32_t duplicate_n = 10000;
    views.resize(duplicate_n);
    for (uint32_t sid = 0; sid < duplicate_n; ++sid)
        for (uint32_t b = 0; b < bins; ++b) values[sid][b] = b * 13;
    auto duplicate = BI::buildDirectIndex(views, duplicate_n, 4, options);
    auto compressed_duplicate = BI::buildCompressedDirectIndex(views, duplicate_n, 4, compressed_options);
    check(compressed_duplicate.encodedPosts.size() == size_t(duplicate_n) * bins,
          "consecutive IDs did not encode in one byte each");
    check(duplicate.totalPostings == size_t(duplicate_n) * bins, "duplicate posting count changed");
    for (uint32_t b = 0; b < bins; ++b) {
        check(duplicate.bins[b].present.empty(), "sparse value population wasted a bitmap");
        size_t off;
        uint32_t count;
        check(duplicate.findPostingRange((b << 16) | (b * 13), off, count), "duplicate key missing");
        check(count == duplicate_n, "duplicate range count changed");
        for (uint32_t sid = 0; sid < duplicate_n; ++sid)
            check(duplicate.csrPosts[off + sid] == sid, "duplicate posting ID changed");
        check(compressed_duplicate.findPostingRange((b << 16) | (b * 13), off, count), "compressed duplicate key missing");
        uint32_t next = 0;
        compressed_duplicate.forEachPosting(off, count, [&](uint32_t sid) {
            check(sid == next++, "compressed consecutive ID changed");
        });
        check(next == duplicate_n, "compressed consecutive IDs lost");
    }
    duplicate.bins[0].post_base = (size_t(1) << 32) + 123;
    size_t off;
    uint32_t count;
    check(duplicate.findPostingRange(uint64_t(0), off, count) && off == (size_t(1) << 32) + 123,
          "global posting offset was truncated to 32 bits");
    views.clear();
    check(BI::buildDirectIndex(views, 0, 4).csrPosts.empty(), "empty index has postings");
    check(BI::buildCompressedDirectIndex(views, 0, 4).encodedPosts.empty(), "empty compressed index has postings");
    // Codec boundaries, including a five-byte difference, and 64-bit byte offsets.
    BI::CompressedInvertedIndex codec;
    const std::vector<uint32_t> ids = {0, 127, 255, 16639, 2297151, UINT32_MAX};
    uint32_t previous = 0;
    for (uint32_t sid : ids) {
        uint8_t bytes[5];
        const auto size = BI::compressed_detail::encode(bytes, sid - previous);
        check(size == BI::compressed_detail::encodedSize(sid - previous), "encoded size disagrees with bytes");
        codec.encodedPosts.insert(codec.encodedPosts.end(), bytes, bytes + size);
        previous = sid;
    }
    std::vector<uint32_t> decoded;
    codec.forEachPosting(0, codec.encodedPosts.size(), [&](uint32_t sid) { decoded.push_back(sid); });
    check(decoded == ids, "codec changed an integer at a byte boundary");
    compressed_duplicate.bins[0].post_base = (size_t(1) << 32) + 999;
    check(compressed_duplicate.findPostingRange(uint32_t(0), off, count) && off == (size_t(1) << 32) + 999,
          "compressed global offset was truncated");
}
} // namespace

int main() {
    try {
        for (uint32_t bits : {1u, 4u, 8u, 16u, 17u, 32u, 58u, 59u, 64u}) cases(64, bits);
        cases(1024, 16);
        cases(65536, 16); // Regression: default uint16_t posting counters overflow.
        histogramCases();
        std::cout << "bindash_index_test: PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "bindash_index_test: FAIL: " << e.what() << '\n';
        return 1;
    }
}
