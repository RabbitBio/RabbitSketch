#ifndef RABBITSKETCH_BINDASH_INVERTED_INDEX_H
#define RABBITSKETCH_BINDASH_INVERTED_INDEX_H

#include "BinDashIndex.h"
#include "SortedInvertedIndex.h"
#include <climits>
#include <cstdint>
#include <limits>
#include <vector>

namespace Sketch { namespace BinDashIndex {

enum class LookupMode { Automatic, SortedValues, Bitmap };

struct DirectIndexBuildOptions {
    size_t workspace_bytes = size_t(1) << 30;
    uint32_t bins_per_block = 32;
    LookupMode lookup_mode = LookupMode::Automatic;
};

// A genome has exactly one value per bin, so each bin has at most N postings.
// Relative offsets therefore fit in uint32_t even when the whole index exceeds
// 2^32 postings. The bin's global base remains size_t.
struct BinPostingTable {
    size_t post_base = 0;
    std::vector<uint32_t> offsets;
    std::vector<uint16_t> values;
    std::vector<uint64_t> present;
    std::vector<uint32_t> ranks;

    uint32_t position(uint32_t value) const {
        if (!present.empty()) {
            const uint32_t word = value >> 6;
            const uint64_t mask = UINT64_C(1) << (value & 63);
            const uint64_t bits = present[word];
            if (!(bits & mask)) return UINT32_MAX;
            return ranks[word] + static_cast<uint32_t>(__builtin_popcountll(bits & (mask - 1)));
        }
        const auto it = std::lower_bound(values.begin(), values.end(), value);
        if (it == values.end() || *it != value) return UINT32_MAX;
        return static_cast<uint32_t>(it - values.begin());
    }

    template<class Fn>
    void eachValue(Fn visit) const {
        if (present.empty()) {
            for (uint32_t i = 0; i < values.size(); ++i) visit(values[i], i);
        } else {
            uint32_t rank = 0;
            for (uint32_t w = 0; w < present.size(); ++w) {
                uint64_t bits = present[w];
                while (bits) {
                    const uint32_t value = (w << 6) + __builtin_ctzll(bits);
                    visit(value, rank++);
                    bits &= bits - 1;
                }
            }
        }
    }
};

struct DirectInvertedIndex {
    uint32_t signature_bits = 0;
    std::vector<BinPostingTable> bins;
    std::vector<uint32_t> csrPosts;
    size_t totalPostings = 0;
    size_t totalUnique = 0;
    size_t workspace_bound_bytes = 0;

    template<class Visit>
    void forEachPosting(size_t off, uint32_t count, Visit visit) const {
        const uint32_t* posts = csrPosts.data() + off;
        for (uint32_t p = 0; p < count; ++p) visit(posts[p]);
    }

    template<class Key>
    bool findPostingRange(Key key, size_t& off, uint32_t& count) const {
        static_assert(std::is_unsigned<Key>::value, "direct BinDash keys must be unsigned");
        const size_t bin = key >> signature_bits;
        if (bin >= bins.size()) return false;
        const uint32_t value = key & ((UINT32_C(1) << signature_bits) - 1);
        const BinPostingTable& table = bins[bin];
        const uint32_t p = table.position(value);
        if (p == UINT32_MAX) return false;
        off = table.post_base + table.offsets[p];
        count = table.offsets[p + 1] - table.offsets[p];
        return true;
    }

    size_t lookupBytes() const {
        size_t bytes = bins.capacity() * sizeof(BinPostingTable);
        for (const auto& bin : bins)
            bytes += bin.offsets.capacity() * sizeof(uint32_t)
                   + bin.values.capacity() * sizeof(uint16_t)
                   + bin.present.capacity() * sizeof(uint64_t)
                   + bin.ranks.capacity() * sizeof(uint32_t);
        return bytes;
    }
};

namespace direct_detail {

// enumerate visits values in ascending order. Count first to allocate exact
// lookup sizes: no geometric growth or second copy of all posting lists.
template<class Enumerate>
void prepareBin(BinPostingTable& table, uint32_t domain, LookupMode mode, Enumerate enumerate) {
    uint32_t unique = 0;
    enumerate([&](uint32_t, uint32_t count) { unique += count > 1; });
    table.offsets.resize(size_t(unique) + 1);
    const size_t words = (size_t(domain) + 63) / 64;
    const bool bitmap = unique && (mode == LookupMode::Bitmap ||
        (mode == LookupMode::Automatic && words * (sizeof(uint64_t) + sizeof(uint32_t))
            < size_t(unique) * sizeof(uint16_t)));
    if (bitmap) {
        table.present.resize(words, 0);
        table.ranks.resize(words);
    } else {
        table.values.resize(unique);
    }
    uint32_t p = 0, offset = 0;
    enumerate([&](uint32_t value, uint32_t count) {
        if (count < 2) return; // Singletons cannot contribute to any pair.
        table.offsets[p] = offset;
        offset += count;
        if (bitmap) table.present[value >> 6] |= UINT64_C(1) << (value & 63);
        else table.values[p] = static_cast<uint16_t>(value);
        ++p;
    });
    table.offsets[p] = offset;
    uint32_t rank = 0;
    for (size_t w = 0; w < table.present.size(); ++w) {
        table.ranks[w] = rank;
        rank += __builtin_popcountll(table.present[w]);
    }
}

} // namespace direct_detail

// Two input passes, with bounded per-worker counters. Pass one counts values
// within each bin and constructs compact lookup tables. Pass two writes genome
// IDs directly into the final CSR array in ascending order. No global (key, ID)
// array, hash map, or intermediate collection of posting lists is created.
template<class Key, class Value>
DirectInvertedIndex buildDirectIndex(const std::vector<ValueKeyView<Key, Value>>& views,
                                    int n, int nThreads,
                                    DirectIndexBuildOptions options = {}) {
    static_assert(std::is_unsigned<Value>::value && sizeof(Value) <= 2,
                  "direct BinDash indexing supports signatures up to 16 bits");
    DirectInvertedIndex index;
    if (n < 0 || size_t(n) != views.size())
        throw std::invalid_argument("BinDash index sample count does not match its input");
    if (!n) return index;
    const uint32_t bits = views[0].bits, bin_count = views[0].bins;
    if (!bits || bits > 16 || bits > sizeof(Value) * 8 || !bin_count ||
        keyBits(bin_count, bits) > sizeof(Key) * 8)
        throw std::invalid_argument("invalid direct BinDash index dimensions");
    for (const auto& view : views)
        if (view.bits != bits || view.bins != bin_count || !view.values)
            throw std::invalid_argument("inconsistent BinDash signature views");
    index.signature_bits = bits;
    index.bins.resize(bin_count);
    const uint32_t domain = UINT32_C(1) << bits;
    const size_t row_bytes = size_t(domain) * sizeof(uint32_t);
    const size_t workspace = std::max(row_bytes, options.workspace_bytes);
    // Sparse samples enumerate only touched values; scanning all 65536 values
    // for every bin would be wasteful for small test sets and small databases.
    const bool sparse_input = uint32_t(n) < domain / 8;
    const uint32_t block_bins = sparse_input ? 1 : static_cast<uint32_t>(std::min<size_t>(
        std::min(bin_count, std::max(1u, options.bins_per_block)), workspace / row_bytes));
    const uint32_t blocks = (uint64_t(bin_count) + block_bins - 1) / block_bins;
    const size_t scratch_bytes = row_bytes * block_bins
        + (sparse_input ? size_t(n) * sizeof(Value) : 0);
    const int workers = static_cast<int>(std::min<size_t>(
        std::min<size_t>(std::max(1, nThreads), blocks), std::max<size_t>(1, workspace / scratch_bytes)));
    index.workspace_bound_bytes = scratch_bytes * workers;

    #pragma omp parallel num_threads(workers)
    {
        std::vector<uint32_t> counts(size_t(block_bins) * domain, 0);
        std::vector<Value> touched;
        if (sparse_input) touched.reserve(n);
        #pragma omp for schedule(dynamic, 1)
        for (uint32_t block = 0; block < blocks; ++block) {
            const uint32_t first = block * block_bins;
            const uint32_t length = std::min(block_bins, bin_count - first);
            if (sparse_input) {
                touched.clear();
                for (int sid = 0; sid < n; ++sid) {
                    const Value value = views[sid].values[first];
                    if (!counts[value]++) touched.push_back(value);
                }
                std::sort(touched.begin(), touched.end());
                direct_detail::prepareBin(index.bins[first], domain, options.lookup_mode,
                    [&](auto visit) { for (Value v : touched) visit(v, counts[v]); });
                for (Value value : touched) counts[value] = 0;
            } else {
                std::fill(counts.begin(), counts.end(), 0);
                for (int sid = 0; sid < n; ++sid) {
                    const Value* values = views[sid].values + first;
                    for (uint32_t b = 0; b < length; ++b)
                        ++counts[size_t(b) * domain + values[b]];
                }
                for (uint32_t b = 0; b < length; ++b) {
                    const uint32_t* row = counts.data() + size_t(b) * domain;
                    direct_detail::prepareBin(index.bins[first + b], domain, options.lookup_mode,
                        [&](auto visit) { for (uint32_t v = 0; v < domain; ++v) visit(v, row[v]); });
                }
            }
        }
    }
    for (auto& bin : index.bins) {
        bin.post_base = index.totalPostings;
        if (bin.offsets.back() > std::numeric_limits<size_t>::max() - index.totalPostings)
            throw std::overflow_error("BinDash posting count exceeds addressable memory");
        index.totalPostings += bin.offsets.back();
        index.totalUnique += bin.offsets.size() - 1;
    }
    index.csrPosts.resize(index.totalPostings);
    #pragma omp parallel num_threads(workers)
    {
        std::vector<uint32_t> cursor(size_t(block_bins) * domain, UINT32_MAX);
        #pragma omp for schedule(dynamic, 1)
        for (uint32_t block = 0; block < blocks; ++block) {
            const uint32_t first = block * block_bins;
            const uint32_t length = std::min(block_bins, bin_count - first);
            if (!sparse_input) std::fill(cursor.begin(), cursor.end(), UINT32_MAX);
            for (uint32_t b = 0; b < length; ++b) {
                const auto& bin = index.bins[first + b];
                bin.eachValue([&](uint32_t value, uint32_t rank) {
                    cursor[size_t(b) * domain + value] = bin.offsets[rank];
                });
            }
            for (int sid = 0; sid < n; ++sid) {
                const Value* values = views[sid].values + first;
                for (uint32_t b = 0; b < length; ++b) {
                    uint32_t& pos = cursor[size_t(b) * domain + values[b]];
                    if (pos != UINT32_MAX)
                        index.csrPosts[index.bins[first + b].post_base + pos++] = sid;
                }
            }
            if (sparse_input) index.bins[first].eachValue([&](uint32_t value, uint32_t) {
                cursor[value] = UINT32_MAX;
            });
        }
    }
    std::cerr << "merge index: " << index.totalUnique << " unique keys\n"
              << "singleton removal: -> " << index.totalUnique << '\n'
              << "CSR flatten: " << index.totalPostings << " postings\n"
              << "direct BinDash CSR: lookup=" << index.lookupBytes()
              << " bytes, postings=" << index.csrPosts.size() * sizeof(uint32_t)
              << " bytes, temporary counters <= " << index.workspace_bound_bytes
              << " bytes, " << workers << " workers\n";
    return index;
}

// Wider signatures retain all their bits through the generic sorted builder.
template<class Key, class Value>
auto buildIndex(const std::vector<ValueKeyView<Key, Value>>& views, int n, int nThreads) {
    if constexpr (sizeof(Value) <= 2) return buildDirectIndex(views, n, nThreads);
    else return buildCSRIndexFromSortedKeys<Key>(views, n, nThreads);
}

}} // namespace Sketch::BinDashIndex
#endif
