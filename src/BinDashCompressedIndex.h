#ifndef RABBITSKETCH_BINDASH_COMPRESSED_INDEX_H
#define RABBITSKETCH_BINDASH_COMPRESSED_INDEX_H

#include "BinDashInvertedIndex.h"
#include "PostingCodec.h"

namespace Sketch { namespace BinDashIndex {

namespace compressed_detail = PostingCodec;

// Read-only view of the direct lookup. All arrays, including the bin metadata,
// live in the same node-local arena as the encoded posting lists.
struct NumaCompressedIndexReplica {
    struct Bin {
        size_t post_base = 0;
        Numa::Span<uint32_t> offsets, ranks;
        Numa::Span<uint16_t> values;
        Numa::Span<uint64_t> present;
        uint32_t position(uint32_t value) const {
            if (!present.empty()) {
                const uint32_t word = value >> 6;
                const uint64_t mask = UINT64_C(1) << (value & 63);
                const uint64_t bits = present[word];
                if (!(bits & mask)) return UINT32_MAX;
                return ranks[word] + static_cast<uint32_t>(__builtin_popcountll(bits & (mask - 1)));
            }
            if (values.empty()) return UINT32_MAX;
            const auto it = std::lower_bound(values.begin(), values.end(), value);
            return it == values.end() || *it != value ? UINT32_MAX : static_cast<uint32_t>(it - values.begin());
        }
    };
    std::shared_ptr<Numa::Arena> arena;
    uint32_t signature_bits;
    Numa::Span<Bin> bins;
    Numa::Span<uint8_t> encodedPosts;
    template<class Source>
    NumaCompressedIndexReplica(const Source& source, std::shared_ptr<Numa::Arena> owner)
        : arena(std::move(owner)), signature_bits(source.signature_bits) {
        bins = {static_cast<Bin*>(arena->allocate(source.bins.size() * sizeof(Bin), alignof(Bin))), source.bins.size()};
        for (size_t i = 0; i < bins.size(); ++i) {
            new (bins.data() + i) Bin();
            const auto& from = source.bins[i];
            auto& to = bins[i];
            to.post_base = from.post_base;
            to.offsets = {arena->copy(from.offsets.data(), from.offsets.size()), from.offsets.size()};
            to.values = {arena->copy(from.values.data(), from.values.size()), from.values.size()};
            to.present = {arena->copy(from.present.data(), from.present.size()), from.present.size()};
            to.ranks = {arena->copy(from.ranks.data(), from.ranks.size()), from.ranks.size()};
        }
        encodedPosts = {arena->copy(source.encodedPosts.data(), source.encodedPosts.size()), source.encodedPosts.size()};
    }
    template<class Key>
    bool findPostingRange(Key key, size_t& off, uint32_t& bytes) const {
        const size_t bin = key >> signature_bits;
        if (bin >= bins.size()) return false;
        const uint32_t value = key & ((UINT32_C(1) << signature_bits) - 1);
        const auto& table = bins[bin];
        const uint32_t p = table.position(value);
        if (p == UINT32_MAX) return false;
        off = table.post_base + table.offsets[p];
        bytes = table.offsets[p + 1] - table.offsets[p];
        return true;
    }
    template<class Visit>
    void forEachPosting(size_t off, uint32_t bytes, Visit visit) const {
        PostingCodec::forEachBytes(encodedPosts.data() + off, bytes, visit);
    }
};

// Sorted sample IDs are encoded as exact unsigned differences. Offsets refer
// to bytes, not posting counts; decoding streams directly into the caller's
// accumulator, without reconstructing a second global posting array.
struct CompressedInvertedIndex {
    using NumaReplica = NumaCompressedIndexReplica;
    uint32_t signature_bits = 0;
    std::vector<BinPostingTable> bins;
    std::vector<uint8_t> encodedPosts;
    size_t totalPostings = 0;
    size_t totalUnique = 0;
    size_t workspace_bound_bytes = 0;

    size_t numaKeyCount() const { return totalUnique; }
    size_t numaStorageBytes() const { return lookupBytes() + encodedPosts.size(); }

    template<class Key>
    bool findPostingRange(Key key, size_t& off, uint32_t& bytes) const {
        const size_t bin = key >> signature_bits;
        if (bin >= bins.size()) return false;
        const uint32_t value = key & ((UINT32_C(1) << signature_bits) - 1);
        const auto& table = bins[bin];
        const uint32_t p = table.position(value);
        if (p == UINT32_MAX) return false;
        off = table.post_base + table.offsets[p];
        bytes = table.offsets[p + 1] - table.offsets[p];
        return true;
    }

    template<class Visit>
    void forEachPosting(size_t off, uint32_t bytes, Visit visit) const {
        PostingCodec::forEachBytes(encodedPosts.data() + off, bytes, visit);
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

// Count both postings and encoded bytes in the first pass. The second pass
// writes directly into final, exactly sized storage. Its scratch allocation
// is bounded independently of N and the total number of postings.
template<class Key, class Value>
CompressedInvertedIndex buildCompressedDirectIndex(
    const std::vector<ValueKeyView<Key, Value>>& views, int n, int nThreads,
    DirectIndexBuildOptions options = {}) {
    static_assert(std::is_unsigned<Value>::value && sizeof(Value) <= 2,
                  "compressed direct indexing supports at most 16 signature bits");
    CompressedInvertedIndex index;
    if (n < 0 || size_t(n) != views.size())
        throw std::invalid_argument("BinDash index sample count does not match its input");
    if (!n) return index;
    const uint32_t bits = views[0].bits, bin_count = views[0].bins;
    if (!bits || bits > 16 || bits > sizeof(Value) * 8 || !bin_count ||
        keyBits(bin_count, bits) > sizeof(Key) * 8)
        throw std::invalid_argument("invalid compressed BinDash index dimensions");
    for (const auto& view : views)
        if (view.bits != bits || view.bins != bin_count || !view.values)
            throw std::invalid_argument("inconsistent BinDash signature views");
    // Five bytes is the worst case for one uint32_t delta. Relative byte
    // offsets deliberately retain the direct index's compact 32-bit layout.
    if (uint64_t(n) * 5 > UINT32_MAX)
        throw std::length_error("compressed BinDash bin exceeds 32-bit byte offsets");
    index.signature_bits = bits;
    index.bins.resize(bin_count);
    std::vector<size_t> logical_posts(bin_count, 0);
    const uint32_t domain = UINT32_C(1) << bits;
    const size_t row_bytes = size_t(domain) * sizeof(uint32_t) * 3;
    const size_t workspace = std::max(row_bytes, options.workspace_bytes);
    const bool sparse = uint32_t(n) < domain / 8;
    const uint32_t block_bins = sparse ? 1 : static_cast<uint32_t>(std::min<size_t>(
        std::min(bin_count, std::min(8u, std::max(1u, options.bins_per_block))), workspace / row_bytes));
    const uint32_t blocks = (uint64_t(bin_count) + block_bins - 1) / block_bins;
    const size_t scratch_bytes = row_bytes * block_bins + (sparse ? size_t(n) * sizeof(Value) : 0);
    const int workers = static_cast<int>(std::min<size_t>(
        std::min<size_t>(std::max(1, nThreads), blocks), std::max<size_t>(1, workspace / scratch_bytes)));
    index.workspace_bound_bytes = scratch_bytes * workers;

    #pragma omp parallel num_threads(workers)
    {
        const size_t cells = size_t(block_bins) * domain;
        std::vector<uint32_t> counts(cells, 0), bytes(cells, 0), last(cells, 0);
        std::vector<Value> touched;
        if (sparse) touched.reserve(n);
        #pragma omp for schedule(dynamic, 1)
        for (uint32_t block = 0; block < blocks; ++block) {
            const uint32_t first = block * block_bins;
            const uint32_t length = std::min(block_bins, bin_count - first);
            if (sparse) touched.clear();
            else {
                std::fill(counts.begin(), counts.end(), 0);
                std::fill(bytes.begin(), bytes.end(), 0);
                std::fill(last.begin(), last.end(), 0);
            }
            for (int sid = 0; sid < n; ++sid) {
                const Value* values = views[sid].values + first;
                for (uint32_t b = 0; b < length; ++b) {
                    const size_t cell = size_t(b) * domain + values[b];
                    if (!counts[cell]++ && sparse) touched.push_back(values[b]);
                    bytes[cell] += compressed_detail::encodedSize(uint32_t(sid) - last[cell]);
                    last[cell] = sid;
                }
            }
            if (sparse) std::sort(touched.begin(), touched.end());
            for (uint32_t b = 0; b < length; ++b) {
                const size_t base = size_t(b) * domain;
                auto visit_counts = [&](auto visit) {
                    auto one = [&](uint32_t value) {
                        const size_t cell = base + value;
                        visit(value, counts[cell] > 1 ? bytes[cell] : 0);
                    };
                    if (sparse) { for (Value value : touched) one(value); }
                    else { for (uint32_t value = 0; value < domain; ++value) one(value); }
                };
                direct_detail::prepareBin(index.bins[first + b], domain, options.lookup_mode, visit_counts);
                auto add = [&](uint32_t value) {
                    if (counts[base + value] > 1) logical_posts[first + b] += counts[base + value];
                };
                if (sparse) { for (Value value : touched) add(value); }
                else { for (uint32_t value = 0; value < domain; ++value) add(value); }
            }
            if (sparse) for (Value value : touched) counts[value] = bytes[value] = last[value] = 0;
        }
    }
    size_t encoded_size = 0;
    for (size_t b = 0; b < index.bins.size(); ++b) {
        auto& bin = index.bins[b];
        bin.post_base = encoded_size;
        if (bin.offsets.back() > SIZE_MAX - encoded_size)
            throw std::overflow_error("BinDash encoded index exceeds addressable memory");
        encoded_size += bin.offsets.back();
        index.totalPostings += logical_posts[b];
        index.totalUnique += bin.offsets.size() - 1;
    }
    std::vector<size_t>().swap(logical_posts);
    index.encodedPosts.resize(encoded_size);
    #pragma omp parallel num_threads(workers)
    {
        const size_t cells = size_t(block_bins) * domain;
        std::vector<uint32_t> cursor(cells, UINT32_MAX), last(cells, 0);
        #pragma omp for schedule(dynamic, 1)
        for (uint32_t block = 0; block < blocks; ++block) {
            const uint32_t first = block * block_bins;
            const uint32_t length = std::min(block_bins, bin_count - first);
            if (!sparse) std::fill(cursor.begin(), cursor.end(), UINT32_MAX);
            for (uint32_t b = 0; b < length; ++b)
                index.bins[first + b].eachValue([&](uint32_t value, uint32_t rank) {
                    const size_t cell = size_t(b) * domain + value;
                    cursor[cell] = index.bins[first + b].offsets[rank];
                    last[cell] = 0;
                });
            for (int sid = 0; sid < n; ++sid) {
                const Value* values = views[sid].values + first;
                for (uint32_t b = 0; b < length; ++b) {
                    const size_t cell = size_t(b) * domain + values[b];
                    uint32_t& pos = cursor[cell];
                    if (pos != UINT32_MAX) {
                        pos += compressed_detail::encode(index.encodedPosts.data() + index.bins[first + b].post_base + pos,
                                                         uint32_t(sid) - last[cell]);
                        last[cell] = sid;
                    }
                }
            }
            if (sparse) index.bins[first].eachValue([&](uint32_t value, uint32_t) { cursor[value] = UINT32_MAX; });
        }
    }
    std::cerr << "merge index: " << index.totalUnique << " unique keys\n"
              << "singleton removal: -> " << index.totalUnique << '\n'
              << "CSR flatten: " << index.totalPostings << " postings\n"
              << "compressed BinDash CSR: lookup=" << index.lookupBytes()
              << " bytes, postings=" << index.encodedPosts.size()
              << " bytes, temporary counters <= " << index.workspace_bound_bytes
              << " bytes, " << workers << " workers\n";
    return index;
}

template<class Key, class Value>
auto buildCompressedIndex(const std::vector<ValueKeyView<Key, Value>>& views, int n, int threads) {
    if constexpr (sizeof(Value) <= 2) return buildCompressedDirectIndex(views, n, threads);
    else return buildCSRIndexFromSortedKeys<Key>(views, n, threads);
}

}} // namespace Sketch::BinDashIndex
#endif
