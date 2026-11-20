#ifndef RABBITSKETCH_PARTITIONED_INDEX_H
#define RABBITSKETCH_PARTITIONED_INDEX_H

#include "phmap.h"
#include "PostingCodec.h"
#include "NumaIndex.h"
#include "api/OpenMPCompat.h"
#include <algorithm>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Sketch {

template<class Key>
struct IndexKeyView {
    const Key* values = nullptr;
    size_t count = 0;
    size_t size() const { return count; }
    const Key& operator[](size_t i) const { return values[i]; }
};

template<class Key>
struct InvertedIndex {
    // off is an encoded byte offset; cnt is a logical occurrence count.
    struct PostRange { size_t off; uint32_t cnt; };
    using NumaReplica = NumaHashIndexReplica<Key, PostRange>;
    phmap::flat_hash_map<Key, PostRange> postIdx;
    std::vector<uint8_t> encodedPosts;
    size_t totalPostings = 0;
    size_t numaKeyCount() const { return postIdx.size(); }
    size_t numaStorageBytes() const {
        return encodedPosts.size() + postIdx.capacity() * (sizeof(typename decltype(postIdx)::value_type) + 1) + 64;
    }

    bool findPostingRange(const Key& key, size_t& off, uint32_t& count) const {
        const auto it = postIdx.find(key);
        if (it == postIdx.end()) return false;
        off = it->second.off;
        count = it->second.cnt;
        return true;
    }
    template<class Visit>
    void forEachPosting(size_t off, uint32_t count, Visit visit) const {
        if (count) PostingCodec::forEachCount(encodedPosts.data() + off, count, visit);
    }
};

struct IndexBuildOptions {
    size_t max_buckets = 1024;
    // Hard limit on concurrent temporary pair/merge buffers. Input sketches,
    // partition metadata, and the final compressed index are separate storage.
    size_t workspace_bytes = size_t(8) << 30;
};
struct IndexBuildStats {
    size_t input_entries = 0, buckets = 0, largest_bucket_entries = 0;
    size_t sort_workspace_bound_bytes = 0, partition_metadata_bytes = 0;
    size_t streaming_partitions = 0, unique_before = 0, unique_after = 0;
    int sort_workers = 0;
};

// Common engine for native sorted views, mutable key arrays, and read-only
// unsorted arrays. No complete (key, sample) copy or decoded posting copy is
// retained. Equal keys share a partition; occurrence multiplicities survive.
// Mutable callers sort their owned arrays first. Read-only unsorted callers
// use bounded rescans instead of making a second full key array.
template<class Key, class KeyRanges>
InvertedIndex<Key> buildPartitionedIndex(
    const KeyRanges& keys, int n, int threads, bool sorted,
    IndexBuildOptions options = {}, IndexBuildStats* stats = nullptr) {
    InvertedIndex<Key> result;
    if (stats) *stats = {};
    if (n < 0 || static_cast<size_t>(n) > keys.size() || !options.max_buckets ||
        !options.workspace_bytes)
        throw std::invalid_argument("invalid partitioned index input or workspace");
    if (!n) return result;
    size_t total = 0;
    for (int i = 0; i < n; ++i) {
        if (keys[i].size() > UINT32_MAX || keys[i].size() > SIZE_MAX - total)
            throw std::length_error("index input is too large");
        total += keys[i].size();
    }
    if (stats) stats->input_entries = total;
    if (!total) return result;
    const int workers = std::max(1, std::min(threads, n));
    const size_t wanted = std::min(options.max_buckets, std::max<size_t>(1, total / 65536));
    const size_t sample_count = std::min<size_t>(65536, total);
    std::vector<Key> sample;
    sample.reserve(sample_count);
    int genome = 0;
    size_t prefix = 0;
    for (size_t s = 0; s < sample_count; ++s) {
        const size_t lo = (total / sample_count) * s + (total % sample_count) * s / sample_count;
        const size_t hi = (total / sample_count) * (s + 1) +
                          (total % sample_count) * (s + 1) / sample_count;
        uint64_t random = s + UINT64_C(0x9e3779b97f4a7c15);
        random = (random ^ (random >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        random = (random ^ (random >> 27)) * UINT64_C(0x94d049bb133111eb);
        const size_t position = lo + (random ^ (random >> 31)) % (hi - lo);
        while (position - prefix >= keys[genome].size()) {
            prefix += keys[genome].size();
            ++genome;
        }
        sample.push_back(keys[genome][position - prefix]);
    }
    std::sort(sample.begin(), sample.end());
    std::vector<Key> cuts;
    for (size_t b = 1; b < wanted; ++b) cuts.push_back(sample[sample_count * b / wanted]);
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    std::vector<Key>().swap(sample);
    const size_t buckets = cuts.size() + 1;
    if (static_cast<size_t>(n) > SIZE_MAX / (buckets + 1) / sizeof(uint32_t))
        throw std::length_error("index partition offsets overflow");
    std::vector<uint32_t> bounds(sorted ? static_cast<size_t>(n) * (buckets + 1) : 0);
    std::vector<size_t> lengths(buckets, 0);
    if (sorted) {
        #pragma omp parallel for num_threads(workers) schedule(static)
        for (int i = 0; i < n; ++i) {
            uint32_t* row = bounds.data() + static_cast<size_t>(i) * (buckets + 1);
            row[0] = 0;
            for (size_t b = 0; b < cuts.size(); ++b) {
                size_t lo = row[b], hi = keys[i].size();
                while (lo < hi) {
                    const size_t mid = lo + (hi - lo) / 2;
                    if (keys[i][mid] < cuts[b]) lo = mid + 1;
                    else hi = mid;
                }
                row[b + 1] = static_cast<uint32_t>(lo);
            }
            row[buckets] = static_cast<uint32_t>(keys[i].size());
        }
        for (int i = 0; i < n; ++i) {
            const uint32_t* row = bounds.data() + static_cast<size_t>(i) * (buckets + 1);
            for (size_t b = 0; b < buckets; ++b) lengths[b] += row[b + 1] - row[b];
        }
    } else {
        for (int i = 0; i < n; ++i) for (size_t k = 0; k < keys[i].size(); ++k)
            ++lengths[std::upper_bound(cuts.begin(), cuts.end(), Key(keys[i][k])) - cuts.begin()];
    }
    struct Pair { Key key; uint32_t genome; };
    struct Cursor { Key key; uint32_t genome, offset; };
    const auto less = [](const auto& a, const auto& b) {
        return a.key < b.key || (!(b.key < a.key) && a.genome < b.genome);
    };
    const size_t largest = *std::max_element(lengths.begin(), lengths.end());
    const size_t heap_bytes = static_cast<size_t>(n) <= SIZE_MAX / sizeof(Cursor)
        ? static_cast<size_t>(n) * sizeof(Cursor) : SIZE_MAX;
    std::vector<size_t> scratch(buckets, 0);
    size_t largest_scratch = 0, streaming = 0;
    for (size_t b = 0; b < buckets; ++b) {
        if (lengths[b] <= options.workspace_bytes / sizeof(Pair))
            scratch[b] = lengths[b] * sizeof(Pair);
        else {
            ++streaming;
            // A skewed partition must not break the workspace limit. Merge
            // sorted runs with one cursor per sample, or rescan with O(1)
            // scratch when even that heap cannot fit the requested budget.
            if (sorted && heap_bytes <= options.workspace_bytes) scratch[b] = heap_bytes;
        }
        largest_scratch = std::max(largest_scratch, scratch[b]);
    }
    const int sort_workers = static_cast<int>(std::min<size_t>(
        std::min<size_t>(workers, buckets), largest_scratch
        ? std::max<size_t>(1, options.workspace_bytes / largest_scratch) : workers));
    if (stats) {
        stats->buckets = buckets;
        stats->largest_bucket_entries = largest;
        stats->sort_workers = sort_workers;
        stats->sort_workspace_bound_bytes = largest_scratch * sort_workers;
        stats->partition_metadata_bytes = bounds.size() * sizeof(uint32_t);
        stats->streaming_partitions = streaming;
    }
    auto each = [&](size_t b, auto visit) {
        for (int i = 0; i < n; ++i) {
            const size_t lo = sorted ? bounds[size_t(i) * (buckets + 1) + b] : 0;
            const size_t hi = sorted ? bounds[size_t(i) * (buckets + 1) + b + 1] : keys[i].size();
            for (size_t k = lo; k < hi; ++k) {
                const Key key = keys[i][k];
                if (sorted || ((b == 0 || !(key < cuts[b - 1])) &&
                               (b == cuts.size() || key < cuts[b]))) visit(key, uint32_t(i));
            }
        }
    };
    struct Partition {
        std::vector<Key> keys;
        std::vector<uint32_t> counts;
        std::vector<size_t> bytes;
        std::vector<uint8_t> posts;
        size_t unique_before = 0, logical_posts = 0;
    };
    std::vector<Partition> partitions(buckets);
    std::exception_ptr error;
    #pragma omp parallel for num_threads(sort_workers) schedule(dynamic, 1)
    for (size_t b = 0; b < buckets; ++b) {
        try {
            if (!lengths[b]) continue;
            std::vector<Pair> pairs;
            const bool pair_sort = lengths[b] <= options.workspace_bytes / sizeof(Pair);
            if (pair_sort) {
                pairs.reserve(lengths[b]);
                each(b, [&](const Key& key, uint32_t sid) { pairs.push_back({key, sid}); });
                std::sort(pairs.begin(), pairs.end(), less);
            }
            auto traverse = [&](auto visit) {
                if (pair_sort) {
                    for (const auto& pair : pairs) visit(pair.key, pair.genome);
                } else if (sorted && heap_bytes <= options.workspace_bytes) {
                    std::vector<Cursor> heap;
                    heap.reserve(n);
                    const auto greater = [&](const Cursor& a, const Cursor& c) { return less(c, a); };
                    for (int i = 0; i < n; ++i) {
                        const uint32_t* row = bounds.data() + size_t(i) * (buckets + 1);
                        if (row[b] != row[b + 1]) heap.push_back({keys[i][row[b]], uint32_t(i), row[b]});
                    }
                    std::make_heap(heap.begin(), heap.end(), greater);
                    while (!heap.empty()) {
                        std::pop_heap(heap.begin(), heap.end(), greater);
                        Cursor cursor = heap.back();
                        heap.pop_back();
                        visit(cursor.key, cursor.genome);
                        if (++cursor.offset != bounds[size_t(cursor.genome) * (buckets + 1) + b + 1]) {
                            cursor.key = keys[cursor.genome][cursor.offset];
                            heap.push_back(cursor);
                            std::push_heap(heap.begin(), heap.end(), greater);
                        }
                    }
                } else {
                    Key previous{}, next{};
                    bool have_previous = false;
                    for (;;) {
                        bool found = false;
                        each(b, [&](const Key& key, uint32_t) {
                            if ((!have_previous || previous < key) && (!found || key < next)) {
                                next = key; found = true;
                            }
                        });
                        if (!found) break;
                        each(b, [&](const Key& key, uint32_t sid) { if (key == next) visit(key, sid); });
                        previous = next; have_previous = true;
                    }
                }
            };
            auto& out = partitions[b];
            size_t unique = 0, encoded = 0;
            // Both traversals are sorted by (key, sample), including duplicate
            // occurrences. Count exactly before allocating the encoded output.
            for (bool writing : {false, true}) {
                if (writing) {
                    out.keys.reserve(unique); out.counts.reserve(unique); out.bytes.reserve(unique);
                    out.posts.resize(encoded);
                }
                bool have_key = false;
                Key key{};
                uint64_t count = 0;
                uint32_t first = 0, last = 0;
                size_t cursor = 0, group_start = 0;
                auto append = [&](uint32_t delta) {
                    if (writing) cursor += PostingCodec::encode(out.posts.data() + cursor, delta);
                    else cursor += PostingCodec::encodedSize(delta);
                };
                auto finish = [&] {
                    if (!have_key) return;
                    if (!writing) ++out.unique_before;
                    if (count <= 1) return;
                    if (writing) {
                        out.keys.push_back(key); out.counts.push_back(static_cast<uint32_t>(count));
                        out.bytes.push_back(cursor - group_start); out.logical_posts += count;
                    } else ++unique;
                };
                traverse([&](const Key& current, uint32_t sid) {
                    if (!have_key || !(current == key)) {
                        finish(); key = current; have_key = true; count = 0; group_start = cursor;
                    }
                    if (++count > UINT32_MAX) throw std::length_error("one index posting list exceeds UINT32_MAX occurrences");
                    if (count == 1) first = sid;
                    else {
                        if (count == 2) append(first);
                        append(sid - last);
                    }
                    last = sid;
                });
                finish();
                if (!writing) encoded = cursor;
            }
        } catch (...) {
            #pragma omp critical(rabbitsketch_partition_build_error)
            { if (!error) error = std::current_exception(); }
        }
    }
    if (error) std::rethrow_exception(error);
    std::vector<uint32_t>().swap(bounds);
    size_t total_unique = 0, encoded_size = 0, unique_before = 0;
    for (const auto& part : partitions) {
        if (part.posts.size() > SIZE_MAX - encoded_size) throw std::length_error("compressed index too large");
        encoded_size += part.posts.size(); total_unique += part.keys.size();
        result.totalPostings += part.logical_posts; unique_before += part.unique_before;
    }
    result.encodedPosts.reserve(encoded_size);
    result.postIdx.reserve(total_unique);
    for (auto& part : partitions) {
        size_t off = result.encodedPosts.size();
        result.encodedPosts.insert(result.encodedPosts.end(), part.posts.begin(), part.posts.end());
        std::vector<uint8_t>().swap(part.posts);
        for (size_t i = 0; i < part.keys.size(); ++i) {
            result.postIdx.emplace(part.keys[i], typename InvertedIndex<Key>::PostRange{off, part.counts[i]});
            off += part.bytes[i];
        }
        std::vector<Key>().swap(part.keys);
        std::vector<uint32_t>().swap(part.counts);
        std::vector<size_t>().swap(part.bytes);
    }
    if (stats) { stats->unique_before = unique_before; stats->unique_after = total_unique; }
    std::cerr << "merge index: " << total_unique << " unique keys\n"
              << "singleton removal: -> " << total_unique << '\n'
              << "CSR flatten: " << result.totalPostings << " postings\n"
              << "compressed CSR: " << result.encodedPosts.size() << " posting bytes, "
              << buckets << " partitions, " << sort_workers << " workers, temporary sorting <= "
              << largest_scratch * sort_workers << " bytes, " << streaming << " streaming partitions\n";
    return result;
}

} // namespace Sketch
#endif
