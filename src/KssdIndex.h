#ifndef RABBITSKETCH_KSSD_INDEX_H
#define RABBITSKETCH_KSSD_INDEX_H

#include "api/OpenMPCompat.h"
#include "phmap.h"
#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace Sketch {
namespace KssdBuild {

struct Index {
    struct Range { size_t off; uint32_t cnt; };
    std::vector<size_t> offsets;
    std::vector<uint32_t> postings;
    phmap::flat_hash_map<uint64_t, Range> sparse;
    size_t unique_before = 0, unique_after = 0;
};

// Two-pass scatter into contiguous radix buckets replaces millions of tiny
// posting vectors and locked map merges. Keys must be unique within a sketch.
// Effective key bits (not sizeof(Key)) balance KSSD's reduced integer space.
template<typename Key, typename GetKeys>
Index buildIndex(int genomes, int threads, unsigned key_bits, GetKeys&& keys) {
    Index result;
    if (genomes <= 0) return result;
    constexpr size_t BUCKETS = 1024;
    const unsigned shift = key_bits > 10 ? key_bits - 10 : 0;
    const int workers = std::max(1, std::min(threads, genomes));
    struct Pair { Key key; uint32_t genome; };
    std::vector<size_t> counts(static_cast<size_t>(workers) * BUCKETS, 0);
    #pragma omp parallel for num_threads(workers) schedule(static)
    for (int part = 0; part < workers; ++part) {
        size_t* local = counts.data() + part * BUCKETS;
        for (int i = static_cast<int>(int64_t(genomes) * part / workers);
             i < int64_t(genomes) * (part + 1) / workers; ++i)
            for (Key key : keys(i)) ++local[static_cast<uint64_t>(key) >> shift];
    }
    std::vector<size_t> begin(BUCKETS + 1, 0);
    std::vector<size_t> cursor(counts.size(), 0);
    for (size_t bucket = 0; bucket < BUCKETS; ++bucket) {
        size_t position = begin[bucket];
        for (int worker = 0; worker < workers; ++worker) {
            const size_t index = worker * BUCKETS + bucket;
            cursor[index] = position;
            position += counts[index];
        }
        begin[bucket + 1] = position;
    }
    std::vector<Pair> pairs(begin.back());
    #pragma omp parallel for num_threads(workers) schedule(static)
    for (int part = 0; part < workers; ++part) {
        size_t* local = cursor.data() + part * BUCKETS;
        for (int i = static_cast<int>(int64_t(genomes) * part / workers);
             i < int64_t(genomes) * (part + 1) / workers; ++i)
            for (Key key : keys(i))
                pairs[local[static_cast<uint64_t>(key) >> shift]++] =
                    {key, static_cast<uint32_t>(i)};
    }
    std::vector<size_t> retained(BUCKETS, 0), unique(BUCKETS, 0), all_unique(BUCKETS, 0);
    #pragma omp parallel for num_threads(workers) schedule(dynamic, 1)
    for (size_t bucket = 0; bucket < BUCKETS; ++bucket) {
        const size_t lo = begin[bucket], hi = begin[bucket + 1];
        std::sort(pairs.begin() + lo, pairs.begin() + hi,
                  [](const Pair& a, const Pair& b) {
                      return a.key < b.key || (a.key == b.key && a.genome < b.genome);
                  });
        for (size_t i = lo; i < hi;) {
            size_t j = i + 1;
            while (j < hi && pairs[j].key == pairs[i].key) ++j;
            ++all_unique[bucket];
            if (j - i > 1) { ++unique[bucket]; retained[bucket] += j - i; }
            i = j;
        }
    }
    std::vector<size_t> posting_begin(BUCKETS + 1, 0);
    for (size_t bucket = 0; bucket < BUCKETS; ++bucket) {
        result.unique_before += all_unique[bucket];
        result.unique_after += unique[bucket];
        posting_begin[bucket + 1] = posting_begin[bucket] + retained[bucket];
    }
    result.postings.resize(posting_begin.back());
    // Avoid allocating a 2 GiB offset table for a handful of genomes. Wider
    // reduced tuples always use sparse offsets, including 32-bit drlevel=2.
    const uint64_t space = key_bits <= 28 ? UINT64_C(1) << key_bits : 0;
    const bool dense = space && space <= std::max<uint64_t>(1 << 20, pairs.size() * 8ULL);
    if (dense) result.offsets.resize(static_cast<size_t>(space) + 1, 0);
    struct Entry { uint64_t key; Index::Range range; };
    std::vector<std::vector<Entry>> entries(dense ? 0 : BUCKETS);
    #pragma omp parallel for num_threads(workers) schedule(dynamic, 1)
    for (size_t bucket = 0; bucket < BUCKETS; ++bucket) {
        const size_t hi = begin[bucket + 1];
        size_t out = posting_begin[bucket];
        if (!dense) entries[bucket].reserve(unique[bucket]);
        for (size_t i = begin[bucket]; i < hi;) {
            size_t j = i + 1;
            while (j < hi && pairs[j].key == pairs[i].key) ++j;
            if (j - i > 1) {
                if (dense) result.offsets[static_cast<size_t>(pairs[i].key) + 1] = j - i;
                else entries[bucket].push_back({static_cast<uint64_t>(pairs[i].key),
                    {out, static_cast<uint32_t>(j - i)}});
                for (size_t k = i; k < j; ++k) result.postings[out++] = pairs[k].genome;
            }
            i = j;
        }
    }
    if (dense) {
        // Parallel exclusive prefix, preserving the key-sorted posting order.
        std::vector<size_t> sums(workers + 1, 0);
        #pragma omp parallel for num_threads(workers)
        for (int part = 0; part < workers; ++part) {
            const size_t lo = 1 + space * part / workers;
            const size_t hi = 1 + space * (part + 1) / workers;
            size_t sum = 0;
            for (size_t i = lo; i < hi; ++i) {
                sum += result.offsets[i]; result.offsets[i] = sum;
            }
            sums[part + 1] = sum;
        }
        for (int part = 1; part <= workers; ++part) sums[part] += sums[part - 1];
        #pragma omp parallel for num_threads(workers)
        for (int part = 0; part < workers; ++part) {
            const size_t lo = 1 + space * part / workers;
            const size_t hi = 1 + space * (part + 1) / workers;
            for (size_t i = lo; i < hi; ++i) result.offsets[i] += sums[part];
        }
    } else {
        result.sparse.reserve(result.unique_after);
        for (const auto& bucket : entries)
            for (const auto& entry : bucket) result.sparse.emplace(entry.key, entry.range);
    }
    return result;
}

} // namespace KssdBuild
} // namespace Sketch
#endif
