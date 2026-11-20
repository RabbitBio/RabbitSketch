#ifndef RABBITSKETCH_KSSD_INDEX_H
#define RABBITSKETCH_KSSD_INDEX_H

#include "SortedInvertedIndex.h"

namespace Sketch { namespace KssdBuild {

struct Index : InvertedIndex<uint64_t> {
    size_t unique_before = 0, unique_after = 0;
};

// KSSD's finalized keys are already sorted. View their native width without
// copying; share partition bounds, compressed postings and the query visitor
// with the other sketch families. The lookup retains exact 64-bit keys.
template<class Key, class GetKeys>
Index buildIndex(int genomes, int threads, unsigned key_bits, GetKeys&& keys,
                 IndexBuildOptions options = {}, IndexBuildStats* stats = nullptr) {
    if (genomes < 0 || key_bits > sizeof(Key) * 8)
        throw std::invalid_argument("invalid KSSD index dimensions");
    std::vector<IndexKeyView<Key>> views;
    views.reserve(genomes);
    for (int i = 0; i < genomes; ++i) {
        const auto& values = keys(i);
        views.push_back({values.data(), values.size()});
    }
    IndexBuildStats measured;
    Index result;
    static_cast<InvertedIndex<uint64_t>&>(result) = buildCSRIndexFromSortedKeys<uint64_t>(
        views, genomes, threads, options, &measured);
    result.unique_before = measured.unique_before;
    result.unique_after = measured.unique_after;
    if (stats) *stats = measured;
    return result;
}

}} // namespace Sketch::KssdBuild
#endif
