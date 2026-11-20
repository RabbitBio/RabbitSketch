#ifndef RABBITSKETCH_SORTED_INVERTED_INDEX_H
#define RABBITSKETCH_SORTED_INVERTED_INDEX_H

#include "InvertedIndex.h"

namespace Sketch {
using SortedIndexBuildOptions = IndexBuildOptions;
using SortedIndexBuildStats = IndexBuildStats;

// Sorted native sketch views avoid copying or modifying the sketch payload.
template<class Key, class KeyRanges>
InvertedIndex<Key> buildCSRIndexFromSortedKeys(
    const KeyRanges& keys, int n, int threads,
    SortedIndexBuildOptions options = {}, SortedIndexBuildStats* stats = nullptr) {
    return buildPartitionedIndex<Key>(keys, n, threads, true, options, stats);
}
} // namespace Sketch
#endif
