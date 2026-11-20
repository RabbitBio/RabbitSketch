#ifndef RABBITSKETCH_BINDASH_INDEX_H
#define RABBITSKETCH_BINDASH_INDEX_H

#include "BinDash.h"
#include "estimators/SimilarityThreshold.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace Sketch { namespace BinDashIndex {

inline uint64_t mix(uint64_t x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

// Preserve every signature bit even when bin position + value exceed 64 bits.
// Hash collisions select map buckets; equality compares the complete key.
struct WideKey {
    uint64_t value;
    uint32_t bin;
    bool operator==(const WideKey& other) const {
        return bin == other.bin && value == other.value;
    }
    bool operator<(const WideKey& other) const {
        return bin < other.bin || (bin == other.bin && value < other.value);
    }
    friend size_t hash_value(const WideKey& key) {
        return static_cast<size_t>(mix(key.value) ^ mix(key.bin));
    }
};

// Bin positions are implicit in the array index. Store only the signature
// value and synthesize the exact composite key when it is used. Values may
// use 8/16/32/64 bits; no signature bits are discarded.
template<class Key, class Value>
struct ValueKeyView {
    const Value* values = nullptr;
    uint32_t bins = 0, bits = 0;
    size_t size() const { return bins; }
    Key operator[](size_t bin) const {
        if constexpr (std::is_same<Key, WideKey>::value)
            return WideKey{static_cast<uint64_t>(values[bin]), static_cast<uint32_t>(bin)};
        else
            return (static_cast<Key>(bin) << bits) | static_cast<Key>(values[bin]);
    }
};

inline unsigned keyBits(uint32_t bins, uint32_t bbits) {
    unsigned position_bits = 0;
    for (uint32_t last = bins - 1; last; last >>= 1) ++position_bits;
    return bbits + position_bits;
}

template<class Key>
std::vector<Key> keys(const BinDash& sketch) {
    const uint32_t bins = sketch.getNumBins(), bits = sketch.getBbits();
    if constexpr (!std::is_same<Key, WideKey>::value) {
        static_assert(std::is_unsigned<Key>::value, "index keys must be unsigned");
        if (keyBits(bins, bits) > sizeof(Key) * 8)
            throw std::invalid_argument("BinDash index key would truncate the bin or signature");
    }
    sketch.packedSignatures();
    std::vector<Key> result(bins);
    for (uint32_t bin = 0; bin < bins; ++bin) {
        const uint64_t value = sketch.getBinValue(bin);
        if constexpr (std::is_same<Key, WideKey>::value)
            result[bin] = WideKey{value, bin};
        else
            result[bin] = (static_cast<Key>(bin) << bits) | static_cast<Key>(value);
    }
    return result;
}

struct Bucket {
    size_t operator()(uint64_t key) const { return static_cast<size_t>(mix(key) >> 54); }
    size_t operator()(const WideKey& key) const {
        return static_cast<size_t>((mix(key.value) ^ mix(key.bin)) >> 54);
    }
};

inline double jaccard(uint32_t common, uint32_t bins, uint32_t bbits) {
    const double random = std::ldexp(1.0, -static_cast<int>(bbits));
    const double matched = static_cast<double>(common) / bins;
    return std::max(0.0, std::min(1.0, (matched - random) / (1.0 - random)));
}

inline double distance(uint32_t common, uint32_t bins, uint32_t bbits, int kmer) {
    const double j = jaccard(common, bins, bbits);
    if (j <= 0.0) return std::numeric_limits<double>::infinity();
    if (j >= 1.0) return 0.0;
    return -std::log(2.0 * j / (1.0 + j)) / static_cast<double>(kmer);
}

// Compute once per search using the native forward transform. Inverse + ceil
// can discard a qualifying pair at a floating-point boundary. bins+1 means
// even an identical pair cannot pass, e.g. strict distance < 0.
inline uint64_t minimumMatches(uint32_t bins, uint32_t bbits, int kmer, double cutoff) {
    detail::validateDistanceThreshold(cutoff);
    uint64_t lo = 0, hi = static_cast<uint64_t>(bins) + 1;
    while (lo < hi) {
        const uint64_t mid = lo + (hi - lo) / 2;
        if (distance(static_cast<uint32_t>(mid), bins, bbits, kmer) < cutoff)
            hi = mid;
        else
            lo = mid + 1;
    }
    return lo;
}

}} // namespace Sketch::BinDashIndex
#endif
