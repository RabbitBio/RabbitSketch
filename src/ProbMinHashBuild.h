#ifndef RABBITSKETCH_PROBMINHASH_BUILD_H
#define RABBITSKETCH_PROBMINHASH_BUILD_H

#include "api/RuntimeInfo.h"
#include "BoundedRadixSort.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Sketch {
namespace ProbMinHashBuild {

// Owned by one worker, reused across input files, explicitly released after
// construction. No counts or k-mer identities survive between samples.
struct Workspace {
    std::vector<uint64_t> occurrences;
    std::unique_ptr<uint64_t[]> scratch;
    size_t scratch_capacity = 0;

    // Keep ordinary counting buffers across files without retaining arbitrarily
    // large outliers. This is a retention limit, not an input or counting limit:
    // an active genome may require more space, including when its hint is low.
    // 128 workers retain at most 12 GiB of occurrence capacity between ordinary
    // inputs, plus the separately bounded sort scratch below.
    static constexpr size_t retained_limit = size_t(12) << 20; // 96 MiB/worker

    void begin(size_t reserve_hint = 0) {
        const size_t maximum = occurrences.max_size();
        if (reserve_hint > maximum)
            throw std::length_error("ProbMinHash occurrence hint exceeds vector maximum");
        const size_t previous_size = occurrences.size();
        occurrences.clear();
        const size_t quantum = size_t(1) << 16;
        size_t wanted = std::max(reserve_hint, size_t(1) << 20);
        const size_t remainder = wanted % quantum;
        if (remainder)
            wanted += std::min(quantum - remainder, maximum - wanted);

        // Small/large ordinary genomes commonly alternate. Shrinking merely
        // because the next hint is smaller repeatedly discards reusable pages.
        // Only release an outlier once the next hint fits the retained range.
        // Compressed input hints can underestimate actual size; push_back still
        // grows as necessary and no occurrence is omitted.
        if (occurrences.capacity() > retained_limit && wanted <= retained_limit) {
            std::vector<uint64_t>().swap(occurrences);
            // A low compressed-file hint may make vector growth overshoot the
            // retention limit even for an ordinary input. Learn its actual
            // size so the next input does not repeat that overshoot.
            if (previous_size <= retained_limit)
                wanted = std::max(wanted, previous_size);
        }
        if (wanted <= occurrences.capacity()) return;

        size_t capacity = wanted;
        if (wanted <= retained_limit) {
            capacity = std::max(occurrences.capacity(), size_t(1) << 20);
            while (capacity < wanted) {
                const size_t growth = std::max(capacity / 2, quantum);
                capacity += std::min(growth, retained_limit - capacity);
            }
            const size_t tail = capacity % quantum;
            if (tail)
                capacity += std::min(quantum - tail, retained_limit - capacity);
        }
        occurrences.reserve(capacity);
    }

    static constexpr size_t scratch_limit = detail::BoundedRadixSort::scratch_limit;

    void sort() {
        const size_t n = occurrences.size();
        const size_t wanted = n < 128 ? 0 : std::min(n, scratch_limit);
        if (scratch_capacity < wanted) {
            scratch.reset(new uint64_t[wanted]);
            scratch_capacity = wanted;
        }
        detail::BoundedRadixSort::sort(occurrences.data(), n, scratch.get(), scratch_capacity);
    }

    void release() {
        std::vector<uint64_t>().swap(occurrences);
        scratch.reset();
        scratch_capacity = 0;
    }
};

// Returns the exact count for an accepted pair. A rejection is possible only
// when matching every still-unread register cannot meet the requested count.
// Byte fingerprints are an optional upper bound: equal nonzero uint64 winners
// must have equal low bytes. Collisions and zero winners can only add candidates,
// never remove a pair accepted by the exact comparison below.
inline bool winnerByteFilterAllows(const uint8_t* left, const uint8_t* right,
                                   size_t m, size_t minimum) {
    return Runtime::countEqualBytes(left, right, m) >= minimum;
}

inline bool countWinnersAtLeast(const uint64_t* left, const uint64_t* right,
                                size_t m, size_t minimum, size_t& common) {
    common = 0;
    if (minimum > m) return false;
    if (m <= 2048 || minimum == 0) {
        common = Runtime::countEqualNonZeroU64(left, right, m);
        return common >= minimum;
    }
    constexpr size_t BLOCK = 512;
    for (size_t offset = 0; offset < m;) {
        const size_t length = std::min(BLOCK, m - offset);
        common += Runtime::countEqualNonZeroU64(left + offset, right + offset, length);
        offset += length;
        if (common + (m - offset) < minimum) return false;
    }
    return common >= minimum;
}

} // namespace ProbMinHashBuild
} // namespace Sketch
#endif
