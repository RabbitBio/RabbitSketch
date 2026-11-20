#ifndef RABBITSKETCH_PROBMINHASH_BUILD_H
#define RABBITSKETCH_PROBMINHASH_BUILD_H

#include "api/RuntimeInfo.h"
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

    // Retain useful capacity, but do not pin the largest genome seen by every
    // worker for the rest of the run. Empty vectors can be replaced without
    // copying any k-mer data. Plain-input byte hints bound the required size.
    void begin(size_t reserve_hint = 0) {
        occurrences.clear();
        const size_t quantum = size_t(1) << 16;
        size_t wanted = std::max(reserve_hint, size_t(1) << 20);
        if (wanted <= std::numeric_limits<size_t>::max() - (quantum - 1))
            wanted = (wanted + quantum - 1) / quantum * quantum;
        if (occurrences.capacity() / 2 > wanted)
            std::vector<uint64_t>().swap(occurrences);
        if (wanted > occurrences.capacity()) occurrences.reserve(wanted);
    }

    static constexpr size_t scratch_limit = size_t(1) << 16;

    // Small buckets use cache-resident LSD sorting. Large buckets are first
    // partitioned in place on their most significant byte (American flag
    // sorting). All values and duplicate frequencies are retained exactly;
    // the final integer order is the same as std::sort.
    void sortRange(uint64_t* values, size_t n, unsigned highest_bit) {
        if (n < 128) {
            std::sort(values, values + n);
            return;
        }
        if (n <= scratch_limit) {
            constexpr unsigned BITS = 11;
            constexpr size_t BUCKETS = size_t(1) << BITS;
            std::array<size_t, BUCKETS> positions;
            uint64_t* source = values;
            uint64_t* target = scratch.get();
            for (unsigned shift = 0; shift <= highest_bit; shift += BITS) {
                positions.fill(0);
                for (size_t i = 0; i < n; ++i)
                    ++positions[(source[i] >> shift) & (BUCKETS - 1)];
                size_t next = 0;
                for (size_t bucket = 0; bucket < BUCKETS; ++bucket) {
                    const size_t count = positions[bucket];
                    positions[bucket] = next;
                    next += count;
                }
                for (size_t i = 0; i < n; ++i) {
                    const uint64_t key = source[i];
                    target[positions[(key >> shift) & (BUCKETS - 1)]++] = key;
                }
                std::swap(source, target);
            }
            if (source != values) std::copy_n(source, n, values);
            return;
        }
        constexpr size_t BUCKETS = 256;
        const unsigned shift = highest_bit > 7 ? highest_bit - 7 : 0;
        std::array<size_t, BUCKETS> end{}, next{};
        for (size_t i = 0; i < n; ++i) ++end[(values[i] >> shift) & 255];
        size_t offset = 0;
        for (size_t b = 0; b < BUCKETS; ++b) {
            next[b] = offset;
            offset += end[b];
            end[b] = offset;
        }
        for (size_t b = 0; b < BUCKETS; ++b) {
            while (next[b] < end[b]) {
                const size_t destination = (values[next[b]] >> shift) & 255;
                if (destination == b) ++next[b];
                else std::swap(values[next[b]], values[next[destination]++]);
            }
        }
        if (shift == 0) return;
        offset = 0;
        for (size_t b = 0; b < BUCKETS; ++b) {
            const size_t count = end[b] - offset;
            if (count > 1) sortRange(values + offset, count, shift - 1);
            offset = end[b];
        }
    }

    void sort() {
        const size_t n = occurrences.size();
        if (n < 128) {
            std::sort(occurrences.begin(), occurrences.end());
            return;
        }
        uint64_t differences = 0;
        const uint64_t first = occurrences.front();
        for (uint64_t key : occurrences) differences |= key ^ first;
        if (!differences) return;
        const size_t wanted = std::min(n, scratch_limit);
        if (scratch_capacity < wanted) {
            scratch.reset(new uint64_t[wanted]);
            scratch_capacity = wanted;
        }
        sortRange(occurrences.data(), n, 63u - static_cast<unsigned>(__builtin_clzll(differences)));
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
