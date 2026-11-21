#ifndef RABBITSKETCH_BOUNDED_RADIX_SORT_H
#define RABBITSKETCH_BOUNDED_RADIX_SORT_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace Sketch {
namespace detail {

// Exact unsigned-integer sorting with bounded caller-owned scratch. Large
// ranges are partitioned in place; small partitions use stable LSD passes.
// Equal values and their multiplicities are preserved. No heap allocation is
// performed here. For n >= 128, scratch must be disjoint from the input range.
class BoundedRadixSort {
public:
    static constexpr size_t scratch_limit = size_t(1) << 16;

    static void sort(uint64_t* values, size_t n,
                     uint64_t* scratch, size_t scratch_capacity) {
        sortWithPartitionBits<9>(values, n, scratch, scratch_capacity);
    }

    // Both widths preserve the same exact order and scratch bound. Keeping the
    // selection explicit permits controlled benchmarks without runtime state.
    template<unsigned PartitionBits>
    static void sortWithPartitionBits(uint64_t* values, size_t n,
                                      uint64_t* scratch, size_t scratch_capacity) {
        static_assert(PartitionBits == 8 || PartitionBits == 9,
                      "The bounded radix partition width must be 8 or 9 bits");
        if (n < 2) return;
        if (!values) throw std::invalid_argument("Null bounded radix input");
        if (n < 128) {
            std::sort(values, values + n);
            return;
        }
        if (!scratch || scratch_capacity < std::min(n, scratch_limit))
            throw std::invalid_argument("Insufficient bounded radix scratch");

        uint64_t differences = 0;
        const uint64_t first = values[0];
        for (size_t i = 1; i < n; ++i) differences |= values[i] ^ first;
        if (!differences) return;
        const unsigned highest_bit = 63u - static_cast<unsigned>(__builtin_clzll(differences));
        sortRange<PartitionBits>(values, n, highest_bit, scratch);
    }

private:
    static constexpr unsigned digit_bits = 11;
    static constexpr size_t digit_buckets = size_t(1) << digit_bits;
    static constexpr uint64_t digit_mask = digit_buckets - 1;

    template<unsigned Digits>
    static void sortLeaf(uint64_t* values, size_t n, uint64_t* scratch) {
        static_assert(Digits >= 1 && Digits <= 6, "Invalid uint64 radix digit count");
        // At most 65,536 values reach a leaf. uint32_t holds both histogram
        // counts and exclusive prefix positions. The largest histogram is
        // 6 * 2048 * 4 = 48 KiB, independent of the complete input length.
        std::array<std::array<uint32_t, digit_buckets>, Digits> positions{};
        // Build every digit histogram in one input scan. Digit frequencies do
        // not change when another digit is sorted, so subsequent scatter passes
        // can reuse these counts without rescanning to rebuild each histogram.
        for (size_t i = 0; i < n; ++i) {
            uint64_t key = values[i];
            for (unsigned digit = 0; digit < Digits; ++digit) {
                ++positions[digit][key & digit_mask];
                key >>= digit_bits;
            }
        }

        uint64_t* source = values;
        uint64_t* target = scratch;
        for (unsigned digit = 0; digit < Digits; ++digit) {
            auto& histogram = positions[digit];
            // A constant digit cannot change order. Skipping it is also safe
            // for a leaf whose upper bits are fixed by preceding MSD buckets.
            bool constant = false;
            uint32_t next = 0;
            for (uint32_t& count : histogram) {
                const uint32_t frequency = count;
                constant |= frequency == n;
                count = next;
                next += frequency;
            }
            if (constant) continue;
            const unsigned shift = digit * digit_bits;
            for (size_t i = 0; i < n; ++i) {
                const uint64_t key = source[i];
                target[histogram[(key >> shift) & digit_mask]++] = key;
            }
            std::swap(source, target);
        }
        if (source != values) std::copy_n(source, n, values);
    }

    static void sortLeaf(uint64_t* values, size_t n,
                         unsigned highest_bit, uint64_t* scratch) {
        switch (highest_bit / digit_bits) {
            case 0: sortLeaf<1>(values, n, scratch); break;
            case 1: sortLeaf<2>(values, n, scratch); break;
            case 2: sortLeaf<3>(values, n, scratch); break;
            case 3: sortLeaf<4>(values, n, scratch); break;
            case 4: sortLeaf<5>(values, n, scratch); break;
            default: sortLeaf<6>(values, n, scratch); break;
        }
    }

    template<unsigned PartitionBits>
    static void sortRange(uint64_t* values, size_t n,
                          unsigned highest_bit, uint64_t* scratch) {
        if (n < 128) {
            std::sort(values, values + n);
            return;
        }
        if (n <= scratch_limit) {
            sortLeaf(values, n, highest_bit, scratch);
            return;
        }
        constexpr size_t buckets = size_t(1) << PartitionBits;
        constexpr uint64_t mask = buckets - 1;
        const unsigned shift = highest_bit >= PartitionBits ? highest_bit + 1 - PartitionBits : 0;
        std::array<size_t, buckets> end{}, next{};
        for (size_t i = 0; i < n; ++i) ++end[(values[i] >> shift) & mask];
        size_t offset = 0;
        for (size_t bucket = 0; bucket < buckets; ++bucket) {
            next[bucket] = offset;
            offset += end[bucket];
            end[bucket] = offset;
        }
        for (size_t bucket = 0; bucket < buckets; ++bucket) {
            while (next[bucket] < end[bucket]) {
                const size_t destination = (values[next[bucket]] >> shift) & mask;
                if (destination == bucket) ++next[bucket];
                else std::swap(values[next[bucket]], values[next[destination]++]);
            }
        }
        if (shift == 0) return;
        offset = 0;
        for (size_t bucket = 0; bucket < buckets; ++bucket) {
            const size_t count = end[bucket] - offset;
            if (count > 1) sortRange<PartitionBits>(values + offset, count, shift - 1, scratch);
            offset = end[bucket];
        }
    }
};

} // namespace detail
} // namespace Sketch
#endif
