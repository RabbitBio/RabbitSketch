#ifndef RABBITSKETCH_API_RUNTIME_INFO_H
#define RABBITSKETCH_API_RUNTIME_INFO_H

#include <cstddef>
#include <cstdint>
#include <string>

namespace Sketch {
namespace Runtime {

struct RuntimeInfo {
    std::string architecture;
    std::string compiler;
    bool portable_baseline = true;
    bool runtime_dispatch_available = false;
    bool cpu_sse2 = false;
    bool cpu_sse41 = false;
    bool cpu_avx2 = false;
    bool cpu_avx512f = false;
    bool cpu_avx512dq = false;
    bool cpu_avx512bw = false;
    std::string selected_byte_path = "scalar";
    std::string selected_u64_path = "scalar";
    std::string selected_rank_path = "scalar";
    std::string environment_override;
    bool environment_override_honored = true;
};

struct ByteRelations {
    size_t equal = 0;
    size_t greater = 0;
    size_t less = 0;
};

/** CPU features and the immutable dispatch decision for this process. */
const RuntimeInfo& runtimeInfo() noexcept;

/** Runtime-dispatched primitives used by register and weighted queries. */
ByteRelations countByteRelations(const uint8_t* left,
                                 const uint8_t* right,
                                 size_t size) noexcept;
size_t countEqualBytes(const uint8_t* left,
                       const uint8_t* right,
                       size_t size) noexcept;
size_t countEqualNonZeroU64(const uint64_t* left,
                            const uint64_t* right,
                            size_t size) noexcept;

/** Count matches between two sorted, duplicate-free uint64 arrays. */
size_t countSortedIntersectionU64(const uint64_t* left,
                                  size_t left_size,
                                  const uint64_t* right,
                                  size_t right_size) noexcept;

/**
 * Compute RankStream fmix64 for eight canonical values and return the lanes
 * whose unsigned fingerprint is <= threshold. Only valid-mask lanes are
 * considered; output values for invalid lanes are unspecified.
 */
uint8_t filterFmix64x8(const uint64_t* canonical_values,
                       uint8_t valid_mask,
                       uint64_t seed,
                       uint64_t threshold,
                       uint64_t* fingerprints) noexcept;

} // namespace Runtime
} // namespace Sketch

#endif // RABBITSKETCH_API_RUNTIME_INFO_H
