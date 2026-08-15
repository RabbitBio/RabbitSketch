#include "api/RuntimeInfo.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>

#if (defined(__GNUC__) || defined(__clang__)) && \
    (defined(__x86_64__) || defined(__i386__))
#define RABBITSKETCH_X86_FUNCTION_DISPATCH 1
#define RABBITSKETCH_TARGET(value) __attribute__((target(value)))
#include <immintrin.h>
#else
#define RABBITSKETCH_X86_FUNCTION_DISPATCH 0
#define RABBITSKETCH_TARGET(value)
#endif

namespace Sketch {
namespace Runtime {
namespace {

using ByteFunction = ByteRelations (*)(const uint8_t*, const uint8_t*, size_t);
using U64Function = size_t (*)(const uint64_t*, const uint64_t*, size_t);
using IntersectionFunction = size_t (*)(const uint64_t*, size_t,
                                        const uint64_t*, size_t);
using RankFilterFunction = uint8_t (*)(const uint64_t*, uint8_t, uint64_t,
                                      uint64_t, uint64_t*);

ByteRelations scalarByteRelations(const uint8_t* left,
                                  const uint8_t* right,
                                  size_t size) noexcept {
    ByteRelations result;
    for (size_t index = 0; index < size; ++index) {
        if (left[index] == right[index]) ++result.equal;
        else if (left[index] > right[index]) ++result.greater;
        else ++result.less;
    }
    return result;
}

size_t scalarEqualNonZeroU64(const uint64_t* left,
                             const uint64_t* right,
                             size_t size) noexcept {
    size_t result = 0;
    for (size_t index = 0; index < size; ++index)
        result += left[index] != 0 && left[index] == right[index];
    return result;
}

size_t scalarSortedIntersectionU64(const uint64_t* left,
                                   size_t left_size,
                                   const uint64_t* right,
                                   size_t right_size) noexcept {
    size_t result = 0;
    size_t left_index = 0;
    size_t right_index = 0;
    while (left_index < left_size && right_index < right_size) {
        if (left[left_index] < right[right_index]) ++left_index;
        else if (right[right_index] < left[left_index]) ++right_index;
        else {
            ++result;
            ++left_index;
            ++right_index;
        }
    }
    return result;
}

inline uint64_t fmix64(uint64_t value, uint64_t seed) noexcept {
    value ^= seed;
    value ^= value >> 33;
    value *= UINT64_C(0xff51afd7ed558ccd);
    value ^= value >> 33;
    value *= UINT64_C(0xc4ceb9fe1a85ec53);
    value ^= value >> 33;
    return value;
}

uint8_t scalarFilterFmix64x8(const uint64_t* canonical_values,
                             uint8_t valid_mask,
                             uint64_t seed,
                             uint64_t threshold,
                             uint64_t* fingerprints) noexcept {
    uint8_t result = 0;
    for (unsigned lane = 0; lane < 8; ++lane) {
        if ((valid_mask & (UINT8_C(1) << lane)) == 0) continue;
        const uint64_t fingerprint = fmix64(canonical_values[lane], seed);
        fingerprints[lane] = fingerprint;
        if (fingerprint <= threshold)
            result |= static_cast<uint8_t>(UINT8_C(1) << lane);
    }
    return result;
}

#if RABBITSKETCH_X86_FUNCTION_DISPATCH
RABBITSKETCH_TARGET("sse2")
ByteRelations sse2ByteRelations(const uint8_t* left,
                                const uint8_t* right,
                                size_t size) noexcept {
    ByteRelations result;
    size_t index = 0;
    const __m128i ones = _mm_set1_epi8(static_cast<char>(0xff));
    for (; index + 16 <= size; index += 16) {
        const __m128i a = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(left + index));
        const __m128i b = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(right + index));
        const __m128i equal = _mm_cmpeq_epi8(a, b);
        const __m128i a_is_max = _mm_cmpeq_epi8(a, _mm_max_epu8(a, b));
        const __m128i greater = _mm_andnot_si128(equal, a_is_max);
        const __m128i less = _mm_andnot_si128(
            equal, _mm_andnot_si128(a_is_max, ones));
        result.equal += static_cast<size_t>(__builtin_popcount(
            static_cast<unsigned>(_mm_movemask_epi8(equal))));
        result.greater += static_cast<size_t>(__builtin_popcount(
            static_cast<unsigned>(_mm_movemask_epi8(greater))));
        result.less += static_cast<size_t>(__builtin_popcount(
            static_cast<unsigned>(_mm_movemask_epi8(less))));
    }
    const ByteRelations tail = scalarByteRelations(
        left + index, right + index, size - index);
    result.equal += tail.equal;
    result.greater += tail.greater;
    result.less += tail.less;
    return result;
}

RABBITSKETCH_TARGET("avx2")
ByteRelations avx2ByteRelations(const uint8_t* left,
                                const uint8_t* right,
                                size_t size) noexcept {
    ByteRelations result;
    size_t index = 0;
    const __m256i ones = _mm256_set1_epi8(static_cast<char>(0xff));
    for (; index + 32 <= size; index += 32) {
        const __m256i a = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(left + index));
        const __m256i b = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(right + index));
        const __m256i equal = _mm256_cmpeq_epi8(a, b);
        const __m256i a_is_max = _mm256_cmpeq_epi8(a, _mm256_max_epu8(a, b));
        const __m256i greater = _mm256_andnot_si256(equal, a_is_max);
        const __m256i less = _mm256_andnot_si256(
            equal, _mm256_andnot_si256(a_is_max, ones));
        result.equal += static_cast<size_t>(__builtin_popcount(
            static_cast<unsigned>(_mm256_movemask_epi8(equal))));
        result.greater += static_cast<size_t>(__builtin_popcount(
            static_cast<unsigned>(_mm256_movemask_epi8(greater))));
        result.less += static_cast<size_t>(__builtin_popcount(
            static_cast<unsigned>(_mm256_movemask_epi8(less))));
    }
    const ByteRelations tail = scalarByteRelations(
        left + index, right + index, size - index);
    result.equal += tail.equal;
    result.greater += tail.greater;
    result.less += tail.less;
    return result;
}

RABBITSKETCH_TARGET("avx512f,avx512bw")
ByteRelations avx512ByteRelations(const uint8_t* left,
                                  const uint8_t* right,
                                  size_t size) noexcept {
    ByteRelations result;
    size_t index = 0;
    for (; index + 64 <= size; index += 64) {
        const __m512i a = _mm512_loadu_si512(left + index);
        const __m512i b = _mm512_loadu_si512(right + index);
        const uint64_t equal = _mm512_cmpeq_epu8_mask(a, b);
        const uint64_t greater = _mm512_cmp_epu8_mask(
            a, b, _MM_CMPINT_GT);
        const uint64_t less = _mm512_cmp_epu8_mask(
            a, b, _MM_CMPINT_LT);
        result.equal += static_cast<size_t>(__builtin_popcountll(equal));
        result.greater += static_cast<size_t>(__builtin_popcountll(greater));
        result.less += static_cast<size_t>(__builtin_popcountll(less));
    }
    const ByteRelations tail = scalarByteRelations(
        left + index, right + index, size - index);
    result.equal += tail.equal;
    result.greater += tail.greater;
    result.less += tail.less;
    return result;
}

RABBITSKETCH_TARGET("sse4.1")
size_t sse41EqualNonZeroU64(const uint64_t* left,
                            const uint64_t* right,
                            size_t size) noexcept {
    size_t result = 0;
    size_t index = 0;
    const __m128i zero = _mm_setzero_si128();
    for (; index + 2 <= size; index += 2) {
        const __m128i a = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(left + index));
        const __m128i b = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(right + index));
        const __m128i equal = _mm_cmpeq_epi64(a, b);
        const __m128i nonzero = _mm_andnot_si128(_mm_cmpeq_epi64(a, zero),
                                                 equal);
        result += static_cast<size_t>(__builtin_popcount(
            static_cast<unsigned>(_mm_movemask_pd(
                _mm_castsi128_pd(nonzero)))));
    }
    return result + scalarEqualNonZeroU64(
        left + index, right + index, size - index);
}

RABBITSKETCH_TARGET("avx2")
size_t avx2EqualNonZeroU64(const uint64_t* left,
                           const uint64_t* right,
                           size_t size) noexcept {
    size_t result = 0;
    size_t index = 0;
    const __m256i zero = _mm256_setzero_si256();
    for (; index + 4 <= size; index += 4) {
        const __m256i a = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(left + index));
        const __m256i b = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(right + index));
        const __m256i equal = _mm256_cmpeq_epi64(a, b);
        const __m256i nonzero = _mm256_andnot_si256(
            _mm256_cmpeq_epi64(a, zero), equal);
        result += static_cast<size_t>(__builtin_popcount(
            static_cast<unsigned>(_mm256_movemask_pd(
                _mm256_castsi256_pd(nonzero)))));
    }
    return result + scalarEqualNonZeroU64(
        left + index, right + index, size - index);
}

RABBITSKETCH_TARGET("avx512f")
size_t avx512EqualNonZeroU64(const uint64_t* left,
                             const uint64_t* right,
                             size_t size) noexcept {
    size_t result = 0;
    size_t index = 0;
    const __m512i zero = _mm512_setzero_si512();
    for (; index + 8 <= size; index += 8) {
        const __m512i a = _mm512_loadu_si512(left + index);
        const __m512i b = _mm512_loadu_si512(right + index);
        const uint8_t mask = _mm512_cmpeq_epu64_mask(a, b) &
            static_cast<uint8_t>(~_mm512_cmpeq_epu64_mask(a, zero));
        result += static_cast<size_t>(__builtin_popcount(mask));
    }
    return result + scalarEqualNonZeroU64(
        left + index, right + index, size - index);
}

RABBITSKETCH_TARGET("avx2")
__m256i avx2MulloEpi64(__m256i left, __m256i right) noexcept {
    const __m256i left_high = _mm256_srli_epi64(left, 32);
    const __m256i right_high = _mm256_srli_epi64(right, 32);
    const __m256i low = _mm256_mul_epu32(left, right);
    const __m256i cross = _mm256_add_epi64(
        _mm256_mul_epu32(left_high, right),
        _mm256_mul_epu32(left, right_high));
    return _mm256_add_epi64(low, _mm256_slli_epi64(cross, 32));
}

RABBITSKETCH_TARGET("avx2")
size_t avx2SortedIntersectionU64(const uint64_t* left,
                                 size_t left_size,
                                 const uint64_t* right,
                                 size_t right_size) noexcept {
    size_t result = 0;
    size_t left_index = 0;
    size_t right_index = 0;
    while (left_index + 4 <= left_size &&
           right_index + 4 <= right_size) {
        const __m256i values = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(left + left_index));
        __m256i matches = _mm256_setzero_si256();
        for (size_t lane = 0; lane < 4; ++lane) {
            const __m256i candidate = _mm256_set1_epi64x(
                static_cast<long long>(right[right_index + lane]));
            matches = _mm256_or_si256(
                matches, _mm256_cmpeq_epi64(values, candidate));
        }
        result += static_cast<size_t>(__builtin_popcount(
            static_cast<unsigned>(_mm256_movemask_pd(
                _mm256_castsi256_pd(matches)))));
        const uint64_t left_max = left[left_index + 3];
        const uint64_t right_max = right[right_index + 3];
        if (left_max <= right_max) left_index += 4;
        if (right_max <= left_max) right_index += 4;
    }
    return result + scalarSortedIntersectionU64(
        left + left_index, left_size - left_index,
        right + right_index, right_size - right_index);
}

RABBITSKETCH_TARGET("avx512f")
size_t avx512SortedIntersectionU64(const uint64_t* left,
                                   size_t left_size,
                                   const uint64_t* right,
                                   size_t right_size) noexcept {
    size_t result = 0;
    size_t left_index = 0;
    size_t right_index = 0;
    while (left_index + 8 <= left_size &&
           right_index + 8 <= right_size) {
        const __m512i values = _mm512_loadu_si512(left + left_index);
        uint8_t matches = 0;
        for (size_t lane = 0; lane < 8; ++lane) {
            const __m512i candidate = _mm512_set1_epi64(
                static_cast<long long>(right[right_index + lane]));
            matches |= static_cast<uint8_t>(
                _mm512_cmpeq_epu64_mask(values, candidate));
        }
        result += static_cast<size_t>(__builtin_popcount(
            static_cast<unsigned>(matches)));
        const uint64_t left_max = left[left_index + 7];
        const uint64_t right_max = right[right_index + 7];
        if (left_max <= right_max) left_index += 8;
        if (right_max <= left_max) right_index += 8;
    }
    return result + scalarSortedIntersectionU64(
        left + left_index, left_size - left_index,
        right + right_index, right_size - right_index);
}

RABBITSKETCH_TARGET("avx2")
uint8_t avx2FilterFmix64x8(const uint64_t* canonical_values,
                           uint8_t valid_mask,
                           uint64_t seed,
                           uint64_t threshold,
                           uint64_t* fingerprints) noexcept {
    uint8_t result = 0;
    const __m256i seed_vector = _mm256_set1_epi64x(
        static_cast<long long>(seed));
    const __m256i first_constant = _mm256_set1_epi64x(
        static_cast<long long>(UINT64_C(0xff51afd7ed558ccd)));
    const __m256i second_constant = _mm256_set1_epi64x(
        static_cast<long long>(UINT64_C(0xc4ceb9fe1a85ec53)));
    const __m256i sign = _mm256_set1_epi64x(
        static_cast<long long>(UINT64_C(0x8000000000000000)));
    const __m256i threshold_vector = _mm256_xor_si256(
        _mm256_set1_epi64x(static_cast<long long>(threshold)), sign);
    for (unsigned base = 0; base < 8; base += 4) {
        __m256i value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(canonical_values + base));
        value = _mm256_xor_si256(value, seed_vector);
        value = _mm256_xor_si256(value, _mm256_srli_epi64(value, 33));
        value = avx2MulloEpi64(value, first_constant);
        value = _mm256_xor_si256(value, _mm256_srli_epi64(value, 33));
        value = avx2MulloEpi64(value, second_constant);
        value = _mm256_xor_si256(value, _mm256_srli_epi64(value, 33));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(fingerprints + base), value);
        const __m256i greater = _mm256_cmpgt_epi64(
            _mm256_xor_si256(value, sign), threshold_vector);
        const unsigned rejected = static_cast<unsigned>(
            _mm256_movemask_pd(_mm256_castsi256_pd(greater)));
        const unsigned lanes = (valid_mask >> base) & 0x0fu;
        result |= static_cast<uint8_t>((lanes & ~rejected) << base);
    }
    return result;
}

RABBITSKETCH_TARGET("avx512f,avx512dq")
uint8_t avx512FilterFmix64x8(const uint64_t* canonical_values,
                             uint8_t valid_mask,
                             uint64_t seed,
                             uint64_t threshold,
                             uint64_t* fingerprints) noexcept {
    __m512i value = _mm512_loadu_si512(canonical_values);
    value = _mm512_xor_epi64(
        value, _mm512_set1_epi64(static_cast<long long>(seed)));
    value = _mm512_xor_epi64(value, _mm512_srli_epi64(value, 33));
    value = _mm512_mullo_epi64(value, _mm512_set1_epi64(
        static_cast<long long>(UINT64_C(0xff51afd7ed558ccd))));
    value = _mm512_xor_epi64(value, _mm512_srli_epi64(value, 33));
    value = _mm512_mullo_epi64(value, _mm512_set1_epi64(
        static_cast<long long>(UINT64_C(0xc4ceb9fe1a85ec53))));
    value = _mm512_xor_epi64(value, _mm512_srli_epi64(value, 33));
    _mm512_storeu_si512(fingerprints, value);
    return static_cast<uint8_t>(_mm512_mask_cmp_epu64_mask(
        static_cast<__mmask8>(valid_mask), value,
        _mm512_set1_epi64(static_cast<long long>(threshold)),
        _MM_CMPINT_LE));
}
#endif

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
        [](unsigned char character) {
            return static_cast<char>(std::tolower(character));
        });
    return value;
}

std::string architectureName() {
#if defined(__x86_64__) || defined(_M_X64)
    return "x86_64";
#elif defined(__i386__) || defined(_M_IX86)
    return "x86";
#elif defined(__aarch64__) || defined(_M_ARM64)
    return "aarch64";
#elif defined(__arm__) || defined(_M_ARM)
    return "arm";
#else
    return "unknown";
#endif
}

std::string compilerName() {
#if defined(__clang__)
    return std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("gcc ") + __VERSION__;
#elif defined(_MSC_VER)
    return std::string("msvc ") + std::to_string(_MSC_VER);
#else
    return "unknown";
#endif
}

struct Dispatch {
    RuntimeInfo info;
    ByteFunction bytes = scalarByteRelations;
    U64Function u64 = scalarEqualNonZeroU64;
    IntersectionFunction intersection = scalarSortedIntersectionU64;
    RankFilterFunction rank_filter = scalarFilterFmix64x8;
};

Dispatch makeDispatch() noexcept {
    Dispatch result;
    result.info.architecture = architectureName();
    result.info.compiler = compilerName();
#ifdef RABBITSKETCH_NATIVE_ARCH
    result.info.portable_baseline = false;
#endif
#if RABBITSKETCH_X86_FUNCTION_DISPATCH
    result.info.runtime_dispatch_available = true;
    __builtin_cpu_init();
    result.info.cpu_sse2 = __builtin_cpu_supports("sse2");
    result.info.cpu_sse41 = __builtin_cpu_supports("sse4.1");
    result.info.cpu_avx2 = __builtin_cpu_supports("avx2");
    result.info.cpu_avx512f = __builtin_cpu_supports("avx512f");
    result.info.cpu_avx512dq = __builtin_cpu_supports("avx512dq");
    result.info.cpu_avx512bw = __builtin_cpu_supports("avx512bw");
#endif
    const char* raw_override = std::getenv("RABBITSKETCH_SIMD");
    const std::string override_value = raw_override == nullptr
        ? "auto" : lower(raw_override);
    if (raw_override != nullptr)
        result.info.environment_override = raw_override;

    const bool explicit_scalar = override_value == "scalar";
    const bool allow_sse = override_value == "auto" ||
                           override_value == "sse2" ||
                           override_value == "sse4.1";
    const bool allow_avx2 = override_value == "auto" ||
                            override_value == "avx2";
    const bool allow_avx512 = override_value == "auto" ||
                              override_value == "avx512" ||
                              override_value == "avx512bw";
    const bool known_override = override_value == "auto" || explicit_scalar ||
        override_value == "sse2" || override_value == "sse4.1" ||
        override_value == "avx2" || override_value == "avx512" ||
        override_value == "avx512bw";

#if RABBITSKETCH_X86_FUNCTION_DISPATCH
    if (!explicit_scalar && allow_avx512 && result.info.cpu_avx512f &&
        result.info.cpu_avx512bw) {
        result.bytes = avx512ByteRelations;
        result.u64 = avx512EqualNonZeroU64;
        result.intersection = avx512SortedIntersectionU64;
        result.info.selected_byte_path = "avx512bw";
        result.info.selected_u64_path = "avx512f";
        if (override_value != "auto" && result.info.cpu_avx512dq) {
            result.rank_filter = avx512FilterFmix64x8;
            result.info.selected_rank_path = "avx512dq";
        } else if (result.info.cpu_avx2) {
            result.rank_filter = avx2FilterFmix64x8;
            result.info.selected_rank_path = "avx2";
        } else if (result.info.cpu_avx512dq) {
            result.rank_filter = avx512FilterFmix64x8;
            result.info.selected_rank_path = "avx512dq";
        }
    } else if (!explicit_scalar && allow_avx2 && result.info.cpu_avx2) {
        result.bytes = avx2ByteRelations;
        result.u64 = avx2EqualNonZeroU64;
        result.intersection = avx2SortedIntersectionU64;
        result.rank_filter = avx2FilterFmix64x8;
        result.info.selected_byte_path = "avx2";
        result.info.selected_u64_path = "avx2";
        result.info.selected_rank_path = "avx2";
    } else if (!explicit_scalar && allow_sse && result.info.cpu_sse2) {
        result.bytes = sse2ByteRelations;
        result.info.selected_byte_path = "sse2";
        if (result.info.cpu_sse41) {
            result.u64 = sse41EqualNonZeroU64;
            result.info.selected_u64_path = "sse4.1";
        }
    }
#endif
    if (raw_override != nullptr && override_value != "auto") {
        result.info.environment_override_honored = known_override &&
            ((explicit_scalar && result.info.selected_byte_path == "scalar") ||
             (override_value == "sse2" &&
                  result.info.selected_byte_path == "sse2") ||
             (override_value == "sse4.1" &&
                  result.info.selected_u64_path == "sse4.1") ||
             (override_value == "avx2" &&
                  result.info.selected_byte_path == "avx2") ||
             ((override_value == "avx512" ||
               override_value == "avx512bw") &&
                  result.info.selected_byte_path == "avx512bw"));
    }
    return result;
}

const Dispatch& dispatch() noexcept {
    static const Dispatch selected = makeDispatch();
    return selected;
}

} // namespace

const RuntimeInfo& runtimeInfo() noexcept {
    return dispatch().info;
}

ByteRelations countByteRelations(const uint8_t* left,
                                 const uint8_t* right,
                                 size_t size) noexcept {
    if (size == 0) return ByteRelations();
    if (left == nullptr || right == nullptr) return ByteRelations();
    return dispatch().bytes(left, right, size);
}

size_t countEqualBytes(const uint8_t* left,
                       const uint8_t* right,
                       size_t size) noexcept {
    return countByteRelations(left, right, size).equal;
}

size_t countEqualNonZeroU64(const uint64_t* left,
                            const uint64_t* right,
                            size_t size) noexcept {
    if (size == 0 || left == nullptr || right == nullptr) return 0;
    return dispatch().u64(left, right, size);
}

size_t countSortedIntersectionU64(const uint64_t* left,
                                  size_t left_size,
                                  const uint64_t* right,
                                  size_t right_size) noexcept {
    if (left_size == 0 || right_size == 0 ||
        left == nullptr || right == nullptr)
        return 0;
    if (std::min(left_size, right_size) < 64)
        return scalarSortedIntersectionU64(
            left, left_size, right, right_size);

    size_t common = 0;
    size_t left_index = 0;
    size_t right_index = 0;
    size_t steps = 0;
    constexpr size_t SAMPLE_STEPS = 256;
    while (left_index < left_size && right_index < right_size &&
           steps < SAMPLE_STEPS) {
        if (left[left_index] < right[right_index]) ++left_index;
        else if (right[right_index] < left[left_index]) ++right_index;
        else {
            ++common;
            ++left_index;
            ++right_index;
        }
        ++steps;
    }
    const bool high_but_shifted_overlap = common != steps &&
        common * 2 >= steps;
    const IntersectionFunction remainder = high_but_shifted_overlap
        ? scalarSortedIntersectionU64 : dispatch().intersection;
    return common + remainder(
        left + left_index, left_size - left_index,
        right + right_index, right_size - right_index);
}

uint8_t filterFmix64x8(const uint64_t* canonical_values,
                       uint8_t valid_mask,
                       uint64_t seed,
                       uint64_t threshold,
                       uint64_t* fingerprints) noexcept {
    if (valid_mask == 0 || canonical_values == nullptr ||
        fingerprints == nullptr)
        return 0;
    return dispatch().rank_filter(canonical_values, valid_mask, seed,
                                  threshold, fingerprints);
}

} // namespace Runtime
} // namespace Sketch
