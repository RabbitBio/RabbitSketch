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
        result.info.selected_byte_path = "avx512bw";
        result.info.selected_u64_path = "avx512f";
    } else if (!explicit_scalar && allow_avx2 && result.info.cpu_avx2) {
        result.bytes = avx2ByteRelations;
        result.u64 = avx2EqualNonZeroU64;
        result.info.selected_byte_path = "avx2";
        result.info.selected_u64_path = "avx2";
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

} // namespace Runtime
} // namespace Sketch
