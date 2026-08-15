#include "api/RuntimeInfo.h"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

uint64_t fmix64(uint64_t value, uint64_t seed) {
    value ^= seed;
    value ^= value >> 33;
    value *= UINT64_C(0xff51afd7ed558ccd);
    value ^= value >> 33;
    value *= UINT64_C(0xc4ceb9fe1a85ec53);
    value ^= value >> 33;
    return value;
}

size_t scalarIntersection(const std::vector<uint64_t>& left,
                          const std::vector<uint64_t>& right) {
    size_t common = 0;
    size_t left_index = 0;
    size_t right_index = 0;
    while (left_index < left.size() && right_index < right.size()) {
        if (left[left_index] < right[right_index]) ++left_index;
        else if (right[right_index] < left[left_index]) ++right_index;
        else {
            ++common;
            ++left_index;
            ++right_index;
        }
    }
    return common;
}

void testDispatchedPrimitives() {
    std::vector<uint8_t> left(259);
    std::vector<uint8_t> right(259);
    Sketch::Runtime::ByteRelations expected;
    for (size_t index = 0; index < left.size(); ++index) {
        left[index] = static_cast<uint8_t>((index * 17) % 63);
        right[index] = static_cast<uint8_t>((index * 29) % 63);
        if (left[index] == right[index]) ++expected.equal;
        else if (left[index] > right[index]) ++expected.greater;
        else ++expected.less;
    }
    const Sketch::Runtime::ByteRelations observed =
        Sketch::Runtime::countByteRelations(
            left.data(), right.data(), left.size());
    check(observed.equal == expected.equal &&
              observed.greater == expected.greater &&
              observed.less == expected.less &&
              Sketch::Runtime::countEqualBytes(
                  left.data(), right.data(), left.size()) == expected.equal,
          "runtime-dispatched byte relations differ from scalar truth");

    std::vector<uint64_t> first(67);
    std::vector<uint64_t> second(67);
    size_t equal_nonzero = 0;
    for (size_t index = 0; index < first.size(); ++index) {
        first[index] = index % 11 == 0 ? 0 : index * UINT64_C(1000003);
        second[index] = index % 3 == 0 ? first[index] : first[index] + 1;
        equal_nonzero += first[index] != 0 && first[index] == second[index];
    }
    check(Sketch::Runtime::countEqualNonZeroU64(
              first.data(), second.data(), first.size()) == equal_nonzero,
          "runtime-dispatched uint64 equality differs from scalar truth");

    std::vector<uint64_t> sorted_left(513);
    std::vector<uint64_t> sorted_right(389);
    for (size_t index = 0; index < sorted_left.size(); ++index)
        sorted_left[index] = index * 3;
    for (size_t index = 0; index < sorted_right.size(); ++index)
        sorted_right[index] = index * 5;
    size_t expected_intersection = 0;
    for (uint64_t value : sorted_left)
        expected_intersection += value % 5 == 0 &&
            value <= sorted_right.back();
    check(Sketch::Runtime::countSortedIntersectionU64(
              sorted_left.data(), sorted_left.size(),
              sorted_right.data(), sorted_right.size()) ==
                  expected_intersection,
          "runtime-dispatched sorted intersection differs from scalar truth");

    uint64_t state = UINT64_C(0x9e3779b97f4a7c15);
    for (size_t trial = 0; trial < 128; ++trial) {
        std::vector<uint64_t> random_left;
        std::vector<uint64_t> random_right;
        uint64_t value = trial;
        state = state * UINT64_C(6364136223846793005) + 1;
        const size_t universe_size = 8 + state % 800;
        for (size_t index = 0; index < universe_size; ++index) {
            state = state * UINT64_C(6364136223846793005) + 1;
            value += 1 + (state & 15);
            if (state & 0x100) random_left.push_back(value);
            if (state & 0x200) random_right.push_back(value);
        }
        check(Sketch::Runtime::countSortedIntersectionU64(
                  random_left.data(), random_left.size(),
                  random_right.data(), random_right.size()) ==
                      scalarIntersection(random_left, random_right),
              "runtime-dispatched sorted intersection failed fuzz case");
    }

    const uint64_t canonical[8] = {
        0, 1, 2, 3, UINT64_C(0x123456789abcdef0),
        UINT64_MAX, UINT64_C(0xaaaaaaaaaaaaaaaa),
        UINT64_C(0x5555555555555555)};
    const uint64_t equality_threshold = fmix64(canonical[3], 42);
    for (uint8_t valid_mask : {UINT8_C(0xff), UINT8_C(0x55), UINT8_C(0x80)}) {
        for (uint64_t threshold : {UINT64_C(0), UINT64_MAX / 3,
                                   equality_threshold, UINT64_MAX}) {
            uint64_t observed_hashes[8] = {};
            const uint8_t observed_mask = Sketch::Runtime::filterFmix64x8(
                canonical, valid_mask, 42, threshold, observed_hashes);
            uint8_t expected_mask = 0;
            for (unsigned lane = 0; lane < 8; ++lane) {
                if ((valid_mask & (UINT8_C(1) << lane)) == 0) continue;
                const uint64_t expected_hash = fmix64(canonical[lane], 42);
                check(observed_hashes[lane] == expected_hash,
                      "runtime-dispatched fmix changed a fingerprint");
                if (expected_hash <= threshold)
                    expected_mask |= static_cast<uint8_t>(UINT8_C(1) << lane);
            }
            check(observed_mask == expected_mask,
                  "runtime-dispatched fmix threshold mask is incorrect");
        }
    }
}

void testRuntimeReport() {
    const Sketch::Runtime::RuntimeInfo& info =
        Sketch::Runtime::runtimeInfo();
    check(!info.architecture.empty() && !info.compiler.empty() &&
              !info.selected_byte_path.empty() &&
              !info.selected_u64_path.empty() &&
              !info.selected_rank_path.empty(),
          "runtime capability report is incomplete");
    const char* override_value = std::getenv("RABBITSKETCH_SIMD");
    if (override_value != nullptr && std::string(override_value) == "scalar")
        check(info.environment_override_honored &&
                  info.selected_byte_path == "scalar" &&
                  info.selected_u64_path == "scalar" &&
                  info.selected_rank_path == "scalar",
              "RABBITSKETCH_SIMD=scalar override was not honored");
}

} // namespace

int main() {
    try {
        testDispatchedPrimitives();
        testRuntimeReport();
        std::cout << "runtime_dispatch_test: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "runtime_dispatch_test: FAIL: " << error.what() << '\n';
        return 1;
    }
}
