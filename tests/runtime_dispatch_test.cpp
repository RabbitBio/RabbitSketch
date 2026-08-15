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
}

void testRuntimeReport() {
    const Sketch::Runtime::RuntimeInfo& info =
        Sketch::Runtime::runtimeInfo();
    check(!info.architecture.empty() && !info.compiler.empty() &&
              !info.selected_byte_path.empty() &&
              !info.selected_u64_path.empty(),
          "runtime capability report is incomplete");
    const char* override_value = std::getenv("RABBITSKETCH_SIMD");
    if (override_value != nullptr && std::string(override_value) == "scalar")
        check(info.environment_override_honored &&
                  info.selected_byte_path == "scalar" &&
                  info.selected_u64_path == "scalar",
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
