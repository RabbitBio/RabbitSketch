#include "BoundedRadixSort.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Sorter = Sketch::detail::BoundedRadixSort;
constexpr uint64_t GUARD = UINT64_C(0xd3c57e91a64fb208);
constexpr size_t GUARDS = 8;

void check(bool okay, const std::string& message) {
    if (!okay) throw std::runtime_error(message);
}

void call_sort(unsigned mode, uint64_t* values, size_t n,
               uint64_t* scratch, size_t capacity) {
    if (mode == 9) Sorter::sortWithPartitionBits<9>(values, n, scratch, capacity);
    else if (mode == 8) Sorter::sortWithPartitionBits<8>(values, n, scratch, capacity);
    else Sorter::sort(values, n, scratch, capacity);
}

void check_guards(const std::vector<uint64_t>& allocation, size_t size,
                  const std::string& name) {
    for (size_t i = 0; i < GUARDS; ++i) {
        check(allocation[i] == GUARD, name + " prefix guard was overwritten");
        check(allocation[GUARDS + size + i] == GUARD,
              name + " suffix guard was overwritten");
    }
}

void verify(const std::vector<uint64_t>& input, const std::string& label, unsigned mode) {
    auto expected = input;
    std::sort(expected.begin(), expected.end());
    const size_t n = input.size();
    const size_t capacity = n < 128 ? 0 : std::min(n, Sorter::scratch_limit);
    std::vector<uint64_t> values(n + 2 * GUARDS, GUARD);
    std::vector<uint64_t> scratch(capacity + 2 * GUARDS, GUARD);
    std::copy(input.begin(), input.end(), values.begin() + GUARDS);
    call_sort(mode, values.data() + GUARDS, n,
              capacity ? scratch.data() + GUARDS : nullptr, capacity);
    check(std::equal(expected.begin(), expected.end(), values.begin() + GUARDS),
          label + ": exact integer order or duplicate multiplicity changed");
    check_guards(values, n, label + " input");
    check_guards(scratch, capacity, label + " scratch");
}

std::vector<uint64_t> make_input(size_t n, unsigned pattern, uint64_t seed) {
    std::mt19937_64 random(seed);
    std::vector<uint64_t> input(n);
    for (size_t i = 0; i < n; ++i) {
        const uint64_t value = random();
        switch (pattern) {
        case 0: input[i] = value; break;
        case 1: input[i] = value & UINT32_MAX; break;
        case 2: input[i] = value & ((UINT64_C(1) << 42) - 1); break;
        case 3: input[i] = value & ((UINT64_C(1) << 21) - 1); break;
        case 4: input[i] = UINT64_C(0xffffffffffe00000) | (value & ((UINT64_C(1) << 21) - 1)); break;
        case 5: input[i] = (i % 2) ? UINT64_C(1) << 63 : 0; break;
        case 6: input[i] = UINT64_MAX; break;
        case 7: input[i] = i % 19; break;
        case 8: input[i] = i; break;
        case 9: input[i] = n - i; break;
        case 10: input[i] = i % 997 ? (value & 31) : value; break;
        default:
            // Almost all values repeatedly enter one MSD bucket, with one
            // far-away bit keeping the enclosing ranges much wider.
            input[i] = i % 4093 ? (value & 15) : UINT64_C(1) << 40;
            break;
        }
    }
    return input;
}

size_t check_contract() {
    size_t checks = 0;
    for (unsigned mode : {0u, 8u, 9u}) {
        call_sort(mode, nullptr, 0, nullptr, 0);
        uint64_t singleton = 17;
        call_sort(mode, &singleton, 1, nullptr, 0);
        check(singleton == 17, "Singleton changed without scratch");
        ++checks;
        for (size_t n : {size_t(2), size_t(128)}) {
            std::vector<uint64_t> scratch(n, GUARD);
            bool rejected = false;
            try { call_sort(mode, nullptr, n, scratch.data(), scratch.size()); }
            catch (const std::invalid_argument&) { rejected = true; }
            check(rejected, "Nonempty null input was not rejected");
            ++checks;
        }
        for (size_t n : {size_t(128), size_t(65536), size_t(65537)}) {
            auto input = make_input(n, 9, 17);
            const auto original = input;
            const size_t needed = std::min(n, Sorter::scratch_limit);
            std::vector<uint64_t> scratch(needed + 2 * GUARDS, GUARD);
            for (bool null_pointer : {false, true}) {
                bool rejected = false;
                try {
                    call_sort(mode, input.data(), n,
                              null_pointer ? nullptr : scratch.data() + GUARDS,
                              null_pointer ? needed : needed - 1);
                } catch (const std::invalid_argument&) { rejected = true; }
                check(rejected, "Insufficient/null scratch was not rejected");
                check(input == original, "Invalid scratch modified input before rejection");
                check(std::all_of(scratch.begin(), scratch.end(),
                                  [](uint64_t x) { return x == GUARD; }),
                      "Invalid scratch was written before rejection");
                ++checks;
            }
        }
    }
    return checks;
}
} // namespace

int main() {
    try {
        const size_t contract_checks = check_contract();
        size_t sort_checks = 0;
        for (size_t n : {size_t(0), size_t(1), size_t(127), size_t(128), size_t(129),
                         size_t(65535), size_t(65536), size_t(65537), size_t(1000000)}) {
            for (unsigned pattern = 0; pattern < 12; ++pattern) {
                const auto input = make_input(n, pattern, 20261121);
                for (unsigned mode : {0u, 9u}) {
                    verify(input, "n=" + std::to_string(n) + " pattern=" + std::to_string(pattern), mode);
                    ++sort_checks;
                }
            }
        }
        for (uint64_t seed : {0ULL, 42ULL, ~0ULL}) {
            for (unsigned pattern : {0u, 1u, 2u, 3u, 10u, 11u}) {
                for (unsigned mode : {8u, 9u}) {
                    verify(make_input(65539, pattern, seed), "multiple-seed distribution", mode);
                    ++sort_checks;
                }
            }
        }
        std::cout << "{\"passed\":true,\"sort_checks\":" << sort_checks
                  << ",\"contract_checks\":" << contract_checks << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "bounded_radix_sort_test: " << error.what() << '\n';
        return 1;
    }
}
