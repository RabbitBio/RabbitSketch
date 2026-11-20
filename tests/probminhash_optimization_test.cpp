#include "probmh.h"
#include "ProbMinHashBuild.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

static void require(bool okay, const char* message) {
    if (!okay) throw std::runtime_error(message);
}

static void equal(const Sketch::ProbMinHash4& a, const Sketch::ProbMinHash4& b) {
    require(a.getM() == b.getM(), "Register count changed");
    require(std::memcmp(a.getRegisters(), b.getRegisters(), a.getM() * sizeof(double)) == 0,
            "Register bits changed");
    require(std::memcmp(a.getWinners(), b.getWinners(), a.getM() * sizeof(uint64_t)) == 0,
            "Winner bits changed");
    const double aw = a.total_weight(), bw = b.total_weight();
    require(std::memcmp(&aw, &bw, sizeof(double)) == 0, "Weight accumulation changed");
    require(a.weightSemantics() == b.weightSemantics(), "Weight metadata changed");
    std::vector<uint64_t> ak, bk;
    a.getWinnerIndexKeys(ak); b.getWinnerIndexKeys(bk);
    require(ak == bk, "Index keys changed");
}

int main() {
    std::mt19937_64 random(20261119);
    Sketch::ProbMinHashBuild::Workspace workspace;
    size_t sort_checks = 0, batch_checks = 0, threshold_checks = 0;
    for (size_t n : {size_t(0), size_t(1), size_t(8191), size_t(8192), size_t(100001), size_t(1048579)}) {
        for (int pattern = 0; pattern < 7; ++pattern) {
            workspace.begin(n);
            for (size_t i = 0; i < n; ++i) {
                uint64_t value = pattern == 0 ? 0 : pattern == 1 ? UINT64_MAX :
                    pattern == 2 ? i % 19 : pattern == 3 ? n - i : random();
                if (pattern == 5) value &= (UINT64_C(1) << 42) - 1;
                if (pattern == 6) value &= (UINT64_C(1) << 21) - 1;
                workspace.occurrences.push_back(value);
            }
            auto expected = workspace.occurrences;
            std::sort(expected.begin(), expected.end());
            workspace.sort();
            require(workspace.occurrences == expected, "Radix sorting changed integers or counts");
            require(workspace.scratch_capacity <= workspace.scratch_limit, "Sort scratch exceeded its bound");
            ++sort_checks;
        }
    }
    workspace.release();
    require(workspace.occurrences.capacity() == 0 && !workspace.scratch,
            "Construction workspace was not released");
    for (uint32_t m : {2u, 7u, 1024u, 4096u, 16384u}) {
        for (uint64_t seed : {0ULL, 42ULL, ~0ULL}) {
            Sketch::ProbMinHash4 a(m, 21, seed, 0), b(m, 21, seed, 0);
            for (size_t length : {size_t(0), size_t(1), size_t(7), size_t(8), size_t(9), size_t(256), size_t(65537)}) {
                std::vector<uint64_t> hashes(length);
                std::vector<double> weights(length);
                for (size_t i = 0; i < length; ++i) {
                    hashes[i] = i % 7 ? random() : uint64_t(i % 11);
                    weights[i] = i % 19 ? 1.0 / 1234567.0 : static_cast<double>(i % 17 + 1) / 1234567.0;
                    if (i % 251 == 0) weights[i] = 0.0;
                    if (i % 257 == 0) weights[i] = -1.0;
                    if (i % 263 == 0) weights[i] = std::numeric_limits<double>::quiet_NaN();
                    a.addHash(hashes[i], weights[i]);
                }
                b.addHashBatch(hashes.data(), weights.data(), length);
                equal(a, b); ++batch_checks;
            }
            const uint64_t extremes[] = {0, 1, UINT64_MAX, 17, 99};
            const double extreme_weights[] = {1.0, 1e-300, 1e300, 2.0, std::numeric_limits<double>::infinity()};
            for (size_t i = 0; i < 5; ++i) a.addHash(extremes[i], extreme_weights[i]);
            b.addHashBatch(extremes, extreme_weights, 5);
            equal(a, b); ++batch_checks;
            a.finalize(); b.finalize(); equal(a, b);
            bool threw = false;
            try { b.addHashBatch(extremes, extreme_weights, 5); }
            catch (const std::logic_error&) { threw = true; }
            require(threw, "Batch update bypassed sealing");
        }
        for (size_t overlap : {size_t(0), size_t(m / 5), size_t(m / 2), size_t(m)}) {
            std::vector<uint64_t> a(m), b(m);
            size_t expected = 0;
            for (size_t i = 0; i < m; ++i) {
                a[i] = i % 31 ? random() : 0;
                b[i] = i < overlap ? a[i] : random();
                expected += a[i] != 0 && a[i] == b[i];
            }
            for (size_t minimum : {size_t(0), size_t(1), expected, expected + 1, size_t(m), size_t(m + 1)}) {
                std::vector<uint8_t> ab(a.begin(), a.end()), bb(b.begin(), b.end());
                const bool allowed = Sketch::ProbMinHashBuild::winnerByteFilterAllows(ab.data(), bb.data(), m, minimum);
                require(allowed || expected < minimum, "Byte filter removed a true winner match");
                // Worst-case collisions must pass the filter but still be
                // rejected by the full winner comparison when appropriate.
                std::vector<uint8_t> collisions(m, 0);
                require(Sketch::ProbMinHashBuild::winnerByteFilterAllows(
                            collisions.data(), collisions.data(), m, minimum) == (minimum <= m),
                        "Byte filter mishandled collisions or an impossible cutoff");
                size_t common = 0;
                const bool passed = Sketch::ProbMinHashBuild::countWinnersAtLeast(a.data(), b.data(), m, minimum, common);
                require(passed == (expected >= minimum), "Threshold pruning lost or added a pair");
                if (passed) require(common == expected, "Accepted pair count changed");
                ++threshold_checks;
            }
        }
    }
    std::cout << "{\"passed\":true,\"sort_checks\":" << sort_checks
              << ",\"batch_checks\":" << batch_checks
              << ",\"threshold_checks\":" << threshold_checks << "}\n";
}
