#include "probmh.h"
#include "ProbMinHashBuild.h"
#include "rank/RankStream.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

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
    std::sort(ak.begin(), ak.end()); std::sort(bk.begin(), bk.end());
    require(ak == bk, "Sorted inverted-index keys changed");
}

static size_t check_batch_hash() {
    using Sketch::Runtime::hashFmix64Batch;
    constexpr uint64_t guard = UINT64_C(0x73ac9d5e6081fb24);
    std::mt19937_64 random(20261121);
    size_t checks = 0;
    for (size_t count : {size_t(0), size_t(1), size_t(7), size_t(8), size_t(9),
                         size_t(31), size_t(32), size_t(255), size_t(256),
                         size_t(257), size_t(1003)}) {
        for (size_t offset : {size_t(1), size_t(2), size_t(3)}) {
            for (uint64_t seed : {0ULL, 42ULL, ~0ULL}) {
                // Offsets exercise valid uint64_t alignment without assuming
                // 16/32/64-byte SIMD alignment, including every short tail.
                std::vector<uint64_t> input(count + offset + 4, guard);
                std::vector<uint64_t> output(input.size(), guard);
                for (size_t i = 0; i < count; ++i)
                    input[offset + i] = i % 13 ? random() : (i % 2 ? 0 : UINT64_MAX);
                const auto original = input;
                auto expected = output;
                for (size_t i = 0; i < count; ++i)
                    expected[offset + i] = Sketch::Rank::RankStream::fmix64(input[offset + i], seed);
                hashFmix64Batch(input.data() + offset, count, seed, output.data() + offset);
                require(output == expected, "Batch hash differs from scalar truth or wrote outside output range");
                require(input == original, "Separate-output batch hash changed its input");
                hashFmix64Batch(input.data() + offset, count, seed, input.data() + offset);
                auto expected_in_place = original;
                std::copy(expected.begin() + offset, expected.begin() + offset + count,
                          expected_in_place.begin() + offset);
                require(input == expected_in_place, "Exact in-place batch hash changed values or guards");
                ++checks;
            }
        }
    }
    std::array<uint64_t, 8> untouched{};
    untouched.fill(guard);
    hashFmix64Batch(nullptr, 0, 42, nullptr);
    hashFmix64Batch(nullptr, untouched.size(), 42, untouched.data());
    hashFmix64Batch(untouched.data(), untouched.size(), 42, nullptr);
    require(std::all_of(untouched.begin(), untouched.end(),
                        [](uint64_t value) { return value == guard; }),
            "Null-input/output batch hash changed caller storage");
    const char* requested = std::getenv("RABBITSKETCH_SIMD");
    if (requested && std::string(requested) == "scalar")
        require(Sketch::Runtime::runtimeInfo().selected_rank_path == "scalar",
                "Batch hash test did not honor scalar runtime selection");
    return checks + 1;
}

static std::vector<uint64_t> canonical_occurrences(const std::vector<std::string>& records) {
    constexpr unsigned k = 21;
    constexpr uint64_t mask = (UINT64_C(1) << (2 * k)) - 1;
    std::vector<uint64_t> result;
    for (const auto& sequence : records) {
        uint64_t forward = 0, reverse = 0;
        unsigned valid = 0;
        for (char base : sequence) {
            const std::string alphabet = "ACGT";
            const size_t value = alphabet.find(base);
            if (value == std::string::npos) { forward = reverse = valid = 0; continue; }
            forward = ((forward << 2) | value) & mask;
            reverse = (reverse >> 2) | (uint64_t(3 - value) << (2 * (k - 1)));
            if (++valid >= k) result.push_back(std::min(forward, reverse));
        }
    }
    return result;
}

static void feed_counts(Sketch::ProbMinHash4& sketch, const std::vector<uint64_t>& sorted,
                        bool batch) {
    if (sorted.empty()) return;
    const double inverse_total = 1.0 / static_cast<double>(sorted.size());
    std::array<uint64_t, 256> keys{};
    std::array<double, 256> weights{};
    size_t used = 0;
    for (size_t i = 0; i < sorted.size();) {
        size_t j = i + 1;
        while (j < sorted.size() && sorted[j] == sorted[i]) ++j;
        const double weight = static_cast<double>(j - i) * inverse_total;
        if (batch) {
            keys[used] = sorted[i]; weights[used++] = weight;
            if (used == keys.size()) {
                sketch.addHashBatch(keys.data(), weights.data(), used);
                used = 0;
            }
        } else sketch.addHash(sorted[i], weight);
        i = j;
    }
    if (used) sketch.addHashBatch(keys.data(), weights.data(), used);
}

static size_t check_counting_pipeline() {
    std::mt19937_64 random(20261121);
    std::string sequence(120037, 'A'), repeat_unit(701, 'A');
    for (size_t i = 0; i < sequence.size(); ++i)
        sequence[i] = i % 509 ? "ACGT"[random() % 4] : 'N';
    for (char& base : repeat_unit) base = "ACGT"[random() % 4];
    std::string repetitive;
    for (unsigned i = 0; i < 23; ++i) repetitive += repeat_unit;
    const std::vector<std::vector<std::string>> inputs = {
        {}, {std::string(20, 'A'), std::string(100, 'N')},
        {std::string(1301, 'A')},
        {sequence, repetitive, repeat_unit, std::string(1003, 'A')},
    };
    size_t checks = 0;
    Sketch::ProbMinHashBuild::Workspace workspace;
    for (const auto& records : inputs) {
        const auto original = canonical_occurrences(records);
        auto expected = original;
        std::sort(expected.begin(), expected.end());
        workspace.begin(original.size());
        workspace.occurrences.insert(workspace.occurrences.end(), original.begin(), original.end());
        workspace.sort();
        require(workspace.occurrences == expected, "DNA occurrence sorting altered frequency counts");
        for (uint32_t m : {1024u, 4096u, 16384u}) {
            Sketch::ProbMinHash4 reference(m, 21, 42, 0), sorted_batch(m, 21, 42, 0), actual(m, 21, 42, 0);
            feed_counts(reference, expected, false);
            feed_counts(sorted_batch, expected, true);
            feed_counts(actual, workspace.occurrences, true);
            equal(reference, sorted_batch);
            equal(reference, actual);
            reference.finalize(); actual.finalize(); equal(reference, actual);
            ++checks;
        }
    }
    workspace.release();
    return checks;
}

static size_t check_workspace_capacity() {
    using Workspace = Sketch::ProbMinHashBuild::Workspace;
    const size_t limit = Workspace::retained_limit;
    Workspace workspace;
    size_t checks = 0;

    // Alternating ordinary genomes must clear the previous sample without
    // throwing away and faulting in the same counting pages again.
    const size_t ordinary = limit * 3 / 4;
    workspace.begin(ordinary);
    const size_t retained = workspace.occurrences.capacity();
    const uint64_t* storage = workspace.occurrences.data();
    require(retained >= ordinary && retained <= limit,
            "Ordinary input exceeded the retained workspace bound");
    for (size_t repeat = 0; repeat < 16; ++repeat) {
        for (size_t hint : {size_t(0), size_t(1), limit / 8, ordinary}) {
            workspace.occurrences.push_back(UINT64_C(123));
            workspace.begin(hint);
            require(workspace.occurrences.empty(), "Previous sample survived workspace.begin");
            require(workspace.occurrences.capacity() == retained &&
                    workspace.occurrences.data() == storage,
                    "Alternating ordinary genomes discarded reusable counting pages");
            ++checks;
        }
    }

    // A gradual increase should need only a handful of allocations, rather
    // than one allocation for almost every slightly larger input file.
    workspace.release();
    workspace.begin(limit / 4);
    size_t capacity = workspace.occurrences.capacity();
    size_t growths = 0;
    const size_t step = std::max(size_t(1), limit / 128);
    for (size_t hint = limit / 4 + step; hint <= ordinary; hint += step) {
        workspace.begin(hint);
        const size_t next = workspace.occurrences.capacity();
        require(next >= hint && next <= limit,
                "Gradual growth failed to reserve input capacity within the retention bound");
        growths += next != capacity;
        capacity = next;
    }
    require(growths <= 8, "Gradual input growth reallocated counting pages too often");
    ++checks;

    // The retention limit is a cache policy, not an input-size limit. A large
    // genome must still fit, but its exceptional capacity must not stay pinned
    // after a return to ordinary input sizes.
    const size_t outlier = limit + limit / 4 + 1;
    workspace.begin(outlier);
    require(workspace.occurrences.capacity() >= outlier,
            "Outlier genome capacity was truncated to the retention limit");
    workspace.occurrences.push_back(UINT64_MAX);
    workspace.begin(limit / 4);
    require(workspace.occurrences.empty(), "Outlier sample survived workspace.begin");
    require(workspace.occurrences.capacity() >= limit / 4 &&
            workspace.occurrences.capacity() <= limit,
            "Outlier counting capacity remained pinned for an ordinary genome");
    ++checks;

    // A compressed file can have an underestimated hint and grow the vector
    // past the retention limit, even though its actual occurrence count fits.
    // Simulate that allocation history without depending on vector's growth
    // factor, then retain enough learned capacity for the next low-hint file.
    workspace.release();
    workspace.occurrences.reserve(limit + limit / 3);
    workspace.occurrences.resize(ordinary, UINT64_C(23));
    workspace.begin(0);
    require(workspace.occurrences.empty(), "Low-hint sample survived workspace.begin");
    require(workspace.occurrences.capacity() >= ordinary &&
            workspace.occurrences.capacity() <= limit,
            "Underestimated hint discarded the previous ordinary sample's capacity");
    const size_t learned = workspace.occurrences.capacity();
    const uint64_t* learned_storage = workspace.occurrences.data();
    for (size_t hint : {size_t(0), size_t(1), limit / 8}) {
        workspace.occurrences.push_back(UINT64_C(31));
        workspace.begin(hint);
        require(workspace.occurrences.empty() &&
                workspace.occurrences.capacity() == learned &&
                workspace.occurrences.data() == learned_storage,
                "Repeated low hints discarded learned ordinary counting capacity");
        ++checks;
    }

    for (size_t hint : {workspace.occurrences.max_size() + 1,
                        std::numeric_limits<size_t>::max()}) {
        bool threw = false;
        try { workspace.begin(hint); }
        catch (const std::length_error&) { threw = true; }
        require(threw, "Impossible reserve hint was not rejected with length_error");
        ++checks;
    }

    // Releasing construction state must include the bounded sorting scratch,
    // and the same worker must remain usable after release or an invalid hint.
    workspace.begin(128);
    for (uint64_t i = 0; i < 128; ++i) workspace.occurrences.push_back(127 - i);
    workspace.sort();
    require(workspace.scratch && workspace.scratch_capacity > 0,
            "Release test did not exercise allocated sorting scratch");
    workspace.release();
    require(workspace.occurrences.empty() && workspace.occurrences.capacity() == 0 &&
            !workspace.scratch && workspace.scratch_capacity == 0,
            "Release retained construction buffers or stale scratch capacity");
    workspace.begin(1);
    workspace.occurrences.push_back(17);
    workspace.begin(1);
    require(workspace.occurrences.empty(), "Released workspace could not be reused");
    workspace.release();
    ++checks;
    return checks;
}

int main(int argc, char** argv) {
    const size_t hash_checks = check_batch_hash();
    if (argc == 2 && std::string(argv[1]) == "--hash-only") {
        std::cout << "{\"passed\":true,\"batch_hash_checks\":" << hash_checks << "}\n";
        return 0;
    }
    std::mt19937_64 random(20261119);
    Sketch::ProbMinHashBuild::Workspace workspace;
    const size_t capacity_checks = check_workspace_capacity();
    const size_t pipeline_checks = check_counting_pipeline();
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
              << ",\"capacity_checks\":" << capacity_checks
              << ",\"batch_hash_checks\":" << hash_checks
              << ",\"counting_pipeline_checks\":" << pipeline_checks
              << ",\"batch_checks\":" << batch_checks
              << ",\"threshold_checks\":" << threshold_checks << "}\n";
}
