#include "Sketch.h"
#include "fastkmv.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double infinity = std::numeric_limits<double>::infinity();
void check(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

// Independent oracle: materialize the union, then count membership in its
// prefix. Unlimited prefix gives ordinary set Jaccard; K gives KMV Jaccard.
double reference(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b,
                 size_t limit = std::numeric_limits<size_t>::max()) {
    if (a.empty() || b.empty()) return 0.0;
    std::vector<uint64_t> joined;
    std::set_union(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(joined));
    if (joined.size() > limit) joined.resize(limit);
    size_t common = 0;
    for (uint64_t h : joined)
        common += std::binary_search(a.begin(), a.end(), h) &&
                  std::binary_search(b.begin(), b.end(), h);
    return static_cast<double>(common) / joined.size();
}

template<class Jaccard, class Distance>
void verify(double expected, Jaccard jaccard, Distance distance) {
    check(jaccard(0.0) == expected, "unfiltered Jaccard differs from independent oracle");
    std::vector<double> thresholds = {0.0, 0.05, 0.2, 0.5, 0.9, 1.0, expected};
    if (expected > 0.0) thresholds.push_back(std::nextafter(expected, 0.0));
    if (expected < 1.0) thresholds.push_back(std::nextafter(expected, infinity));
    for (double threshold : thresholds) {
        const double actual = jaccard(threshold);
        if (expected >= threshold)
            check(actual == expected, "Jaccard pruning lost a qualifying pair");
        else
            check(actual == 0.0 || actual == expected, "pruning returned a partial estimate");
    }
    const double full_distance = distance(infinity);
    std::vector<double> cutoffs = {0.0, 0.02, 0.05, 0.1, 1.0, 2.0, infinity};
    if (std::isfinite(full_distance)) {
        cutoffs.push_back(full_distance);
        cutoffs.push_back(std::nextafter(full_distance, infinity));
        if (full_distance > 0.0) cutoffs.push_back(std::nextafter(full_distance, 0.0));
    }
    for (double cutoff : cutoffs) {
        const double actual = distance(cutoff);
        if (full_distance <= cutoff)
            check(actual == full_distance, "distance pruning lost a qualifying boundary pair");
        else
            check(actual == full_distance || actual == infinity, "distance is a partial estimate");
        // The CLI uses a strict cutoff, including at exact floating boundaries.
        check((actual < cutoff) == (full_distance < cutoff), "cutoff output membership changed");
    }
    check(jaccard(0.0) == expected && distance(infinity) == full_distance,
          "thresholded query changed subsequent unfiltered results");
}

std::vector<uint64_t> sample(std::mt19937_64& rng, size_t size, size_t universe) {
    std::vector<uint64_t> pool(universe);
    for (size_t i = 0; i < universe; ++i) pool[i] = i;
    std::shuffle(pool.begin(), pool.end(), rng);
    pool.resize(size);
    std::sort(pool.begin(), pool.end());
    return pool;
}

void fastkmvCases() {
    // 7/25 * 25 rounds above 7: a raw ceil would wrongly require 8 matches.
    {
        std::vector<uint64_t> a, b;
        for (uint64_t i = 0; i < 25; ++i) {
            a.push_back(i);
            b.push_back(i < 7 ? i : i + 18);
        }
        auto left = Sketch::FastKMV::fromRegisters(25, 21, 42, a);
        auto right = Sketch::FastKMV::fromRegisters(25, 21, 42, b);
        verify(7.0 / 25.0, [&](double t) { return left.jaccard(right, t); },
                            [&](double t) { return left.distance(right, t); });
    }
    std::mt19937_64 rng(73);
    for (uint32_t k : {2u, 3u, 7u, 8u, 9u, 15u, 16u, 17u, 31u, 32u, 33u, 64u, 65u, 1024u}) {
        for (unsigned trial = 0; trial < 160; ++trial) {
            const uint32_t kb = (trial & 1) ? k * 2 : k;
            const size_t na = trial % 3 == 0 ? rng() % (k + 1) : k;
            const size_t nb = trial % 5 == 0 ? rng() % (kb + 1) : kb;
            const size_t universe = kb * (2 + trial % 8);
            auto a = sample(rng, na, universe), b = sample(rng, nb, universe);
            if (trial % 11 == 0) b = a;
            if (trial % 13 == 0) {
                for (auto& h : a) h += Sketch::FastKMV::KEY_MAX - universe;
                for (auto& h : b) h += Sketch::FastKMV::KEY_MAX - universe;
            }
            auto left = Sketch::FastKMV::fromRegisters(k, 21, 42, a);
            auto right = Sketch::FastKMV::fromRegisters(kb, 21, 42, b);
            if (b.size() > k) b.resize(k);
            const double expected = reference(a, b, k);
            verify(expected, [&](double t) { return left.jaccard(right, t); },
                             [&](double t) { return left.distance(right, t); });
            verify(expected, [&](double t) { return right.jaccard(left, t); },
                             [&](double t) { return right.distance(left, t); });
        }
    }
    // Positive overlap confined to the end: a high threshold must reject.
    std::vector<uint64_t> a, b;
    for (uint64_t i = 0; i < 1024; ++i) { a.push_back(i * 2); b.push_back(i * 2 + 1); }
    b.back() = a.back();
    auto left = Sketch::FastKMV::fromRegisters(4096, 21, 42, a);
    auto right = Sketch::FastKMV::fromRegisters(4096, 21, 42, b);
    check(left.jaccard(right) > 0.0 && left.jaccard(right, 0.9) == 0.0,
          "FastKMV threshold path did not reject");
}

void kssdCases() {
    std::mt19937_64 rng(991);
    // Small subspace keeps fixtures cheap; cover both native tuple widths.
    for (bool wide : {false, true}) {
        auto params = Sketch::getKssdParameters(wide ? 16 : 8, 3, 0);
        Sketch::Kssd left(params), right(params);
        for (size_t size : {0u, 1u, 3u, 7u, 8u, 9u, 15u, 16u, 17u, 31u, 32u, 33u, 64u, 65u, 1024u}) {
            for (unsigned trial = 0; trial < 100; ++trial) {
                const size_t nb = trial % 3 ? size : size / 8;
                const size_t universe = std::max(size_t(1), size * (2 + trial % 8));
                auto a = sample(rng, size, universe), b = sample(rng, nb, universe);
                if (trial % 11 == 0) b = a;
                const uint64_t offset = trial % 2 ? (wide ? UINT64_MAX : UINT32_MAX) - universe : 0;
                for (auto& h : a) h += offset;
                for (auto& h : b) h += offset;
                if (wide) { left.loadHashes64(a); right.loadHashes64(b); }
                else {
                    left.loadHashes({a.begin(), a.end()});
                    right.loadHashes({b.begin(), b.end()});
                }
                const double expected = reference(a, b);
                check(left.jaccard(&right) == expected, "legacy KSSD Jaccard overload changed");
                if (expected == 0.0)
                    check(left.distance(&right) == 1.0, "legacy KSSD zero-overlap distance changed");
                verify(expected, [&](double t) { return left.jaccard(&right, t); },
                                 [&](double t) { return left.distance(&right, t); });
                verify(expected, [&](double t) { return right.jaccard(&left, t); },
                                 [&](double t) { return right.distance(&left, t); });
            }
        }
        std::vector<uint64_t> a, b;
        for (uint64_t i = 0; i < 1024; ++i) { a.push_back(i * 2); b.push_back(i * 2 + 1); }
        b.back() = a.back();
        if (wide) { left.loadHashes64(a); right.loadHashes64(b); }
        else { left.loadHashes({a.begin(), a.end()}); right.loadHashes({b.begin(), b.end()}); }
        check(left.jaccard(&right) > 0.0 && left.jaccard(&right, 0.9) == 0.0,
              "KSSD threshold path did not reject");
    }
    // Different reductions must be aligned before deriving cardinality bounds.
    auto dense = Sketch::getKssdParameters(10, 5, 1);
    auto sparse = Sketch::getKssdParameters(10, 5, 2);
    Sketch::Kssd left(dense), right(sparse);
    std::string seq(60000, 'A');
    for (char& c : seq) c = "ACGT"[rng() & 3];
    left.update(seq.data(), seq.size());
    for (size_t i = 5000; i < seq.size(); i += 17) seq[i] = "ACGT"[rng() & 3];
    right.update(seq.data(), seq.size());
    auto projected = left.project(2);
    const auto& a = projected.getHashes();
    const auto& b = right.getHashes();
    const double expected = reference({a.begin(), a.end()}, {b.begin(), b.end()});
    verify(expected, [&](double t) { return left.jaccard(&right, t); },
                     [&](double t) { return left.distance(&right, t); });
}

template<class F> void rejects(F f) {
    bool rejected = false;
    try { f(); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "invalid threshold was accepted");
}
void invalidThresholds() {
    Sketch::FastKMV a(8), b(8);
    auto params = Sketch::getKssdParameters(8, 3, 0);
    Sketch::Kssd c(params), d(params);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (double t : {-0.1, 1.1, nan, infinity, -infinity}) {
        rejects([&] { a.jaccard(b, t); });
        rejects([&] { c.jaccard(&d, t); });
    }
    for (double t : {-0.1, nan, -infinity}) {
        rejects([&] { a.distance(b, t); });
        rejects([&] { c.distance(&d, t); });
    }
}
} // namespace

int main() {
    try {
        fastkmvCases();
        kssdCases();
        invalidThresholds();
        std::cout << "pairwise_threshold_test: PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "pairwise_threshold_test: FAIL: " << e.what() << '\n';
        return 1;
    }
}
