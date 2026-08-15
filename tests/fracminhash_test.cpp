#include "estimators/UnifiedQuery.h"
#include "fracminhash.h"
#include "rank/CanonicalKmer.h"
#include "rank/RankStream.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <iterator>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Exception, typename Function>
void expectThrows(Function function, const std::string& message) {
    try {
        function();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error(message);
}

std::string dna(size_t length, uint64_t seed) {
    static constexpr char bases[] = {'A', 'C', 'G', 'T'};
    std::string result(length, 'A');
    uint64_t state = seed;
    for (size_t index = 0; index < length; ++index) {
        state = state * UINT64_C(6364136223846793005) +
            UINT64_C(1442695040888963407);
        result[index] = bases[state >> 62];
    }
    for (size_t index = 997; index < length; index += 4093)
        result[index] = 'N';
    return result;
}

std::vector<uint64_t> reference(
    const std::vector<std::string>& sequences, uint32_t scaled) {
    std::set<uint64_t> values;
    const uint64_t threshold = scaled == 1
        ? UINT64_MAX : UINT64_MAX / scaled;
    const Sketch::Rank::RankStream ranks(42);
    for (const std::string& sequence : sequences) {
        Sketch::Rank::CanonicalKmerIterator iterator(
            sequence.data(), sequence.size(), 21);
        uint64_t canonical = 0;
        while (iterator.next(canonical)) {
            const uint64_t fingerprint = ranks.fingerprint(canonical);
            if (fingerprint <= threshold) values.insert(fingerprint);
        }
    }
    return std::vector<uint64_t>(values.begin(), values.end());
}

Sketch::FracMinHash build(const std::vector<std::string>& sequences,
                          uint32_t scaled) {
    Sketch::FracMinHash result(scaled, 21, 42);
    for (const std::string& sequence : sequences)
        result.update(sequence.data(), sequence.size());
    result.finalize();
    return result;
}

void testOptimizedUpdateMatchesReference() {
    const std::vector<std::string> sequences = {
        dna(30000, 17), dna(17000, 29), std::string(5000, 'A')};
    for (uint32_t scaled : {1u, 7u, 64u, 65u, 1000u}) {
        const Sketch::FracMinHash sketch = build(sequences, scaled);
        check(sketch.hashes() == reference(sequences, scaled),
              "optimized FracMinHash update changed retained hashes");
        check(sketch.isSealed() &&
                  sketch.memoryBytes() >= sketch.size() * sizeof(uint64_t),
              "FracMinHash finalize lifecycle is incomplete");
    }

    const size_t lengths[] = {
        0, 1, 20, 21, 22, 27, 28, 29, 31, 32, 35, 36, 37};
    for (size_t length : lengths) {
        const std::vector<std::string> tail = {dna(length, length + 101)};
        for (uint32_t scaled : {2u, 65u})
            check(build(tail, scaled).hashes() == reference(tail, scaled),
                  "FracMinHash tail batch changed retained hashes");
    }
}

void testProjectionMergeAndUnifiedQuery() {
    const Sketch::FracMinHash left = build({dna(120000, 41)}, 2);
    const Sketch::FracMinHash right = build({dna(120000, 43)}, 5);
    const auto projected = left.project(5);
    check(projected.hashes().size() == left.retainedAt(5) &&
              projected.getScaled() == 5,
          "FracMinHash projection retained the wrong threshold prefix");

    const auto merged = left.merge(right);
    std::vector<uint64_t> expected;
    expected.reserve(projected.size() + right.size());
    std::set_union(projected.hashes().begin(), projected.hashes().end(),
                   right.hashes().begin(), right.hashes().end(),
                   std::back_inserter(expected));
    check(merged.getScaled() == 5 && merged.hashes() == expected,
          "FracMinHash direct merge differs from projected set union");

    const double native = left.jaccard(right);
    const Sketch::Query::Result unified = Sketch::Query::query(left, right);
    check(unified.estimate_available &&
              std::abs(unified.jaccard - native) < 1e-15 &&
              unified.effective_samples ==
                  std::min(left.retainedAt(5), right.retainedAt(5)) &&
              std::abs(unified.left_cardinality - left.cardinalityAt(5)) <
                  1e-9,
          "no-copy unified FracMinHash query changed its estimand");
}

void testValidationAndExternalHashes() {
    Sketch::FracMinHash sketch(4, 21, 42);
    sketch.addHash(0);
    sketch.addHash(0);
    sketch.addHash(UINT64_MAX / 4);
    sketch.addHash(UINT64_MAX);
    sketch.finalize();
    check(sketch.hashes() ==
              std::vector<uint64_t>({0, UINT64_MAX / 4}),
          "FracMinHash addHash threshold/deduplication is incorrect");
    expectThrows<std::logic_error>(
        [&]() { sketch.addHash(1); },
        "FracMinHash accepted an update after finalize");
    expectThrows<std::invalid_argument>(
        [&]() { (void)sketch.project(2); },
        "FracMinHash projected to an unavailable denser scale");
    expectThrows<std::invalid_argument>(
        [&]() {
            (void)Sketch::FracMinHash::fromHashes(
                4, 21, 42, {2, 1});
        },
        "FracMinHash accepted an unsorted restored state");

    Sketch::FracMinHash wrong_seed(4, 21, 43);
    expectThrows<std::invalid_argument>(
        [&]() { (void)sketch.jaccard(wrong_seed); },
        "FracMinHash compared incompatible rank spaces");
}

} // namespace

int main() {
    try {
        testOptimizedUpdateMatchesReference();
        testProjectionMergeAndUnifiedQuery();
        testValidationAndExternalHashes();
        std::cout << "fracminhash_test: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "fracminhash_test: FAIL: " << error.what() << '\n';
        return 1;
    }
}
