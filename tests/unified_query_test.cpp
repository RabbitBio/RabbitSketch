#include "api/SketchConfig.h"
#include "estimators/UnifiedQuery.h"
#include "fastkmv.h"
#include "fracminhash.h"

#include <cstdint>
#include <exception>
#include <iostream>
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

std::string deterministicDna(size_t length, uint64_t seed) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::string sequence(length, 'A');
    uint64_t state = seed;
    for (size_t index = 0; index < length; ++index) {
        state = state * UINT64_C(6364136223846793005) +
            UINT64_C(1442695040888963407);
        sequence[index] = bases[(state >> 62) & 3u];
    }
    return sequence;
}

void testSameBackendPlanning() {
    const std::string sequence = deterministicDna(5000, 17);
    Sketch::FastKMV large(256, 21, 42);
    Sketch::FastKMV small(64, 21, 42);
    large.update(sequence.data(), sequence.size());
    small.update(sequence.data(), sequence.size());

    const Sketch::Query::Result native = Sketch::Query::query(large, small);
    check(native.estimate_available && native.value == 1.0 &&
              native.plan.method == Sketch::Query::Method::FastKmvNative &&
              native.plan.common_resolution == 64 &&
              native.plan.exact_projection,
          "FastKMV did not use exact common-resolution projection");

}

void testPlannerRejections() {
    const std::string sequence = deterministicDna(2000, 29);
    Sketch::FastKMV fast(64, 21, 42);
    Sketch::FastKMV wrong_seed(64, 21, 43);
    Sketch::FracMinHash frac(2, 21, 42);
    fast.update(sequence.data(), sequence.size());
    wrong_seed.update(sequence.data(), sequence.size());
    frac.update(sequence.data(), sequence.size());

    const Sketch::Query::Result incompatible =
        Sketch::Query::query(fast, wrong_seed);
    check(!incompatible.estimate_available &&
              incompatible.plan.support == Sketch::Query::Support::Incompatible,
          "planner accepted sketches built with different seeds");

    const Sketch::Query::Plan cross = Sketch::Query::planQuery(
        Sketch::Query::describe(fast), Sketch::Query::describe(frac));
    check(!cross.executable() &&
              cross.support == Sketch::Query::Support::Unsupported &&
              cross.reason == "cross-backend comparisons are not supported",
          "planner accepted a cross-backend query");

    Sketch::Query::Request invalid;
    invalid.confidence_level = 1.0;
    expectThrows<std::invalid_argument>(
        [&]() {
            (void)Sketch::Query::planQuery(
                Sketch::Query::describe(fast),
                Sketch::Query::describe(fast), invalid);
        },
        "planner accepted confidence_level=1");
}

void testPublicCapabilityContract() {
    Sketch::API::SketchConfig config;
    config.validate();
    check(Sketch::API::capabilitiesFor(Sketch::API::Algorithm::FastKMV)
              .supports(Sketch::API::Capability::Containment) &&
              Sketch::API::capabilitiesFor(
                  Sketch::API::Algorithm::ProbMinHash)
              .supports(Sketch::API::Capability::WeightedJaccard) &&
              Sketch::API::capabilitiesFor(
                  Sketch::API::Algorithm::OrderMinHash)
              .supports(Sketch::API::Capability::OrderDistance),
          "public capability table omitted a production operation");
    config.kmer_size = 0;
    expectThrows<std::invalid_argument>(
        [&]() { config.validate(); },
        "SketchConfig accepted k=0");
}

} // namespace

int main() {
    try {
        testSameBackendPlanning();
        testPlannerRejections();
        testPublicCapabilityContract();
        std::cout << "unified_query_test: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "unified_query_test: FAIL: " << error.what() << '\n';
        return 1;
    }
}
