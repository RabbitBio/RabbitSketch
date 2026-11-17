#include "api/SketchOperations.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

using namespace Sketch;
using namespace Sketch::API;

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class Exception, class Function>
void rejects(Function function, const char* message) {
    try { function(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}
std::string dna(size_t length) {
    std::mt19937_64 random(20261117);
    std::string result(length, 'A');
    for (char& base : result) base = "ACGT"[random() % 4];
    return result;
}

void testNativeAndUnified() {
    auto sequence = dna(5000);
    FastKMV fast(256, 21);
    fast.update(sequence.data(), sequence.size());
    check(BuiltSketch::fromFastKMV(fast).cardinality() == fast.cardinality(),
          "unified FastKMV changed the native estimate");
    FracMinHash frac(4, 21);
    frac.update(sequence.data(), sequence.size());
    check(BuiltSketch::fromFracMinHash(frac).cardinality() == frac.cardinality(),
          "unified FracMinHash changed the native estimate");
    HyperLogLog hll(8, 21);
    hll.update(sequence.data());
    check(BuiltSketch::fromHyperLogLog(hll).cardinality() == hll.cardinality(),
          "unified HLL changed the native estimate");
    SetSketch set(8, 2.0, 5.0, 21);
    set.update(sequence.data(), sequence.size());
    check(BuiltSketch::fromSetSketch(set).cardinality() == set.cardinality(),
          "unified SetSketch changed the native estimate");
    auto legacy = std::make_unique<MinHash>(21, 256, 42);
    std::string mutable_sequence = sequence;
    legacy->update(mutable_sequence.data());
    const auto legacy_built = BuiltSketch::fromLegacyMinHash(std::move(legacy));
    rejects<std::logic_error>([&] { legacy_built.cardinality(); },
                             "legacy MinHash exposed a new unified entry point");

    for (Algorithm algorithm : {Algorithm::FastKMV,
         Algorithm::FracMinHash, Algorithm::HyperLogLog, Algorithm::SetSketch,
         Algorithm::KSSD, Algorithm::BinDash})
        check(capabilitiesFor(algorithm).supports(Capability::Cardinality),
              "cardinality capability missing");
    check(!capabilitiesFor(Algorithm::OrderMinHash).supports(Capability::Cardinality),
          "OrderMinHash advertised ordinary cardinality");
    OrderMinHash order(sequence);
    const auto built = BuiltSketch::fromOrderMinHash(std::move(order));
    rejects<std::logic_error>([&] { built.cardinality(); },
                             "OrderMinHash returned ordinary cardinality");
}

void testKssd() {
    const auto sequence = dna(12000);
    // Exercise both 64-bit (dr=1) and 32-bit (dr=2) tuples at k=20.
    for (int dr : {1, 2}) {
        auto params = std::make_shared<kssd_parameter_t>(10, 5, dr);
        Kssd sketch(params);
        check(sketch.cardinality() == 0.0, "empty KSSD estimate is nonzero");
        sketch.update(sequence.data(), sequence.size());
        const double before = sketch.cardinality();
        check(before > 0.0 && !sketch.isSealed(), "KSSD query sealed or lost a sample");
        sketch.update(sequence.data(), sequence.size());
        check(sketch.cardinality() == before, "KSSD counted repeated input");
        sketch.finalize();
        const auto retained = sketch.getHashes().size() + sketch.getHashes64().size();
        check(before == retained / sketch.admissionProbability(),
              "KSSD cardinality disagrees with its sampling probability");
        check(sketch.merge(sketch).cardinality() == before,
              "KSSD self merge changed cardinality");
        auto restored = std::make_unique<Kssd>(params);
        if (dr == 1) restored->loadHashes64(sketch.getHashes64());
        else restored->loadHashes(sketch.getHashes());
        check(BuiltSketch::fromKssd(std::move(restored)).cardinality() == before,
              "restored/unified KSSD estimate changed");
        if (dr == 1) {
            auto projected = sketch.project(2);
            Kssd direct(std::make_shared<kssd_parameter_t>(10, 5, 2));
            direct.update(sequence.data(), sequence.size());
            check(projected.cardinality() == direct.cardinality(),
                  "KSSD projection used the original admission probability");
        }
    }
}

void testBinDash() {
    const double infinity = std::numeric_limits<double>::infinity();
    const std::vector<uint64_t> packed(16, 0);
    auto empty = BinDash::fromPacked(1, 21, 16, 42, 0, packed);
    auto partial = BinDash::fromPacked(1, 21, 16, 42, 32, packed);
    auto saturated = BinDash::fromPacked(1, 21, 16, 42, 64, packed);
    check(empty.cardinality() == 0.0, "empty BinDash estimate is nonzero");
    check(std::abs(partial.cardinality() - 64.0 * std::log(2.0)) < 1e-12,
          "BinDash lost occupancy metadata on restore");
    check(saturated.cardinality() == infinity &&
          BuiltSketch::fromBinDash(saturated).cardinality() == infinity,
          "saturated BinDash reported a finite cardinality");
    BinDash compact(64, 21, 16, 42, false);
    const auto sequence = dna(1000);
    compact.update(sequence.data(), sequence.size());
    const double before = compact.cardinality();
    compact.finalize();
    check(std::isfinite(before) && before > 0 && compact.cardinality() == before &&
          BuiltSketch::fromBinDash(compact).cardinality() == before,
          "compact BinDash finalization changed occupancy cardinality");
}

void testUnsupportedProbMinHash() {
    check(!capabilitiesFor(Algorithm::ProbMinHash).supports(Capability::Cardinality),
          "ProbMinHash advertised cardinality");
    ProbMinHash4 sketch(128, 21);
    sketch.addHash(1);
    const auto built = BuiltSketch::fromProbMinHash(std::move(sketch));
    rejects<std::logic_error>([&] { built.cardinality(); },
                             "ProbMinHash exposed a cardinality estimator");
}
} // namespace

int main() {
    try {
        testNativeAndUnified();
        testKssd();
        testBinDash();
        testUnsupportedProbMinHash();
        std::cout << "cardinality tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
