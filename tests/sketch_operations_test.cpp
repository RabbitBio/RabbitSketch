#include "api/SketchOperations.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>

using namespace Sketch;
using namespace Sketch::API;

namespace {
void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
template<class Exception, class Function>
void rejects(Function function, const std::string& message) {
    try { function(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}
std::string dna(size_t n, uint64_t seed) {
    std::mt19937_64 random(seed);
    std::string result(n, 'A');
    for (auto& base : result) base = "ACGT"[random() % 4];
    return result;
}
BuiltSketch build(const SketchConfig& config, const std::vector<std::string>& records) {
    MultiSketchBuilder builder({config});
    for (const auto& sequence : records) {
        IO::FastxRecord record;
        record.name = "sample";
        record.sequence = sequence;
        builder.update(record);
    }
    return builder.finalize("sample").front().sketch;
}

SketchConfig configuration(Algorithm algorithm) {
    SketchConfig c;
    c.algorithm = algorithm;
    c.kmer_size = 15;
    c.resolution = 128;
    switch (algorithm) {
        case Algorithm::FracMinHash: c.sampling = SamplingMode::Scaled; c.scaled = 4; break;
        case Algorithm::HyperLogLog:
        case Algorithm::SetSketch: c.sampling = SamplingMode::RegisterPrecision; c.resolution = 8; break;
        case Algorithm::BinDash: c.sampling = SamplingMode::AlgorithmNative; c.resolution = 256; break;
        case Algorithm::ProbMinHash: c.sampling = SamplingMode::AlgorithmNative; break;
        case Algorithm::KSSD:
            c.sampling = SamplingMode::KssdReduction;
            c.kmer_size = 20; c.kssd_half_subk = 5; c.resolution = 1; break;
        default: break;
    }
    return c;
}

void testSevenBackends() {
    const std::string a = dna(16000, 51);
    const std::string b = a.substr(0, 12000) + dna(4000, 52);
    const std::string c = dna(16000, 53);
    for (const auto algorithm : {Algorithm::FastKMV, Algorithm::FracMinHash,
         Algorithm::HyperLogLog, Algorithm::SetSketch, Algorithm::ProbMinHash,
         Algorithm::BinDash, Algorithm::KSSD}) {
        auto config = configuration(algorithm);
        const std::string name = algorithmName(algorithm);
        const auto left = build(config, {a}), right = build(config, {b});
        const auto united = build(config, {a, b});
        const auto merged = left.merge(right);
        check(merged.storeHashes() == united.storeHashes(), name + " merge differs from direct union build");
        check(left.storeHashes() == build(config, {a}).storeHashes(), name + " merge mutated source");
        check(left.merge(left).storeHashes() == left.storeHashes(), name + " merge is not idempotent");

        const auto caps = capabilitiesFor(algorithm);
        for (auto capability : {Capability::Merge, Capability::Project,
                                Capability::StoreHashes, Capability::IndexDict})
            check(caps.supports(capability), name + " common API capability missing");
        const bool cardinality = algorithm != Algorithm::ProbMinHash;
        check(caps.supports(Capability::Cardinality) == cardinality, name + " cardinality capability mismatch");
        if (cardinality) check(left.cardinality() > 1000, name + " cardinality estimate failed");
        else rejects<std::logic_error>([&] { left.cardinality(); }, name + " exposed cardinality");

        uint32_t target = 64;
        if (algorithm == Algorithm::HyperLogLog || algorithm == Algorithm::SetSketch) target = 6;
        if (algorithm == Algorithm::BinDash) target = 128;
        if (algorithm == Algorithm::FracMinHash) target = 8;
        if (algorithm == Algorithm::KSSD) target = 2;
        const auto projected = left.project(target);
        const auto projected_right = right.project(target);
        check(projected.merge(projected_right).storeHashes() == merged.project(target).storeHashes(),
              name + " merge/project do not commute");
        if (algorithm == Algorithm::ProbMinHash) {
            const auto values = left.storeHashes();
            check(projected.storeHashes() == std::vector<uint64_t>(values.begin(), values.begin() + target),
                  "ProbMinHash projection lost source race prefix");
        } else {
            auto smaller = config;
            if (algorithm == Algorithm::FracMinHash) smaller.scaled = target;
            else smaller.resolution = target;
            check(projected.storeHashes() == build(smaller, {a}).storeHashes(),
                  name + " projection differs from direct construction at target resolution");
        }

        const std::vector<BuiltSketch> sketches{left, right, build(config, {c}), left};
        const auto index = BuiltSketch::indexDict(sketches, 4);
        check(index.size() == 4 && index.keyCount() > 0, name + " index empty");
        const bool approximate = algorithm == Algorithm::HyperLogLog || algorithm == Algorithm::SetSketch;
        check(index.approximateCandidates() == approximate, name + " wrong index recall contract");
        const auto matches = index.pairs(0.05, 4);
        std::set<std::pair<uint32_t, uint32_t>> found;
        for (const auto& match : matches) {
            check(match.left < match.right, name + " duplicate or symmetric index pair");
            check(found.emplace(match.left, match.right).second, name + " duplicate result");
            check(std::abs(match.estimate.value - sketches[match.left].query(sketches[match.right]).value) < 1e-12,
                  name + " index estimate differs from native query");
        }
        check(found.count({0, 3}) == 1, name + " index missed identical nonempty sketches");
        if (!index.approximateCandidates()) {
            for (uint32_t i = 0; i < sketches.size(); ++i)
                for (uint32_t j = i + 1; j < sketches.size(); ++j)
                    check((found.count({i, j}) != 0) == (sketches[i].query(sketches[j]).value >= 0.05),
                          name + " exact index differs from exhaustive threshold query");
        }
        check(index.pairs(0.05, 1).size() == matches.size(), name + " index is thread-count dependent");
        check(!indexDict({projected, projected}, 2).pairs(0.5).empty(), name + " projected state cannot be indexed");
        auto wrong = config; wrong.seed += 1;
        const auto incompatible = build(wrong, {a});
        rejects<std::invalid_argument>([&] { left.merge(incompatible); }, name + " merge accepted a different seed");
        rejects<std::invalid_argument>([&] { indexDict({left, incompatible}); }, name + " index accepted a different seed");
        rejects<std::invalid_argument>([&] { index.search(0, std::numeric_limits<double>::quiet_NaN()); }, "NaN index threshold accepted");
        rejects<std::out_of_range>([&] { index.search(4); }, "invalid sample accepted");
        std::cout << name << " common operations OK\n";
    }
}

void testSetSketchRefinement() {
    // Coarse rank-zero partitions can contribute nonzero ranks after folding.
    for (const double base : {1.1, 2.0, 3.0}) {
        SetSketch fine(10, base, 5.0, 21, true, 42);
        SetSketch coarse(4, base, 5.0, 21, true, 42);
        for (uint64_t i = 1; i < 1024; ++i) {
            const uint64_t fingerprint = (i << 54) | 7;
            fine.addFingerprint(fingerprint);
            coarse.addFingerprint(fingerprint);
        }
        check(fine.project(4).getCore() == coarse.getCore(), "SetSketch lost censored rank-zero partitions");
        check(fine.project(7).project(4).getCore() == coarse.getCore(), "SetSketch chained projection differs");
        auto compact = SetSketch::fromRegisters(10, base, 5.0, 21, 42, fine.getCore());
        rejects<std::logic_error>([&] { compact.project(4); }, "register-only SetSketch incorrectly folded");
        rejects<std::invalid_argument>([&] { fine.project(11); }, "SetSketch increased precision");
    }
}

void testBinDashState() {
    const auto sequence = dna(2000, 81);
    for (const uint32_t bits : {1u, 16u, 64u}) {
        BinDash fine(6, 21, bits, 42), coarse(2, 21, bits, 42);
        fine.update(sequence.data(), sequence.size());
        coarse.update(sequence.data(), sequence.size());
        fine.finalize(); coarse.finalize();
        check(fine.project(128).storeHashes() == coarse.storeHashes(), "BinDash divisor projection differs");
        check(fine.merge(fine).storeHashes() == fine.storeHashes(), "BinDash merge after finalize failed");
        rejects<std::invalid_argument>([&] { fine.project(256); }, "nondivisor BinDash projection accepted");
        const auto packed = BinDash::fromPacked(6, 21, bits, 42, fine.getRawNonempty(), fine.packedSignatures());
        rejects<std::logic_error>([&] { packed.merge(fine); }, "BinDash merged lossy fingerprints as full hashes");
        check(indexDict(std::vector<BinDash>{fine, fine}, 2).pairs(0.9).size() == 1,
              "native BinDash index truncated a wide key");
    }
}

void testProbMinHashProjection() {
    using Semantics = Rank::WeightSemantics;
    ProbMinHash4 a(128, 21, 42, 8, Semantics::UserSupplied);
    ProbMinHash4 b(128, 21, 42, 8, Semantics::UserSupplied);
    for (uint64_t i = 1; i < 2000; ++i) {
        a.addHash(i, 1 + i % 5);
        b.addHash(i, 1 + i % 3);
    }
    a.finalize(); b.finalize();
    auto p = a.project(32), q = b.project(32);
    check(p.getSourceM() == 128 && p.getMaxL() == 8, "projection lost original race metadata");
    for (uint32_t i = 0; i < p.getM(); ++i)
        check(p.getRegisters()[i] == a.getRegisters()[i], "projection changed a race value");
    check(a.project(64).project(32).metadata().parameter_fingerprint == p.metadata().parameter_fingerprint,
          "chained ProbMinHash projection changed its sampling family");
    check(p.merge(q).storeHashes() == a.merge(b).project(32).storeHashes(), "weighted projection/merge mismatch");
    const auto index = indexDict(std::vector<ProbMinHash4>{p, q, p}, 3);
    for (const auto& match : index.pairs(0.0, 2)) {
        check(match.estimate.plan.metric == Query::Metric::WeightedJaccard, "weighted index used ordinary Jaccard");
        if (match.left == 0 && match.right == 1)
            check(match.estimate.value == p.jaccard_weighted(q), "weighted index ignored winner identity");
    }
    ProbMinHash4 independent(32, 21, 42, 8, Semantics::UserSupplied);
    check(!Query::query(p, independent).estimate_available, "projected and independent races compared silently");
    rejects<std::invalid_argument>([&] { p.merge(independent); }, "incompatible projected race merged");
    rejects<std::logic_error>([&] { p.addHash(1); }, "projected race was updated with smaller-size parameters");
    rejects<std::invalid_argument>([&] { p.project(64); }, "ProbMinHash expanded missing registers");
}

template<class T>
void testRegisterFinalize(T sketch, std::string sequence) {
    check(!sketch.isSealed(), "register sketch starts sealed");
    sketch.update(sequence.data());
    const auto core = sketch.getCore();
    const auto witnesses = sketch.getWitnesses();
    const auto projected = sketch.project(6).getCore();
    const double cardinality = sketch.cardinality();
    sketch.finalize();
    sketch.finalize();
    check(sketch.isSealed() && sketch.getCore() == core &&
          sketch.getWitnesses() == witnesses && sketch.cardinality() == cardinality,
          "finalize changed register data or failed to seal");
    rejects<std::logic_error>([&] { sketch.update(sequence.data()); },
                             "register update accepted after finalize");
    rejects<std::logic_error>([&] { sketch.addFingerprint(123); },
                             "fingerprint update accepted after finalize");
    T copy = sketch;
    rejects<std::logic_error>([&] { copy.addFingerprint(123); },
                             "copy lost the sealed state");
    check(sketch.project(6).getCore() == projected, "finalize disabled exact projection");
    check(sketch.merge(sketch).getCore() == core, "finalize disabled register merge");
    auto index = indexDict(std::vector<T>{sketch, sketch});
    check(index.size() == 2 && index.candidates(0).size() == 1,
          "finalize discarded index witnesses");
}

void testFinalizeLifecycle() {
    for (auto sequence : {std::string(), dna(2000, 101)}) {
        testRegisterFinalize(HyperLogLog(8, 21, 42, true), sequence);
        testRegisterFinalize(SetSketch(8, 2.0, 5.0, 21, true, 42), sequence);
    }
    SetSketch set(8);
    set.finalize();
    rejects<std::logic_error>([&] { set.update(nullptr, 0); },
                             "length-aware update bypassed SetSketch sealing");

    const auto config = configuration(Algorithm::FastKMV);
    MultiSketchBuilder modern({config}), legacy({config}), no_label({config});
    IO::FastxRecord record;
    record.name = "record";
    record.sequence = dna(2000, 102);
    modern.update(record);
    legacy.update(record);
    no_label.update(record);
    const auto a = modern.finalize("sample", "memory", record.name);
    const auto b = legacy.finish("sample", "memory", record.name);
    check(a.front().sketch.storeHashes() == b.front().sketch.storeHashes() &&
          a.front().label == b.front().label && a.front().source == b.front().source &&
          a.front().record_name == b.front().record_name,
          "finish compatibility alias changed the result");
    check(modern.finished() && legacy.finished() && a.front().stats.sealed,
          "builder finalize did not complete the lifecycle");
    check(no_label.finalize().front().label.find("sketch:") == 0,
          "default finalize label is missing");
    rejects<std::logic_error>([&] { modern.update(record); }, "builder updated after finalize");
    rejects<std::logic_error>([&] { modern.finish("again"); }, "finish reopened a finalized builder");
    rejects<std::logic_error>([&] { legacy.finalize(); }, "finalize reopened a finished builder");
}

void testEmptyAndInvalid() {
    const auto config = configuration(Algorithm::FastKMV);
    const auto empty = build(config, {});
    const auto empty_index = indexDict({empty, empty}, 2);
    const auto empty_pairs = empty_index.pairs(0.0);
    check(empty_pairs.size() == 1 && empty_pairs.front().estimate.value == empty.query(empty).value,
          "empty/empty index result differs from native empty-set convention");
    check(empty_index.pairs(0.9).size() == (empty.query(empty).value >= 0.9 ? 1u : 0u),
          "empty/empty index threshold differs from native query");
    check(indexDict(std::vector<BuiltSketch>{}).pairs().empty(), "empty index produced matches");
    rejects<std::invalid_argument>([&] { indexDict({empty}, 0); }, "zero index thread count accepted");
    const auto hll = build(configuration(Algorithm::HyperLogLog), {});
    rejects<std::invalid_argument>([&] { empty.merge(hll); }, "cross-backend merge accepted");
    rejects<std::invalid_argument>([&] { indexDict({empty, hll}); }, "mixed index accepted");
}
} // namespace

int main() {
    try {
        testSevenBackends();
        testSetSketchRefinement();
        testBinDashState();
        testProbMinHashProjection();
        testFinalizeLifecycle();
        testEmptyAndInvalid();
        std::cout << "sketch operations tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
