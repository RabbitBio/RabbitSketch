#ifndef RABBITSKETCH_API_SKETCH_CONFIG_H
#define RABBITSKETCH_API_SKETCH_CONFIG_H

#include "rank/RankMetadata.h"

#include <cstdint>
#include <string>

namespace Sketch {
namespace API {

/** Public algorithm names.  These are intentionally broader than
 * Rank::Backend: legacy, compact, and order-sensitive implementations do not
 * all share the RabbitSketch2 rank-space contract. */
enum class Algorithm : uint8_t {
    LegacyMinHash = 1,
    FastKMV = 2,
    FracMinHash = 3,
    HyperLogLog = 4,
    SetSketch = 5,
    KSSD = 6,
    ProbMinHash = 7,
    BinDash = 8,
    OrderMinHash = 9
};

enum class MoleculeType : uint8_t {
    DNA = 1,
    RNA = 2,
    Protein = 3,
    Dayhoff = 4,
    HydrophobicPolar = 5
};

enum class AmbiguousPolicy : uint8_t {
    SkipKmer = 1,
    RejectRecord = 2,
    MapToUnknown = 3
};

enum class AggregationMode : uint8_t {
    OneSketchPerFile = 1,
    OneSketchPerRecord = 2,
    OneSketchPerCollection = 3
};

enum class SamplingMode : uint8_t {
    BottomK = 1,
    Scaled = 2,
    RegisterPrecision = 3,
    KssdReduction = 4,
    AlgorithmNative = 5
};

/** Stable, backend-independent build contract.  Backend-specific constructors
 * remain available, but high-level builders consume this type so C++ and
 * Python cannot silently choose different sequence semantics. */
struct SketchConfig {
    Algorithm algorithm = Algorithm::FastKMV;
    MoleculeType molecule = MoleculeType::DNA;
    AmbiguousPolicy ambiguous_policy = AmbiguousPolicy::SkipKmer;
    AggregationMode aggregation = AggregationMode::OneSketchPerFile;
    SamplingMode sampling = SamplingMode::BottomK;
    Rank::WeightSemantics weight_semantics =
        Rank::WeightSemantics::UnweightedSet;

    uint16_t kmer_size = 21;
    uint64_t seed = 42;
    uint32_t resolution = 1024;
    uint64_t scaled = 0;
    bool canonical = true;
    bool track_abundance = false;
    bool homopolymer_compressed = false;
    uint32_t min_abundance = 1;
    uint32_t max_abundance = 0;  // zero means unbounded
    int16_t minimum_base_quality = -1;  // negative means disabled

    // Backend-specific parameters that affect the statistical
    // contract.  They live here (instead of hidden constructor defaults) so a
    // build can be reproduced from its public configuration alone.
    uint8_t bindash_bits = 16;
    uint32_t probminhash_max_l = 0;
    double setsketch_base = 2.0;
    double setsketch_a = 5.0;
    uint8_t kssd_half_subk = 6;
    uint16_t order_l = 2;
    uint32_t order_m = 500;

    // Optional task-level budgets.  Zero means that the backend-native
    // resolution fields above are authoritative.
    uint64_t memory_budget_bytes = 0;
    double target_relative_error = 0.0;

    void validate() const;
    uint64_t fingerprint() const noexcept;
};

enum class Capability : uint64_t {
    SequenceUpdate = UINT64_C(1) << 0,
    RecordUpdate = UINT64_C(1) << 1,
    HashUpdate = UINT64_C(1) << 2,
    WeightedUpdate = UINT64_C(1) << 3,
    Seal = UINT64_C(1) << 4,
    Cardinality = UINT64_C(1) << 5,
    Jaccard = UINT64_C(1) << 6,
    JaccardDistance = UINT64_C(1) << 7,
    MashDistance = UINT64_C(1) << 8,
    Containment = UINT64_C(1) << 9,
    Intersection = UINT64_C(1) << 10,
    Union = UINT64_C(1) << 11,
    ANI = UINT64_C(1) << 12,
    WeightedJaccard = UINT64_C(1) << 13,
    WeightedContainment = UINT64_C(1) << 14,
    OrderSimilarity = UINT64_C(1) << 15,
    Merge = UINT64_C(1) << 16,
    Project = UINT64_C(1) << 17,
    SearchKeys = UINT64_C(1) << 18,
    OrderDistance = UINT64_C(1) << 19
};

struct Capabilities {
    Algorithm algorithm = Algorithm::FastKMV;
    uint64_t flags = 0;

    bool supports(Capability capability) const noexcept {
        return (flags & static_cast<uint64_t>(capability)) != 0;
    }
};

struct BuildStats {
    uint64_t records = 0;
    uint64_t input_bases = 0;
    uint64_t candidate_kmers = 0;
    uint64_t accepted_kmers = 0;
    uint64_t skipped_ambiguous_kmers = 0;
    uint64_t skipped_low_quality_kmers = 0;
    uint64_t distinct_kmers = 0;
    uint64_t accepted_distinct_kmers = 0;
    uint64_t skipped_abundance_kmers = 0;
    bool sealed = false;
};

Capabilities capabilitiesFor(Algorithm algorithm) noexcept;
const char* algorithmName(Algorithm algorithm) noexcept;
const char* moleculeTypeName(MoleculeType molecule) noexcept;
const char* samplingModeName(SamplingMode sampling) noexcept;

} // namespace API
} // namespace Sketch

#endif // RABBITSKETCH_API_SKETCH_CONFIG_H
