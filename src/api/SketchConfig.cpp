#include "api/SketchConfig.h"

#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace Sketch {
namespace API {
namespace {

uint64_t addFlag(uint64_t flags, Capability capability) noexcept {
    return flags | static_cast<uint64_t>(capability);
}

void hashByte(uint64_t& hash, uint8_t value) noexcept {
    hash ^= value;
    hash *= UINT64_C(1099511628211);
}

void hashInteger(uint64_t& hash, uint64_t value, size_t width) noexcept {
    for (size_t i = 0; i < width; ++i)
        hashByte(hash, static_cast<uint8_t>(value >> (8 * i)));
}

} // namespace

void SketchConfig::validate() const {
    if (kmer_size == 0 || kmer_size > 32)
        throw std::invalid_argument("SketchConfig k-mer size must be in [1, 32]");
    if (resolution == 0 && algorithm != Algorithm::KSSD)
        throw std::invalid_argument("SketchConfig resolution must be positive");
    if (sampling == SamplingMode::BottomK && resolution < 2)
        throw std::invalid_argument("bottom-k resolution must be at least 2");
    if (sampling == SamplingMode::Scaled && scaled == 0)
        throw std::invalid_argument("scaled sampling requires scaled > 0");
    if (sampling != SamplingMode::Scaled && scaled != 0)
        throw std::invalid_argument(
            "scaled is only valid when sampling mode is Scaled");
    if (min_abundance == 0)
        throw std::invalid_argument("minimum abundance must be positive");
    if (max_abundance != 0 && max_abundance < min_abundance)
        throw std::invalid_argument(
            "maximum abundance must be zero or at least minimum abundance");
    if (minimum_base_quality > 93)
        throw std::invalid_argument("minimum base quality must be <= 93");
    if (bindash_bits == 0 || bindash_bits > 64)
        throw std::invalid_argument("BinDash bit width must be in [1, 64]");
    if (!(setsketch_base > 1.0) || !std::isfinite(setsketch_base) ||
        !(setsketch_a > 0.0) || !std::isfinite(setsketch_a))
        throw std::invalid_argument(
            "SetSketch base/a must be finite with base > 1 and a > 0");
    if (kssd_half_subk < 3 || kssd_half_subk > 7)
        throw std::invalid_argument("KSSD half-sub-k must be in [3, 7]");
    if (order_l > 5 || order_m == 0 ||
        order_m > static_cast<uint32_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument(
            "OrderMinHash requires order_l in [0, 5] and order_m in [1, INT_MAX]");
    if (!(target_relative_error >= 0.0) ||
        !std::isfinite(target_relative_error) || target_relative_error >= 1.0)
        throw std::invalid_argument(
            "target relative error must be finite and in [0, 1)");
    if (track_abundance &&
        weight_semantics == Rank::WeightSemantics::UnweightedSet)
        throw std::invalid_argument(
            "abundance tracking requires frequency or user-supplied weight semantics");
    if (!track_abundance &&
        weight_semantics != Rank::WeightSemantics::UnweightedSet &&
        algorithm != Algorithm::ProbMinHash)
        throw std::invalid_argument(
            "weighted semantics require abundance tracking for this algorithm");
    if (track_abundance && algorithm != Algorithm::ProbMinHash)
        throw std::invalid_argument(
            "the unified FASTX builder currently supports abundance only for ProbMinHash");
    if (track_abundance &&
        weight_semantics != Rank::WeightSemantics::Frequency)
        throw std::invalid_argument(
            "FASTX abundance builds require frequency weight semantics");
    if (!track_abundance && (min_abundance != 1 || max_abundance != 0))
        throw std::invalid_argument(
            "abundance filters require track_abundance=true");

    switch (algorithm) {
        case Algorithm::LegacyMinHash:
        case Algorithm::FastKMV:
            if (sampling != SamplingMode::BottomK)
                throw std::invalid_argument(
                    "legacy/FastKMV require bottom-k sampling");
            break;
        case Algorithm::FracMinHash:
            if (sampling != SamplingMode::Scaled || scaled == 0 ||
                scaled > std::numeric_limits<uint32_t>::max())
                throw std::invalid_argument(
                    "FracMinHash requires a 32-bit positive scaled value");
            break;
        case Algorithm::HyperLogLog:
            if (sampling != SamplingMode::RegisterPrecision ||
                resolution < 4 || resolution > 20)
                throw std::invalid_argument(
                    "HLL requires register precision in [4, 20]");
            break;
        case Algorithm::SetSketch:
            if (sampling != SamplingMode::RegisterPrecision ||
                resolution < 4 || resolution > 16)
                throw std::invalid_argument(
                    "SetSketch requires register precision in [4, 16]");
            break;
        case Algorithm::KSSD:
            if (sampling != SamplingMode::KssdReduction ||
                kmer_size % 2 != 0 || kmer_size / 2 < kssd_half_subk ||
                kmer_size / 2 > 16 || resolution > kssd_half_subk - 3)
                throw std::invalid_argument(
                    "KSSD requires even k, half_subk<=k/2<=16, and a valid drlevel");
            break;
        case Algorithm::ProbMinHash:
        case Algorithm::BinDash:
        case Algorithm::OrderMinHash:
            if (sampling != SamplingMode::AlgorithmNative)
                throw std::invalid_argument(
                    "ProbMinHash/BinDash/OrderMinHash require algorithm-native sampling");
            if (algorithm == Algorithm::ProbMinHash && resolution < 2)
                throw std::invalid_argument(
                    "ProbMinHash requires at least two registers");
            if (algorithm == Algorithm::BinDash &&
                (resolution < 64 || resolution % 64 != 0))
                throw std::invalid_argument(
                    "BinDash resolution must be a positive multiple of 64 bins");
            break;
    }
}

uint64_t SketchConfig::fingerprint() const noexcept {
    uint64_t hash = UINT64_C(1469598103934665603);
    hashInteger(hash, static_cast<uint8_t>(algorithm), 1);
    hashInteger(hash, static_cast<uint8_t>(molecule), 1);
    hashInteger(hash, static_cast<uint8_t>(ambiguous_policy), 1);
    hashInteger(hash, static_cast<uint8_t>(aggregation), 1);
    hashInteger(hash, static_cast<uint8_t>(sampling), 1);
    hashInteger(hash, static_cast<uint8_t>(weight_semantics), 1);
    hashInteger(hash, kmer_size, 2);
    hashInteger(hash, seed, 8);
    hashInteger(hash, resolution, 4);
    hashInteger(hash, scaled, 8);
    hashInteger(hash, canonical ? 1 : 0, 1);
    hashInteger(hash, track_abundance ? 1 : 0, 1);
    hashInteger(hash, homopolymer_compressed ? 1 : 0, 1);
    hashInteger(hash, min_abundance, 4);
    hashInteger(hash, max_abundance, 4);
    hashInteger(hash, static_cast<uint16_t>(minimum_base_quality), 2);
    hashInteger(hash, bindash_bits, 1);
    hashInteger(hash, probminhash_max_l, 4);
    uint64_t set_base_bits = 0;
    uint64_t set_a_bits = 0;
    std::memcpy(&set_base_bits, &setsketch_base, sizeof(set_base_bits));
    std::memcpy(&set_a_bits, &setsketch_a, sizeof(set_a_bits));
    hashInteger(hash, set_base_bits, 8);
    hashInteger(hash, set_a_bits, 8);
    hashInteger(hash, kssd_half_subk, 1);
    hashInteger(hash, order_l, 2);
    hashInteger(hash, order_m, 4);
    hashInteger(hash, memory_budget_bytes, 8);
    uint64_t error_bits = 0;
    static_assert(sizeof(error_bits) == sizeof(target_relative_error),
                  "unexpected double width");
    std::memcpy(&error_bits, &target_relative_error, sizeof(error_bits));
    hashInteger(hash, error_bits, 8);
    return hash;
}

Capabilities capabilitiesFor(Algorithm algorithm) noexcept {
    Capabilities result;
    result.algorithm = algorithm;
    uint64_t flags = 0;
    flags = addFlag(flags, Capability::SequenceUpdate);
    flags = addFlag(flags, Capability::RecordUpdate);
    flags = addFlag(flags, Capability::Seal);

    switch (algorithm) {
        case Algorithm::LegacyMinHash:
        case Algorithm::FastKMV:
        case Algorithm::FracMinHash:
            flags = addFlag(flags, Capability::HashUpdate);
            flags = addFlag(flags, Capability::Cardinality);
            flags = addFlag(flags, Capability::Jaccard);
            flags = addFlag(flags, Capability::JaccardDistance);
            flags = addFlag(flags, Capability::MashDistance);
            flags = addFlag(flags, Capability::Containment);
            flags = addFlag(flags, Capability::Intersection);
            flags = addFlag(flags, Capability::Union);
            flags = addFlag(flags, Capability::ANI);
            flags = addFlag(flags, Capability::Merge);
            flags = addFlag(flags, Capability::Project);
            flags = addFlag(flags, Capability::SearchKeys);
            break;
        case Algorithm::HyperLogLog:
            flags = addFlag(flags, Capability::HashUpdate);
            flags = addFlag(flags, Capability::Cardinality);
            flags = addFlag(flags, Capability::Jaccard);
            flags = addFlag(flags, Capability::JaccardDistance);
            flags = addFlag(flags, Capability::MashDistance);
            flags = addFlag(flags, Capability::Containment);
            flags = addFlag(flags, Capability::Intersection);
            flags = addFlag(flags, Capability::Union);
            flags = addFlag(flags, Capability::ANI);
            flags = addFlag(flags, Capability::Merge);
            flags = addFlag(flags, Capability::Project);
            flags = addFlag(flags, Capability::SearchKeys);
            break;
        case Algorithm::SetSketch:
            flags = addFlag(flags, Capability::HashUpdate);
            flags = addFlag(flags, Capability::Cardinality);
            flags = addFlag(flags, Capability::Jaccard);
            flags = addFlag(flags, Capability::JaccardDistance);
            flags = addFlag(flags, Capability::MashDistance);
            flags = addFlag(flags, Capability::Containment);
            flags = addFlag(flags, Capability::Intersection);
            flags = addFlag(flags, Capability::Union);
            flags = addFlag(flags, Capability::ANI);
            flags = addFlag(flags, Capability::Merge);
            flags = addFlag(flags, Capability::SearchKeys);
            break;
        case Algorithm::KSSD:
            flags = addFlag(flags, Capability::Jaccard);
            flags = addFlag(flags, Capability::JaccardDistance);
            flags = addFlag(flags, Capability::MashDistance);
            flags = addFlag(flags, Capability::ANI);
            flags = addFlag(flags, Capability::Project);
            flags = addFlag(flags, Capability::SearchKeys);
            break;
        case Algorithm::ProbMinHash:
            flags = addFlag(flags, Capability::HashUpdate);
            flags = addFlag(flags, Capability::WeightedUpdate);
            flags = addFlag(flags, Capability::Jaccard);
            flags = addFlag(flags, Capability::JaccardDistance);
            flags = addFlag(flags, Capability::WeightedJaccard);
            flags = addFlag(flags, Capability::WeightedContainment);
            flags = addFlag(flags, Capability::MashDistance);
            flags = addFlag(flags, Capability::ANI);
            flags = addFlag(flags, Capability::Merge);
            flags = addFlag(flags, Capability::SearchKeys);
            break;
        case Algorithm::BinDash:
            flags = addFlag(flags, Capability::Cardinality);
            flags = addFlag(flags, Capability::Jaccard);
            flags = addFlag(flags, Capability::JaccardDistance);
            flags = addFlag(flags, Capability::MashDistance);
            flags = addFlag(flags, Capability::Containment);
            flags = addFlag(flags, Capability::Intersection);
            flags = addFlag(flags, Capability::Union);
            flags = addFlag(flags, Capability::ANI);
            flags = addFlag(flags, Capability::SearchKeys);
            break;
        case Algorithm::OrderMinHash:
            flags = addFlag(flags, Capability::OrderSimilarity);
            flags = addFlag(flags, Capability::OrderDistance);
            flags = addFlag(flags, Capability::SearchKeys);
            break;
    }
    result.flags = flags;
    return result;
}

const char* algorithmName(Algorithm algorithm) noexcept {
    switch (algorithm) {
        case Algorithm::LegacyMinHash: return "minhash";
        case Algorithm::FastKMV: return "fastkmv";
        case Algorithm::FracMinHash: return "fracminhash";
        case Algorithm::HyperLogLog: return "hll";
        case Algorithm::SetSketch: return "setsketch";
        case Algorithm::KSSD: return "kssd";
        case Algorithm::ProbMinHash: return "probminhash";
        case Algorithm::BinDash: return "bindash";
        case Algorithm::OrderMinHash: return "orderminhash";
    }
    return "unknown";
}

const char* moleculeTypeName(MoleculeType molecule) noexcept {
    switch (molecule) {
        case MoleculeType::DNA: return "dna";
        case MoleculeType::RNA: return "rna";
        case MoleculeType::Protein: return "protein";
        case MoleculeType::Dayhoff: return "dayhoff";
        case MoleculeType::HydrophobicPolar: return "hp";
    }
    return "unknown";
}

const char* samplingModeName(SamplingMode sampling) noexcept {
    switch (sampling) {
        case SamplingMode::BottomK: return "bottom_k";
        case SamplingMode::Scaled: return "scaled";
        case SamplingMode::RegisterPrecision: return "register_precision";
        case SamplingMode::KssdReduction: return "kssd_reduction";
        case SamplingMode::AlgorithmNative: return "algorithm_native";
    }
    return "unknown";
}

} // namespace API
} // namespace Sketch
