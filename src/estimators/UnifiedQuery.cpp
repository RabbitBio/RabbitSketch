#include "estimators/UnifiedQuery.h"

#include "BinDash.h"
#include "SetSketch.h"
#include "Sketch.h"
#include "fastkmv.h"
#include "fracminhash.h"
#include "probmh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Sketch {
namespace Query {
namespace {

bool isKnownMetric(Metric metric) noexcept {
    const uint8_t value = static_cast<uint8_t>(metric);
    return value >= static_cast<uint8_t>(Metric::Jaccard) &&
           value <= static_cast<uint8_t>(Metric::OrderDistance);
}

bool kssdMetricSupported(Metric metric) noexcept {
    return metric == Metric::Jaccard || metric == Metric::ANI ||
           metric == Metric::JaccardDistance ||
           metric == Metric::MashDistance;
}

bool setMetricSupported(Metric metric) noexcept {
    return metric != Metric::WeightedJaccard &&
           metric != Metric::WeightedContainment &&
           metric != Metric::OrderSimilarity &&
           metric != Metric::OrderDistance;
}

bool probMinHashMetricSupported(
    Metric metric, Rank::WeightSemantics semantics) noexcept {
    if (metric == Metric::WeightedJaccard ||
        metric == Metric::WeightedContainment)
        return true;
    if (semantics != Rank::WeightSemantics::UnweightedSet)
        return false;
    return metric == Metric::Jaccard ||
           metric == Metric::JaccardDistance ||
           metric == Metric::MashDistance ||
           metric == Metric::ANI;
}

bool orderMetricSupported(Metric metric) noexcept {
    return metric == Metric::OrderSimilarity ||
           metric == Metric::OrderDistance;
}

void validateRequest(const Request& request) {
    if (!isKnownMetric(request.metric))
        throw std::invalid_argument("RabbitSketch2 query metric is invalid");
    if (!(request.confidence_level > 0.0 && request.confidence_level < 1.0) ||
        !std::isfinite(request.confidence_level))
        throw std::invalid_argument(
            "RabbitSketch2 confidence level must be finite and in (0, 1)");
}

Plan basePlan(const SketchDescriptor& left,
              const SketchDescriptor& right,
              const Request& request) {
    Plan plan;
    plan.metric = request.metric;
    plan.left_backend = left.metadata.backend;
    plan.right_backend = right.metadata.backend;
    return plan;
}

void setIncompatible(Plan& plan, const std::string& mismatch) {
    plan.support = Support::Incompatible;
    plan.method = Method::None;
    plan.reason = "incompatible Rank-Space metadata: " + mismatch;
}

std::string backendContractMismatch(const SketchDescriptor& descriptor) {
    const Rank::RankMetadata& metadata = descriptor.metadata;
    switch (metadata.backend) {
        case Rank::Backend::FastKMV:
            if (metadata.resolution_kind != Rank::ResolutionKind::BottomK)
                return "FastKMV resolution kind";
            if (metadata.fingerprint_bits != 53)
                return "FastKMV fingerprint width";
            break;
        case Rank::Backend::KSSD:
            if (metadata.resolution_kind != Rank::ResolutionKind::KssdDrLevel)
                return "KSSD resolution kind";
            if (metadata.fingerprint_bits != 64)
                return "KSSD fingerprint width";
            break;
        case Rank::Backend::SetSketch:
        case Rank::Backend::HLL:
            if (metadata.resolution_kind !=
                Rank::ResolutionKind::RegisterPrecisionBits)
                return "register-sketch resolution kind";
            if (metadata.fingerprint_bits != 64)
                return "register-sketch fingerprint width";
            break;
        case Rank::Backend::ProbMinHash:
            if (metadata.resolution_kind != Rank::ResolutionKind::RegisterCount)
                return "ProbMinHash resolution kind";
            if (metadata.fingerprint_bits != 64)
                return "ProbMinHash fingerprint width";
            break;
        case Rank::Backend::BinDash:
            if (metadata.resolution_kind != Rank::ResolutionKind::RegisterCount)
                return "BinDash resolution kind";
            break;
        case Rank::Backend::LegacyMinHash:
            if (metadata.resolution_kind != Rank::ResolutionKind::BottomK)
                return "legacy MinHash resolution kind";
            if (metadata.fingerprint_bits != 64)
                return "legacy MinHash fingerprint width";
            break;
        case Rank::Backend::FracMinHash:
            if (metadata.resolution_kind != Rank::ResolutionKind::Scaled)
                return "FracMinHash resolution kind";
            if (metadata.fingerprint_bits != 64)
                return "FracMinHash fingerprint width";
            break;
        case Rank::Backend::OrderMinHash:
            if (metadata.resolution_kind != Rank::ResolutionKind::OrderSampleCount)
                return "OrderMinHash resolution kind";
            if (metadata.fingerprint_bits != 64)
                return "OrderMinHash fingerprint width";
            break;
    }
    return std::string();
}

double aniFromJaccard(double jaccard, uint16_t kmer_size) {
    if (!(jaccard > 0.0)) return 0.0;
    if (jaccard >= 1.0) return 1.0;
    return std::pow(2.0 * jaccard / (1.0 + jaccard),
                    1.0 / static_cast<double>(kmer_size));
}

double aniFromContainment(double containment, uint16_t kmer_size) {
    if (!(containment > 0.0)) return 0.0;
    if (containment >= 1.0) return 1.0;
    return std::pow(containment, 1.0 / static_cast<double>(kmer_size));
}

double mashDistanceFromJaccard(double jaccard, uint16_t kmer_size) {
    if (!(jaccard > 0.0)) return std::numeric_limits<double>::infinity();
    if (jaccard >= 1.0) return 0.0;
    return -std::log(2.0 * jaccard / (1.0 + jaccard)) /
           static_cast<double>(kmer_size);
}

double bottomKJaccard(const std::vector<uint64_t>& left,
                      const std::vector<uint64_t>& right,
                      uint32_t maximum_union_samples) {
    size_t left_index = 0;
    size_t right_index = 0;
    uint64_t common = 0;
    uint64_t sampled_union = 0;
    while (sampled_union < maximum_union_samples &&
           (left_index < left.size() || right_index < right.size())) {
        if (right_index >= right.size() ||
            (left_index < left.size() &&
             left[left_index] < right[right_index])) {
            ++left_index;
        } else if (left_index >= left.size() ||
                   right[right_index] < left[left_index]) {
            ++right_index;
        } else {
            ++left_index;
            ++right_index;
            ++common;
        }
        ++sampled_union;
    }
    return sampled_union == 0
        ? 1.0 : static_cast<double>(common) / sampled_union;
}

void populateDerivedMetrics(Result& result, uint16_t kmer_size) {
    if (std::isfinite(result.jaccard)) {
        result.jaccard = std::max(0.0, std::min(1.0, result.jaccard));
        result.jaccard_distance = 1.0 - result.jaccard;
        result.mash_distance =
            mashDistanceFromJaccard(result.jaccard, kmer_size);
        result.ani = aniFromJaccard(result.jaccard, kmer_size);
    }
    if (std::isfinite(result.left_cardinality) &&
        std::isfinite(result.right_cardinality) &&
        std::isfinite(result.intersection)) {
        result.union_size = std::max(
            0.0, result.left_cardinality + result.right_cardinality -
                     result.intersection);
    }
    if (std::isfinite(result.left_containment) &&
        std::isfinite(result.right_containment)) {
        result.left_containment = std::max(
            0.0, std::min(1.0, result.left_containment));
        result.right_containment = std::max(
            0.0, std::min(1.0, result.right_containment));
        result.max_containment = std::max(
            result.left_containment, result.right_containment);
        result.average_containment = 0.5 *
            (result.left_containment + result.right_containment);
        result.left_containment_ani =
            aniFromContainment(result.left_containment, kmer_size);
        result.right_containment_ani =
            aniFromContainment(result.right_containment, kmer_size);
        result.max_containment_ani = std::max(
            result.left_containment_ani, result.right_containment_ani);
    }
}

void selectMetricValue(Result& result, Metric metric) {
    switch (metric) {
        case Metric::Jaccard: result.value = result.jaccard; break;
        case Metric::JaccardDistance:
            result.value = result.jaccard_distance;
            break;
        case Metric::MashDistance: result.value = result.mash_distance; break;
        case Metric::Intersection: result.value = result.intersection; break;
        case Metric::Union: result.value = result.union_size; break;
        case Metric::LeftContainment:
            result.value = result.left_containment;
            break;
        case Metric::RightContainment:
            result.value = result.right_containment;
            break;
        case Metric::MaxContainment:
            result.value = result.max_containment;
            break;
        case Metric::AverageContainment:
            result.value = result.average_containment;
            break;
        case Metric::ANI: result.value = result.ani; break;
        case Metric::LeftContainmentANI:
            result.value = result.left_containment_ani;
            break;
        case Metric::RightContainmentANI:
            result.value = result.right_containment_ani;
            break;
        case Metric::MaxContainmentANI:
            result.value = result.max_containment_ani;
            break;
        case Metric::WeightedJaccard:
            result.value = result.weighted_jaccard;
            break;
        case Metric::WeightedContainment:
            result.value = result.weighted_left_containment;
            break;
        case Metric::OrderSimilarity:
            result.value = result.order_similarity;
            break;
        case Metric::OrderDistance:
            result.value = result.order_distance;
            break;
    }
}

void populateCardinalityMetrics(Result& result,
                                double jaccard,
                                double left_cardinality,
                                double right_cardinality,
                                uint16_t kmer_size) {
    result.estimate_available = true;
    result.jaccard = jaccard;
    result.left_cardinality = std::max(0.0, left_cardinality);
    result.right_cardinality = std::max(0.0, right_cardinality);
    if (!std::isfinite(result.left_cardinality) ||
        !std::isfinite(result.right_cardinality)) {
        result.warnings.push_back(
            "cardinality-derived metrics are unavailable because a cardinality estimate saturated");
        populateDerivedMetrics(result, kmer_size);
        selectMetricValue(result, result.plan.metric);
        return;
    }
    const double cardinality_sum =
        result.left_cardinality + result.right_cardinality;
    result.intersection = result.jaccard * cardinality_sum /
        (1.0 + result.jaccard);
    result.intersection = std::max(
        0.0, std::min(result.intersection,
                      std::min(result.left_cardinality,
                               result.right_cardinality)));
    result.left_containment = result.left_cardinality > 0.0
        ? result.intersection / result.left_cardinality : 0.0;
    result.right_containment = result.right_cardinality > 0.0
        ? result.intersection / result.right_cardinality : 0.0;
    populateDerivedMetrics(result, kmer_size);
    selectMetricValue(result, result.plan.metric);
}

} // namespace

SketchDescriptor describe(const FastKMV& sketch) {
    SketchDescriptor descriptor;
    descriptor.metadata = sketch.metadata();
    return descriptor;
}

SketchDescriptor describe(const FracMinHash& sketch) {
    SketchDescriptor descriptor;
    descriptor.metadata = sketch.metadata();
    return descriptor;
}

SketchDescriptor describe(const BinDash& sketch) {
    SketchDescriptor descriptor;
    descriptor.metadata = sketch.metadata();
    return descriptor;
}

SketchDescriptor describe(const HyperLogLog& sketch) {
    SketchDescriptor descriptor;
    descriptor.metadata = sketch.metadata();
    return descriptor;
}

SketchDescriptor describe(const SetSketch& sketch) {
    SketchDescriptor descriptor;
    descriptor.metadata = sketch.metadata();
    return descriptor;
}

SketchDescriptor describe(const Kssd& sketch) {
    SketchDescriptor descriptor;
    descriptor.metadata = sketch.metadata();
    return descriptor;
}

SketchDescriptor describe(const MinHash& sketch) {
    SketchDescriptor descriptor;
    descriptor.metadata = sketch.metadata();
    return descriptor;
}

SketchDescriptor describe(const OrderMinHash& sketch) {
    SketchDescriptor descriptor;
    descriptor.metadata = sketch.metadata();
    return descriptor;
}

SketchDescriptor describe(const ProbMinHash4& sketch) {
    SketchDescriptor descriptor;
    descriptor.metadata = sketch.metadata();
    return descriptor;
}

Plan planQuery(const SketchDescriptor& left,
               const SketchDescriptor& right,
               const Request& request) {
    validateRequest(request);
    Plan plan = basePlan(left, right, request);

    const std::string left_contract = backendContractMismatch(left);
    if (!left_contract.empty()) {
        setIncompatible(plan, "left " + left_contract);
        return plan;
    }
    const std::string right_contract = backendContractMismatch(right);
    if (!right_contract.empty()) {
        setIncompatible(plan, "right " + right_contract);
        return plan;
    }

    const Rank::Backend left_backend = left.metadata.backend;
    const Rank::Backend right_backend = right.metadata.backend;

    const bool involves_native_set_sketch =
        left_backend == Rank::Backend::FastKMV ||
        left_backend == Rank::Backend::BinDash ||
        left_backend == Rank::Backend::LegacyMinHash ||
        left_backend == Rank::Backend::FracMinHash ||
        left_backend == Rank::Backend::HLL ||
        left_backend == Rank::Backend::SetSketch ||
        right_backend == Rank::Backend::FastKMV ||
        right_backend == Rank::Backend::BinDash ||
        right_backend == Rank::Backend::LegacyMinHash ||
        right_backend == Rank::Backend::FracMinHash ||
        right_backend == Rank::Backend::HLL ||
        right_backend == Rank::Backend::SetSketch;
    if (involves_native_set_sketch && !setMetricSupported(request.metric)) {
        plan.reason =
            "requested weighted/order metric is not defined for this set-sketch pair";
        return plan;
    }

    if (left_backend == right_backend) {
        const std::string mismatch = left.metadata.sameFamilyMismatch(
            right.metadata, /*allow_resolution_change=*/true);
        if (!mismatch.empty()) {
            setIncompatible(plan, mismatch);
            return plan;
        }

        switch (left_backend) {
            case Rank::Backend::FastKMV:
                plan.support = Support::Supported;
                plan.method = Method::FastKmvNative;
                plan.common_resolution = std::min(
                    left.metadata.resolution, right.metadata.resolution);
                plan.exact_projection =
                    left.metadata.resolution != right.metadata.resolution;
                plan.reason = plan.exact_projection
                    ? "exact bottom-K prefix projection to common K"
                    : "direct coordinated FastKMV comparison";
                return plan;

            case Rank::Backend::HLL:
                plan.support = Support::Supported;
                plan.method = Method::HllNative;
                plan.common_resolution = std::min(
                    left.metadata.resolution, right.metadata.resolution);
                plan.exact_projection =
                    left.metadata.resolution != right.metadata.resolution;
                plan.reason = plan.exact_projection
                    ? "exact HLL register folding to common precision"
                    : "direct coordinated HLL comparison";
                return plan;

            case Rank::Backend::SetSketch:
                plan.method = Method::SetSketchNative;
                if (left.metadata.resolution != right.metadata.resolution) {
                    plan.reason =
                        "SetSketch lower-precision folding is not distribution-proven";
                    return plan;
                }
                plan.support = Support::Supported;
                plan.common_resolution = left.metadata.resolution;
                plan.reason = "direct same-precision SetSketch comparison";
                return plan;

            case Rank::Backend::KSSD:
                plan.method = Method::KssdNative;
                if (!kssdMetricSupported(request.metric)) {
                    plan.reason =
                        "native KSSD payload exposes Jaccard/ANI but not cardinality-derived estimands";
                    return plan;
                }
                plan.support = Support::Supported;
                plan.common_resolution = std::max(
                    left.metadata.resolution, right.metadata.resolution);
                plan.exact_projection =
                    left.metadata.resolution != right.metadata.resolution;
                plan.reason = plan.exact_projection
                    ? "exact KSSD projection to common sparse drlevel"
                    : "direct coordinated KSSD comparison";
                return plan;

            case Rank::Backend::ProbMinHash:
                plan.method = Method::ProbMinHashNative;
                if (left.metadata.resolution != right.metadata.resolution) {
                    plan.reason =
                        "ProbMinHash comparison requires equal register counts";
                    return plan;
                }
                if (!probMinHashMetricSupported(
                        request.metric, left.metadata.weight_semantics)) {
                    plan.reason =
                        "weighted ProbMinHash does not estimate ordinary set/cardinality metrics";
                    return plan;
                }
                plan.support = Support::Supported;
                plan.common_resolution = left.metadata.resolution;
                plan.reason =
                    left.metadata.weight_semantics ==
                            Rank::WeightSemantics::UnweightedSet
                        ? "direct coordinated ProbMinHash comparison"
                        : "direct coordinated weighted ProbMinHash winner comparison";
                return plan;

            case Rank::Backend::BinDash:
                plan.method = Method::BinDashNative;
                if (left.metadata.resolution != right.metadata.resolution) {
                    plan.reason =
                        "BinDash comparison requires equal bin counts";
                    return plan;
                }
                plan.support = Support::Supported;
                plan.common_resolution = left.metadata.resolution;
                plan.reason = "direct coordinated BinDash comparison";
                return plan;

            case Rank::Backend::LegacyMinHash:
                plan.support = Support::Supported;
                plan.method = Method::LegacyMinHashNative;
                plan.common_resolution = std::min(
                    left.metadata.resolution, right.metadata.resolution);
                plan.exact_projection =
                    left.metadata.resolution != right.metadata.resolution;
                plan.reason = plan.exact_projection
                    ? "exact legacy bottom-K prefix projection to common K"
                    : "direct coordinated legacy MinHash comparison";
                return plan;

            case Rank::Backend::FracMinHash:
                plan.support = Support::Supported;
                plan.method = Method::FracMinHashNative;
                plan.common_resolution = std::max(
                    left.metadata.resolution, right.metadata.resolution);
                plan.exact_projection =
                    left.metadata.resolution != right.metadata.resolution;
                plan.reason = plan.exact_projection
                    ? "exact FracMinHash threshold projection to common scaled"
                    : "direct coordinated FracMinHash comparison";
                return plan;

            case Rank::Backend::OrderMinHash:
                plan.method = Method::OrderMinHashNative;
                if (!orderMetricSupported(request.metric)) {
                    plan.reason =
                        "OrderMinHash exposes order similarity/distance, not set metrics";
                    return plan;
                }
                plan.support = Support::Supported;
                plan.common_resolution = std::min(
                    left.metadata.resolution, right.metadata.resolution);
                plan.exact_projection =
                    left.metadata.resolution != right.metadata.resolution;
                plan.reason = plan.exact_projection
                    ? "exact OrderMinHash prefix comparison at common sample count"
                    : "direct coordinated OrderMinHash comparison";
                return plan;
        }
    }

    plan.reason = "cross-backend comparisons are not supported";
    return plan;
}

Result query(const FastKMV& left,
             const FastKMV& right,
             const Request& request) {
    Result result;
    result.plan = planQuery(describe(left), describe(right), request);
    if (!result.plan.executable()) return result;
    populateCardinalityMetrics(result, left.jaccard(right),
                               left.cardinality(), right.cardinality(),
                               left.metadata().kmer_size);
    result.confidence_level = request.confidence_level;
    result.effective_samples = std::min(left.size(), right.size());
    return result;
}

Result query(const FracMinHash& left,
             const FracMinHash& right,
             const Request& request) {
    Result result;
    result.plan = planQuery(describe(left), describe(right), request);
    if (!result.plan.executable()) return result;
    const FracMinHash projected_left = left.project(
        result.plan.common_resolution);
    const FracMinHash projected_right = right.project(
        result.plan.common_resolution);
    populateCardinalityMetrics(
        result, projected_left.jaccard(projected_right),
        projected_left.cardinality(), projected_right.cardinality(),
        left.metadata().kmer_size);
    result.confidence_level = request.confidence_level;
    result.effective_samples = std::min(
        projected_left.size(), projected_right.size());
    return result;
}

Result query(const BinDash& left,
             const BinDash& right,
             const Request& request) {
    Result result;
    result.plan = planQuery(describe(left), describe(right), request);
    if (!result.plan.executable()) return result;
    populateCardinalityMetrics(result, left.jaccard(right),
                               left.cardinality(), right.cardinality(),
                               left.metadata().kmer_size);
    result.confidence_level = request.confidence_level;
    result.effective_samples = result.plan.common_resolution;
    if (std::isfinite(result.jaccard)) {
        const double p_random = std::ldexp(
            1.0, -static_cast<int>(left.getBbits()));
        const double observed =
            p_random + (1.0 - p_random) * result.jaccard;
        result.plugin_standard_error =
            std::sqrt(std::max(0.0, observed * (1.0 - observed)) /
                      static_cast<double>(result.effective_samples)) /
            (1.0 - p_random);
        result.approximate_ci95_lower = std::max(
            0.0, result.jaccard - 1.96 * result.plugin_standard_error);
        result.approximate_ci95_upper = std::min(
            1.0, result.jaccard + 1.96 * result.plugin_standard_error);
        result.normal_ci95_eligible = result.effective_samples >= 30;
    }
    return result;
}

Result query(const HyperLogLog& left,
             const HyperLogLog& right,
             const Request& request) {
    Result result;
    result.plan = planQuery(describe(left), describe(right), request);
    if (!result.plan.executable()) return result;
    populateCardinalityMetrics(result, left.jaccard_index(right),
                               left.cardinality(), right.cardinality(),
                               left.metadata().kmer_size);
    result.confidence_level = request.confidence_level;
    result.effective_samples = UINT64_C(1) << result.plan.common_resolution;
    return result;
}

Result query(const SetSketch& left,
             const SetSketch& right,
             const Request& request) {
    Result result;
    result.plan = planQuery(describe(left), describe(right), request);
    if (!result.plan.executable()) return result;
    populateCardinalityMetrics(result, left.jaccard_index(right),
                               left.cardinality(), right.cardinality(),
                               left.metadata().kmer_size);
    result.confidence_level = request.confidence_level;
    result.effective_samples = UINT64_C(1) << result.plan.common_resolution;
    return result;
}

Result query(Kssd& left, Kssd& right, const Request& request) {
    Result result;
    result.plan = planQuery(describe(left), describe(right), request);
    if (!result.plan.executable()) return result;
    result.estimate_available = true;
    result.jaccard = left.jaccard(&right);
    result.confidence_level = request.confidence_level;
    populateDerivedMetrics(result, left.metadata().kmer_size);
    selectMetricValue(result, request.metric);
    return result;
}

Result query(MinHash& left,
             MinHash& right,
             const Request& request) {
    Result result;
    result.plan = planQuery(describe(left), describe(right), request);
    if (!result.plan.executable()) return result;
    left.finalize();
    right.finalize();
    const std::vector<uint64_t>& left_hashes = left.getHashesSorted();
    const std::vector<uint64_t>& right_hashes = right.getHashesSorted();
    populateCardinalityMetrics(result,
                               bottomKJaccard(
                                   left_hashes, right_hashes,
                                   result.plan.common_resolution),
                               left.cardinality(), right.cardinality(),
                               left.metadata().kmer_size);
    result.confidence_level = request.confidence_level;
    result.effective_samples = std::min(
        static_cast<uint64_t>(left.getSketchSize()),
        static_cast<uint64_t>(right.getSketchSize()));
    return result;
}

Result query(const OrderMinHash& left,
             const OrderMinHash& right,
             const Request& request) {
    Result result;
    result.plan = planQuery(describe(left), describe(right), request);
    if (!result.plan.executable()) return result;
    result.estimate_available = true;
    result.order_similarity = left.similarity(right);
    result.order_distance = 1.0 - result.order_similarity;
    result.value = request.metric == Metric::OrderDistance
        ? result.order_distance : result.order_similarity;
    result.confidence_level = request.confidence_level;
    result.effective_samples = result.plan.common_resolution;
    result.plugin_standard_error = std::sqrt(
        std::max(0.0, result.order_similarity *
                      (1.0 - result.order_similarity)) /
        static_cast<double>(result.effective_samples));
    result.approximate_ci95_lower = std::max(
        0.0, result.order_similarity - 1.96 * result.plugin_standard_error);
    result.approximate_ci95_upper = std::min(
        1.0, result.order_similarity + 1.96 * result.plugin_standard_error);
    result.normal_ci95_eligible = result.effective_samples >= 30;
    return result;
}

Result query(const ProbMinHash4& left,
             const ProbMinHash4& right,
             const Request& request) {
    Result result;
    result.plan = planQuery(describe(left), describe(right), request);
    if (!result.plan.executable()) return result;

    result.estimate_available = true;
    result.confidence_level = request.confidence_level;
    result.effective_samples = result.plan.common_resolution;
    result.left_weight = left.total_weight();
    result.right_weight = right.total_weight();

    const bool unweighted = left.weightSemantics() ==
        Rank::WeightSemantics::UnweightedSet;
    if (unweighted) {
        result.jaccard = left.jaccard(right);
        populateDerivedMetrics(result, left.metadata().kmer_size);
    }

    result.weighted_jaccard = left.jaccard_weighted(right);
    const double weight_sum = result.left_weight + result.right_weight;
    result.weighted_intersection = result.weighted_jaccard * weight_sum /
        (1.0 + result.weighted_jaccard);
    result.weighted_intersection = std::max(
        0.0, std::min(result.weighted_intersection,
                      std::min(result.left_weight, result.right_weight)));
    result.weighted_union = std::max(
        0.0, weight_sum - result.weighted_intersection);
    result.weighted_left_containment = result.left_weight > 0.0
        ? result.weighted_intersection / result.left_weight : 0.0;
    result.weighted_right_containment = result.right_weight > 0.0
        ? result.weighted_intersection / result.right_weight : 0.0;

    const double selected_jaccard = unweighted
        ? result.jaccard : result.weighted_jaccard;
    result.plugin_standard_error = std::sqrt(
        std::max(0.0, selected_jaccard * (1.0 - selected_jaccard)) /
        static_cast<double>(result.effective_samples));
    result.approximate_ci95_lower = std::max(
        0.0, selected_jaccard - 1.96 * result.plugin_standard_error);
    result.approximate_ci95_upper = std::min(
        1.0, selected_jaccard + 1.96 * result.plugin_standard_error);
    result.normal_ci95_eligible = result.effective_samples >= 30;
    if (!unweighted) {
        result.warnings.push_back(
            "left_weight/right_weight are accumulated input mass, not set cardinality");
    }
    selectMetricValue(result, request.metric);
    return result;
}

const char* metricName(Metric metric) noexcept {
    switch (metric) {
        case Metric::Jaccard: return "jaccard";
        case Metric::Intersection: return "intersection";
        case Metric::LeftContainment: return "left_containment";
        case Metric::RightContainment: return "right_containment";
        case Metric::ANI: return "ani_from_jaccard";
        case Metric::JaccardDistance: return "jaccard_distance";
        case Metric::MashDistance: return "mash_distance";
        case Metric::Union: return "union";
        case Metric::MaxContainment: return "max_containment";
        case Metric::AverageContainment: return "average_containment";
        case Metric::LeftContainmentANI: return "left_containment_ani";
        case Metric::RightContainmentANI: return "right_containment_ani";
        case Metric::MaxContainmentANI: return "max_containment_ani";
        case Metric::WeightedJaccard: return "weighted_jaccard";
        case Metric::WeightedContainment: return "weighted_containment";
        case Metric::OrderSimilarity: return "order_similarity";
        case Metric::OrderDistance: return "order_distance";
    }
    return "unknown";
}

const char* supportName(Support support) noexcept {
    switch (support) {
        case Support::Supported: return "supported";
        case Support::Unsupported: return "unsupported";
        case Support::Incompatible: return "incompatible";
    }
    return "unknown";
}

const char* methodName(Method method) noexcept {
    switch (method) {
        case Method::None: return "none";
        case Method::FastKmvNative: return "fastkmv_native";
        case Method::HllNative: return "hll_native";
        case Method::SetSketchNative: return "setsketch_native";
        case Method::KssdNative: return "kssd_native";
        case Method::ProbMinHashNative: return "probminhash_native";
        case Method::BinDashNative: return "bindash_native";
        case Method::LegacyMinHashNative: return "legacy_minhash_native";
        case Method::FracMinHashNative: return "fracminhash_native";
        case Method::OrderMinHashNative: return "orderminhash_native";
    }
    return "unknown";
}

} // namespace Query
} // namespace Sketch
