#ifndef RABBITSKETCH_UNIFIED_QUERY_H
#define RABBITSKETCH_UNIFIED_QUERY_H

#include "rank/RankMetadata.h"

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace Sketch {

class FastKMV;
class FracMinHash;
class BinDash;
class HyperLogLog;
class Kssd;
class MinHash;
class OrderMinHash;
class ProbMinHash4;
class SetSketch;

namespace Query {

/** Common pairwise estimands exposed by the RabbitSketch2 query layer. */
enum class Metric : uint8_t {
    Jaccard = 1,
    Intersection = 2,
    LeftContainment = 3,
    RightContainment = 4,
    /** Historical RabbitSketch2 ANI: ANI derived from Jaccard. */
    ANI = 5,
    JaccardDistance = 6,
    MashDistance = 7,
    Union = 8,
    MaxContainment = 9,
    AverageContainment = 10,
    LeftContainmentANI = 11,
    RightContainmentANI = 12,
    MaxContainmentANI = 13,
    WeightedJaccard = 14,
    WeightedContainment = 15,
    OrderSimilarity = 16,
    OrderDistance = 17
};

/** Scientific status of a proposed query, independent of implementation speed. */
enum class Support : uint8_t {
    Supported = 1,
    Unsupported = 2,
    Incompatible = 3
};

/** Concrete execution path selected by the planner. */
enum class Method : uint8_t {
    None = 0,
    FastKmvNative = 1,
    HllNative = 2,
    SetSketchNative = 3,
    KssdNative = 4,
    ProbMinHashNative = 5,
    BinDashNative = 6,
    LegacyMinHashNative = 7,
    FracMinHashNative = 8,
    OrderMinHashNative = 9
};

/** Backend-independent query policy. */
struct Request {
    Metric metric = Metric::Jaccard;
    double confidence_level = 0.95;
};

/** Metadata plus the payload capability needed for planning without sequences. */
struct SketchDescriptor {
    Rank::RankMetadata metadata;
};

struct Plan {
    Metric metric = Metric::Jaccard;
    Support support = Support::Unsupported;
    Method method = Method::None;
    Rank::Backend left_backend = Rank::Backend::FastKMV;
    Rank::Backend right_backend = Rank::Backend::FastKMV;
    uint32_t common_resolution = 0;
    bool exact_projection = false;
    std::string reason;

    bool executable() const noexcept {
        return support == Support::Supported;
    }
};

/**
 * Backend-independent result schema.  Unavailable quantities remain NaN and
 * estimate_available is false when the planner rejects a query.
 */
struct Result {
    Plan plan;
    bool estimate_available = false;
    double value = std::numeric_limits<double>::quiet_NaN();
    double jaccard = std::numeric_limits<double>::quiet_NaN();
    double jaccard_distance = std::numeric_limits<double>::quiet_NaN();
    double mash_distance = std::numeric_limits<double>::quiet_NaN();
    double ani = std::numeric_limits<double>::quiet_NaN();
    double intersection = std::numeric_limits<double>::quiet_NaN();
    double union_size = std::numeric_limits<double>::quiet_NaN();
    double left_containment = std::numeric_limits<double>::quiet_NaN();
    double right_containment = std::numeric_limits<double>::quiet_NaN();
    double max_containment = std::numeric_limits<double>::quiet_NaN();
    double average_containment = std::numeric_limits<double>::quiet_NaN();
    double left_containment_ani = std::numeric_limits<double>::quiet_NaN();
    double right_containment_ani = std::numeric_limits<double>::quiet_NaN();
    double max_containment_ani = std::numeric_limits<double>::quiet_NaN();
    double weighted_jaccard = std::numeric_limits<double>::quiet_NaN();
    double weighted_left_containment =
        std::numeric_limits<double>::quiet_NaN();
    double weighted_right_containment =
        std::numeric_limits<double>::quiet_NaN();
    double weighted_intersection =
        std::numeric_limits<double>::quiet_NaN();
    double weighted_union = std::numeric_limits<double>::quiet_NaN();
    double order_similarity = std::numeric_limits<double>::quiet_NaN();
    double order_distance = std::numeric_limits<double>::quiet_NaN();
    double left_weight = std::numeric_limits<double>::quiet_NaN();
    double right_weight = std::numeric_limits<double>::quiet_NaN();
    double left_cardinality = std::numeric_limits<double>::quiet_NaN();
    double right_cardinality = std::numeric_limits<double>::quiet_NaN();
    uint64_t effective_samples = 0;
    double confidence_level = 0.95;
    double plugin_standard_error = std::numeric_limits<double>::quiet_NaN();
    double approximate_ci95_lower = std::numeric_limits<double>::quiet_NaN();
    double approximate_ci95_upper = std::numeric_limits<double>::quiet_NaN();
    bool normal_ci95_eligible = false;
    std::vector<std::string> warnings;
};

SketchDescriptor describe(const FastKMV& sketch);
SketchDescriptor describe(const FracMinHash& sketch);
SketchDescriptor describe(const BinDash& sketch);
SketchDescriptor describe(const HyperLogLog& sketch);
SketchDescriptor describe(const SetSketch& sketch);
SketchDescriptor describe(const Kssd& sketch);
SketchDescriptor describe(const MinHash& sketch);
SketchDescriptor describe(const OrderMinHash& sketch);
SketchDescriptor describe(const ProbMinHash4& sketch);

Plan planQuery(const SketchDescriptor& left,
               const SketchDescriptor& right,
               const Request& request = Request());

Result query(const FastKMV& left,
             const FastKMV& right,
             const Request& request = Request());
Result query(const FracMinHash& left,
             const FracMinHash& right,
             const Request& request = Request());
Result query(const BinDash& left,
             const BinDash& right,
             const Request& request = Request());
Result query(const HyperLogLog& left,
             const HyperLogLog& right,
             const Request& request = Request());
Result query(const SetSketch& left,
             const SetSketch& right,
             const Request& request = Request());
Result query(Kssd& left,
             Kssd& right,
             const Request& request = Request());
Result query(MinHash& left,
             MinHash& right,
             const Request& request = Request());
Result query(const OrderMinHash& left,
             const OrderMinHash& right,
             const Request& request = Request());
Result query(const ProbMinHash4& left,
             const ProbMinHash4& right,
             const Request& request = Request());

const char* metricName(Metric metric) noexcept;
const char* supportName(Support support) noexcept;
const char* methodName(Method method) noexcept;

} // namespace Query
} // namespace Sketch

#endif // RABBITSKETCH_UNIFIED_QUERY_H
