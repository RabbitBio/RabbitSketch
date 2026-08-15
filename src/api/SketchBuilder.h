#ifndef RABBITSKETCH_API_SKETCH_BUILDER_H
#define RABBITSKETCH_API_SKETCH_BUILDER_H

#include "api/SketchConfig.h"
#include "estimators/UnifiedQuery.h"
#include "io/FastxReader.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Sketch {

class BinDash;
class FastKMV;
class FracMinHash;
class HyperLogLog;
class Kssd;
class MinHash;
class OrderMinHash;
class ProbMinHash4;
class SetSketch;

namespace API {

/** Immutable type-erased sketch returned by the unified FASTX builder. */
class BuiltSketch {
public:
    static BuiltSketch fromLegacyMinHash(std::unique_ptr<MinHash> sketch);
    static BuiltSketch fromFastKMV(FastKMV sketch);
    static BuiltSketch fromFracMinHash(FracMinHash sketch);
    static BuiltSketch fromHyperLogLog(HyperLogLog sketch);
    static BuiltSketch fromSetSketch(SetSketch sketch);
    static BuiltSketch fromKssd(std::unique_ptr<Kssd> sketch);
    static BuiltSketch fromProbMinHash(ProbMinHash4 sketch);
    static BuiltSketch fromBinDash(BinDash sketch);
    static BuiltSketch fromOrderMinHash(OrderMinHash sketch);

    Algorithm algorithm() const noexcept;
    const Rank::RankMetadata& metadata() const noexcept;
    Query::SketchDescriptor descriptor() const;
    Query::Result query(
        const BuiltSketch& other,
        const Query::Request& request = Query::Request()) const;

private:
    struct Impl;
    explicit BuiltSketch(std::shared_ptr<Impl> impl);
    std::shared_ptr<Impl> impl_;
};

struct BuildResult {
    BuiltSketch sketch;
    std::string label;
    SketchConfig config;
    BuildStats stats;
    std::string source;
    std::string record_name;

    BuildResult(BuiltSketch sketch_value,
                std::string label_value,
                SketchConfig config_value,
                BuildStats stats_value,
                std::string source_value,
                std::string record_name_value)
        : sketch(std::move(sketch_value)), label(std::move(label_value)),
          config(config_value),
          stats(stats_value), source(std::move(source_value)),
          record_name(std::move(record_name_value)) {}
};

/** One-pass fan-out from normalized FASTX records into multiple algorithms. */
class MultiSketchBuilder {
public:
    explicit MultiSketchBuilder(std::vector<SketchConfig> configs);
    ~MultiSketchBuilder();
    MultiSketchBuilder(MultiSketchBuilder&&) noexcept;
    MultiSketchBuilder& operator=(MultiSketchBuilder&&) noexcept;
    MultiSketchBuilder(const MultiSketchBuilder&) = delete;
    MultiSketchBuilder& operator=(const MultiSketchBuilder&) = delete;

    void update(const IO::FastxRecord& record);
    std::vector<BuildResult> finish(
        const std::string& label_prefix,
        const std::string& source = std::string(),
        const std::string& record_name = std::string());
    size_t lanes() const noexcept;
    bool finished() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/** Build all requested aggregation modes while reading each input file once. */
std::vector<BuildResult> buildFastxFiles(
    const std::vector<std::string>& paths,
    const std::vector<SketchConfig>& configs);

inline std::vector<BuildResult> buildFastx(
    const std::string& path,
    const std::vector<SketchConfig>& configs) {
    return buildFastxFiles({path}, configs);
}

} // namespace API
} // namespace Sketch

#endif // RABBITSKETCH_API_SKETCH_BUILDER_H
