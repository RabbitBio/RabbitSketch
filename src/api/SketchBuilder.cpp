#include "api/SketchBuilder.h"

#include "BinDash.h"
#include "SetSketch.h"
#include "Sketch.h"
#include "fastkmv.h"
#include "fracminhash.h"
#include "probmh.h"
#include "rank/CanonicalKmer.h"
#include "rank/RankStream.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

namespace Sketch {
namespace API {

struct BuiltSketch::Impl {
    using State = std::variant<
        std::shared_ptr<MinHash>,
        std::shared_ptr<FastKMV>,
        std::shared_ptr<FracMinHash>,
        std::shared_ptr<HyperLogLog>,
        std::shared_ptr<SetSketch>,
        std::shared_ptr<Kssd>,
        std::shared_ptr<ProbMinHash4>,
        std::shared_ptr<BinDash>,
        std::shared_ptr<OrderMinHash>>;

    template <typename T>
    Impl(Algorithm algorithm_value, std::shared_ptr<T> state_value)
        : algorithm(algorithm_value), metadata(state_value->metadata()),
          state(std::move(state_value)) {}

    Algorithm algorithm;
    Rank::RankMetadata metadata;
    State state;
};

BuiltSketch::BuiltSketch(std::shared_ptr<Impl> impl) : impl_(std::move(impl)) {}

BuiltSketch BuiltSketch::fromLegacyMinHash(
    std::unique_ptr<MinHash> sketch) {
    if (!sketch) throw std::invalid_argument("MinHash sketch must not be null");
    std::shared_ptr<MinHash> shared(std::move(sketch));
    return BuiltSketch(std::make_shared<Impl>(
        Algorithm::LegacyMinHash, std::move(shared)));
}

BuiltSketch BuiltSketch::fromFastKMV(FastKMV sketch) {
    return BuiltSketch(std::make_shared<Impl>(
        Algorithm::FastKMV, std::make_shared<FastKMV>(std::move(sketch))));
}

BuiltSketch BuiltSketch::fromFracMinHash(FracMinHash sketch) {
    return BuiltSketch(std::make_shared<Impl>(
        Algorithm::FracMinHash,
        std::make_shared<FracMinHash>(std::move(sketch))));
}

BuiltSketch BuiltSketch::fromHyperLogLog(HyperLogLog sketch) {
    return BuiltSketch(std::make_shared<Impl>(
        Algorithm::HyperLogLog,
        std::make_shared<HyperLogLog>(std::move(sketch))));
}

BuiltSketch BuiltSketch::fromSetSketch(SetSketch sketch) {
    return BuiltSketch(std::make_shared<Impl>(
        Algorithm::SetSketch,
        std::make_shared<SetSketch>(std::move(sketch))));
}

BuiltSketch BuiltSketch::fromKssd(std::unique_ptr<Kssd> sketch) {
    if (!sketch) throw std::invalid_argument("KSSD sketch must not be null");
    std::shared_ptr<Kssd> shared(std::move(sketch));
    return BuiltSketch(std::make_shared<Impl>(
        Algorithm::KSSD, std::move(shared)));
}

BuiltSketch BuiltSketch::fromProbMinHash(ProbMinHash4 sketch) {
    return BuiltSketch(std::make_shared<Impl>(
        Algorithm::ProbMinHash,
        std::make_shared<ProbMinHash4>(std::move(sketch))));
}

BuiltSketch BuiltSketch::fromBinDash(BinDash sketch) {
    return BuiltSketch(std::make_shared<Impl>(
        Algorithm::BinDash, std::make_shared<BinDash>(std::move(sketch))));
}

BuiltSketch BuiltSketch::fromOrderMinHash(OrderMinHash sketch) {
    return BuiltSketch(std::make_shared<Impl>(
        Algorithm::OrderMinHash,
        std::make_shared<OrderMinHash>(std::move(sketch))));
}

Algorithm BuiltSketch::algorithm() const noexcept {
    return impl_->algorithm;
}

const Rank::RankMetadata& BuiltSketch::metadata() const noexcept {
    return impl_->metadata;
}

Query::SketchDescriptor BuiltSketch::descriptor() const {
    Query::SketchDescriptor result;
    result.metadata = impl_->metadata;
    return result;
}

Query::Result BuiltSketch::query(
    const BuiltSketch& other, const Query::Request& request) const {
    const auto rejected = [&]() {
        Query::Result result;
        result.plan = Query::planQuery(descriptor(), other.descriptor(), request);
        return result;
    };
    if (impl_->algorithm != other.impl_->algorithm) return rejected();

    return std::visit(
        [&](const auto& left) -> Query::Result {
            return std::visit(
                [&](const auto& right) -> Query::Result {
                    using Left = typename std::decay_t<decltype(left)>::element_type;
                    using Right = typename std::decay_t<decltype(right)>::element_type;
                    if constexpr (std::is_same<Left, Right>::value)
                        return Query::query(*left, *right, request);
                    return rejected();
                }, other.impl_->state);
        }, impl_->state);
}

namespace {

std::string cleanLabel(std::string value) {
    for (char& character : value)
        if (character == '\t' || character == '\r' || character == '\n' ||
            character == '\0')
            character = '_';
    return value;
}

std::string baseName(const std::string& path) {
    if (path == "-") return "stdin";
    const size_t slash = path.find_last_of("/\\");
    return cleanLabel(slash == std::string::npos
        ? path : path.substr(slash + 1));
}

std::string configToken(const SketchConfig& config) {
    std::ostringstream token;
    token << algorithmName(config.algorithm) << '-' << std::hex
          << std::setw(16) << std::setfill('0') << config.fingerprint();
    return token.str();
}

class Lane {
public:
    explicit Lane(SketchConfig config) : config_(config) {}
    virtual ~Lane() = default;
    virtual void updateSegment(const std::string& sequence) = 0;
    virtual BuiltSketch finish(BuildStats& stats) = 0;
    const SketchConfig& config() const noexcept { return config_; }
    BuildStats& stats() noexcept { return stats_; }
protected:
    SketchConfig config_;
    BuildStats stats_;
};

class LegacyMinHashLane final : public Lane {
public:
    explicit LegacyMinHashLane(const SketchConfig& config)
        : Lane(config), sketch_(new MinHash(
              config.kmer_size, config.resolution,
              checkedSeed(config.seed), config.canonical)) {}
    void updateSegment(const std::string& sequence) override {
        std::vector<char> mutable_sequence(sequence.begin(), sequence.end());
        mutable_sequence.push_back('\0');
        sketch_->update(mutable_sequence.data());
    }
    BuiltSketch finish(BuildStats&) override {
        sketch_->finalize();
        return BuiltSketch::fromLegacyMinHash(std::move(sketch_));
    }
private:
    static uint32_t checkedSeed(uint64_t seed) {
        if (seed > std::numeric_limits<uint32_t>::max())
            throw std::invalid_argument(
                "legacy MinHash seed must fit in 32 bits");
        return static_cast<uint32_t>(seed);
    }
    std::unique_ptr<MinHash> sketch_;
};

class FastKmvLane final : public Lane {
public:
    explicit FastKmvLane(const SketchConfig& config)
        : Lane(config), sketch_(config.resolution, config.kmer_size, config.seed) {}
    void updateSegment(const std::string& sequence) override {
        sketch_.update(sequence.data(), sequence.size());
    }
    BuiltSketch finish(BuildStats&) override {
        sketch_.finalize();
        return BuiltSketch::fromFastKMV(std::move(sketch_));
    }
private:
    FastKMV sketch_;
};

class FracMinHashLane final : public Lane {
public:
    explicit FracMinHashLane(const SketchConfig& config)
        : Lane(config), sketch_(static_cast<uint32_t>(config.scaled),
                                config.kmer_size, config.seed) {}
    void updateSegment(const std::string& sequence) override {
        sketch_.update(sequence.data(), sequence.size());
    }
    BuiltSketch finish(BuildStats&) override {
        sketch_.finalize();
        return BuiltSketch::fromFracMinHash(std::move(sketch_));
    }
private:
    FracMinHash sketch_;
};

class HllLane final : public Lane {
public:
    explicit HllLane(const SketchConfig& config)
        : Lane(config), sketch_(config.resolution, config.kmer_size,
                                config.seed, false) {}
    void updateSegment(const std::string& sequence) override {
        std::string mutable_sequence = sequence;
        sketch_.update(mutable_sequence.data());
    }
    BuiltSketch finish(BuildStats&) override {
        return BuiltSketch::fromHyperLogLog(std::move(sketch_));
    }
private:
    HyperLogLog sketch_;
};

class SetSketchLane final : public Lane {
public:
    explicit SetSketchLane(const SketchConfig& config)
        : Lane(config), sketch_(
              config.resolution, config.setsketch_base, config.setsketch_a,
              config.kmer_size, false, config.seed) {}
    void updateSegment(const std::string& sequence) override {
        std::string mutable_sequence = sequence;
        sketch_.update(mutable_sequence.data(), mutable_sequence.size());
    }
    BuiltSketch finish(BuildStats&) override {
        return BuiltSketch::fromSetSketch(std::move(sketch_));
    }
private:
    SetSketch sketch_;
};

class KssdLane final : public Lane {
public:
    explicit KssdLane(const SketchConfig& config)
        : Lane(config), parameters_(std::make_shared<kssd_parameter_t>(
              config.kmer_size / 2, config.kssd_half_subk,
              static_cast<int>(config.resolution))),
          sketch_(parameters_, config.seed) {}
    void updateSegment(const std::string& sequence) override {
        sketch_.update(sequence.c_str());
    }
    BuiltSketch finish(BuildStats&) override {
        sketch_.finalize();
        return BuiltSketch::fromKssd(
            std::unique_ptr<Kssd>(new Kssd(std::move(sketch_))));
    }
private:
    std::shared_ptr<const kssd_parameter_t> parameters_;
    Kssd sketch_;
};

class ProbMinHashLane final : public Lane {
public:
    explicit ProbMinHashLane(const SketchConfig& config)
        : Lane(config), sketch_(config.resolution, config.kmer_size,
              config.seed, config.probminhash_max_l,
              config.weight_semantics) {}
    void updateSegment(const std::string& sequence) override {
        if (!config_.track_abundance) {
            sketch_.update(sequence.data(), sequence.size());
            return;
        }
        Rank::CanonicalKmerIterator iterator(
            sequence.data(), sequence.size(), config_.kmer_size);
        const Rank::RankStream ranks(config_.seed);
        uint64_t canonical = 0;
        while (iterator.next(canonical))
            ++abundance_[ranks.fingerprint(canonical)];
    }
    BuiltSketch finish(BuildStats& stats) override {
        if (config_.track_abundance) {
            stats.distinct_kmers = abundance_.size();
            for (const auto& element : abundance_) {
                const uint64_t count = element.second;
                if (count < config_.min_abundance ||
                    (config_.max_abundance != 0 &&
                     count > config_.max_abundance)) {
                    ++stats.skipped_abundance_kmers;
                    continue;
                }
                sketch_.addHash(element.first, static_cast<double>(count));
                ++stats.accepted_distinct_kmers;
            }
            phmap::flat_hash_map<uint64_t, uint64_t>().swap(abundance_);
        }
        sketch_.finalize();
        return BuiltSketch::fromProbMinHash(std::move(sketch_));
    }
private:
    ProbMinHash4 sketch_;
    phmap::flat_hash_map<uint64_t, uint64_t> abundance_;
};

class BinDashLane final : public Lane {
public:
    explicit BinDashLane(const SketchConfig& config)
        : Lane(config), sketch_(config.resolution / 64,
              config.kmer_size, config.bindash_bits, config.seed) {}
    void updateSegment(const std::string& sequence) override {
        sketch_.update(sequence.data(), sequence.size());
    }
    BuiltSketch finish(BuildStats&) override {
        sketch_.finalize();
        return BuiltSketch::fromBinDash(std::move(sketch_));
    }
private:
    BinDash sketch_;
};

class OrderMinHashLane final : public Lane {
public:
    explicit OrderMinHashLane(const SketchConfig& config) : Lane(config) {}
    void updateSegment(const std::string& sequence) override {
        if (!sequence_.empty())
            throw std::logic_error(
                "OrderMinHash received more than one sequence segment");
        sequence_ = sequence;
    }
    BuiltSketch finish(BuildStats&) override {
        if (sequence_.empty())
            throw std::invalid_argument(
                "OrderMinHash record is shorter than its k-mer size");
        sketch_.setK(config_.kmer_size);
        sketch_.setL(config_.order_l);
        sketch_.setM(config_.order_m);
        sketch_.setSeed(config_.seed);
        sketch_.setReverseComplement(config_.canonical);
        sketch_.buildSketch(sequence_);
        return BuiltSketch::fromOrderMinHash(std::move(sketch_));
    }
private:
    OrderMinHash sketch_;
    std::string sequence_;
};

std::unique_ptr<Lane> makeLane(const SketchConfig& config) {
    config.validate();
    if (config.molecule != MoleculeType::DNA &&
        config.molecule != MoleculeType::RNA)
        throw std::invalid_argument(
            "the genomic FASTX builder currently accepts DNA or RNA");
    if (!config.canonical && config.algorithm != Algorithm::LegacyMinHash &&
        config.algorithm != Algorithm::OrderMinHash)
        throw std::invalid_argument(
            "forward-only mode is currently implemented only for legacy MinHash and OrderMinHash");
    switch (config.algorithm) {
        case Algorithm::LegacyMinHash:
            return std::unique_ptr<Lane>(new LegacyMinHashLane(config));
        case Algorithm::FastKMV:
            return std::unique_ptr<Lane>(new FastKmvLane(config));
        case Algorithm::FracMinHash:
            return std::unique_ptr<Lane>(new FracMinHashLane(config));
        case Algorithm::HyperLogLog:
            return std::unique_ptr<Lane>(new HllLane(config));
        case Algorithm::SetSketch:
            return std::unique_ptr<Lane>(new SetSketchLane(config));
        case Algorithm::KSSD:
            return std::unique_ptr<Lane>(new KssdLane(config));
        case Algorithm::ProbMinHash:
            return std::unique_ptr<Lane>(new ProbMinHashLane(config));
        case Algorithm::BinDash:
            return std::unique_ptr<Lane>(new BinDashLane(config));
        case Algorithm::OrderMinHash:
            if (config.aggregation != AggregationMode::OneSketchPerRecord)
                throw std::invalid_argument(
                    "OrderMinHash requires one-sketch-per-record aggregation");
            if (config.ambiguous_policy != AmbiguousPolicy::RejectRecord)
                throw std::invalid_argument(
                    "OrderMinHash requires RejectRecord ambiguous handling");
            if (config.minimum_base_quality >= 0)
                throw std::invalid_argument(
                    "OrderMinHash does not support gapped quality filtering");
            return std::unique_ptr<Lane>(new OrderMinHashLane(config));
    }
    throw std::invalid_argument("unknown sketch algorithm");
}

struct NormalizedRecord {
    std::string sequence;
    std::string quality;
};

NormalizedRecord normalizeRecord(const IO::FastxRecord& record,
                                 const SketchConfig& config) {
    if (config.minimum_base_quality >= 0 && record.quality.empty())
        throw std::invalid_argument(
            "minimum base quality was requested for a FASTA record: " +
            record.name);
    if (!record.quality.empty() &&
        record.quality.size() != record.sequence.size())
        throw std::invalid_argument(
            "FASTQ sequence/quality length mismatch: " + record.name);
    NormalizedRecord normalized;
    normalized.sequence.reserve(record.sequence.size());
    normalized.quality.reserve(record.quality.size());
    for (size_t index = 0; index < record.sequence.size(); ++index) {
        unsigned char raw = static_cast<unsigned char>(record.sequence[index]);
        char base = static_cast<char>(std::toupper(raw));
        if (config.molecule == MoleculeType::RNA && base == 'U') base = 'T';
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
            if (config.ambiguous_policy == AmbiguousPolicy::RejectRecord)
                throw std::invalid_argument(
                    "ambiguous base in rejected FASTX record: " + record.name);
            base = 'N';
        }
        char quality = 0;
        if (!record.quality.empty()) {
            quality = record.quality[index];
            const unsigned char encoded = static_cast<unsigned char>(quality);
            if (encoded < 33 || encoded > 126)
                throw std::invalid_argument(
                    "FASTQ quality is outside printable Phred+33 range: " +
                    record.name);
        }
        if (config.homopolymer_compressed && !normalized.sequence.empty() &&
            normalized.sequence.back() == base) {
            if (!normalized.quality.empty() && quality < normalized.quality.back())
                normalized.quality.back() = quality;
            continue;
        }
        normalized.sequence.push_back(base);
        if (!record.quality.empty()) normalized.quality.push_back(quality);
    }
    return normalized;
}

void updateLane(Lane& lane, const IO::FastxRecord& record) {
    const SketchConfig& config = lane.config();
    BuildStats& stats = lane.stats();
    const NormalizedRecord normalized = normalizeRecord(record, config);
    ++stats.records;
    stats.input_bases += record.sequence.size();
    const size_t length = normalized.sequence.size();
    const size_t k = config.kmer_size;
    if (length < k) return;
    const size_t window_count = length - k + 1;
    stats.candidate_kmers += window_count;

    size_t ambiguous = 0;
    size_t low_quality = 0;
    const auto low = [&](size_t index) {
        return !normalized.quality.empty() &&
            static_cast<int>(static_cast<unsigned char>(
                normalized.quality[index])) - 33 < config.minimum_base_quality;
    };
    for (size_t index = 0; index < k; ++index) {
        ambiguous += normalized.sequence[index] == 'N';
        low_quality += low(index);
    }
    for (size_t start = 0; start < window_count; ++start) {
        if (ambiguous != 0) ++stats.skipped_ambiguous_kmers;
        else if (low_quality != 0) ++stats.skipped_low_quality_kmers;
        else ++stats.accepted_kmers;
        if (start + k < length) {
            ambiguous -= normalized.sequence[start] == 'N';
            ambiguous += normalized.sequence[start + k] == 'N';
            low_quality -= low(start);
            low_quality += low(start + k);
        }
    }

    size_t segment_begin = 0;
    while (segment_begin < length) {
        while (segment_begin < length &&
               (normalized.sequence[segment_begin] == 'N' || low(segment_begin)))
            ++segment_begin;
        size_t segment_end = segment_begin;
        while (segment_end < length &&
               normalized.sequence[segment_end] != 'N' && !low(segment_end))
            ++segment_end;
        if (segment_end - segment_begin >= k)
            lane.updateSegment(normalized.sequence.substr(
                segment_begin, segment_end - segment_begin));
        segment_begin = segment_end + (segment_end < length ? 1 : 0);
    }
}

std::vector<SketchConfig> configsFor(
    const std::vector<SketchConfig>& configs, AggregationMode aggregation) {
    std::vector<SketchConfig> selected;
    for (const SketchConfig& config : configs)
        if (config.aggregation == aggregation) selected.push_back(config);
    return selected;
}

} // namespace

struct MultiSketchBuilder::Impl {
    std::vector<std::unique_ptr<Lane>> lanes;
    bool finished = false;
};

MultiSketchBuilder::MultiSketchBuilder(std::vector<SketchConfig> configs)
    : impl_(new Impl()) {
    if (configs.empty())
        throw std::invalid_argument(
            "MultiSketchBuilder requires at least one configuration");
    impl_->lanes.reserve(configs.size());
    for (const SketchConfig& config : configs)
        impl_->lanes.push_back(makeLane(config));
}

MultiSketchBuilder::~MultiSketchBuilder() = default;
MultiSketchBuilder::MultiSketchBuilder(MultiSketchBuilder&&) noexcept = default;
MultiSketchBuilder& MultiSketchBuilder::operator=(
    MultiSketchBuilder&&) noexcept = default;

void MultiSketchBuilder::update(const IO::FastxRecord& record) {
    if (!impl_ || impl_->finished)
        throw std::logic_error("MultiSketchBuilder update called after finish");
    for (std::unique_ptr<Lane>& lane : impl_->lanes)
        updateLane(*lane, record);
}

std::vector<BuildResult> MultiSketchBuilder::finish(
    const std::string& label_prefix,
    const std::string& source,
    const std::string& record_name) {
    if (!impl_ || impl_->finished)
        throw std::logic_error("MultiSketchBuilder finish called more than once");
    impl_->finished = true;
    std::vector<BuildResult> result;
    result.reserve(impl_->lanes.size());
    for (std::unique_ptr<Lane>& lane : impl_->lanes) {
        BuildStats stats = lane->stats();
        const std::string label = cleanLabel(label_prefix) + ':' +
            configToken(lane->config());
        BuiltSketch sketch = lane->finish(stats);
        stats.sealed = true;
        result.emplace_back(std::move(sketch), label, lane->config(), stats,
                            source, record_name);
    }
    return result;
}

size_t MultiSketchBuilder::lanes() const noexcept {
    return impl_ ? impl_->lanes.size() : 0;
}

bool MultiSketchBuilder::finished() const noexcept {
    return !impl_ || impl_->finished;
}

std::vector<BuildResult> buildFastxFiles(
    const std::vector<std::string>& paths,
    const std::vector<SketchConfig>& configs) {
    if (paths.empty())
        throw std::invalid_argument("FASTX build requires at least one path");
    if (configs.empty())
        throw std::invalid_argument("FASTX build requires at least one config");
    for (const SketchConfig& config : configs) config.validate();

    const std::vector<SketchConfig> collection_configs = configsFor(
        configs, AggregationMode::OneSketchPerCollection);
    const std::vector<SketchConfig> file_configs = configsFor(
        configs, AggregationMode::OneSketchPerFile);
    const std::vector<SketchConfig> record_configs = configsFor(
        configs, AggregationMode::OneSketchPerRecord);
    std::unique_ptr<MultiSketchBuilder> collection;
    if (!collection_configs.empty())
        collection.reset(new MultiSketchBuilder(collection_configs));

    std::vector<BuildResult> result;
    for (const std::string& path : paths) {
        // Record-level sketches are completed as the stream advances, but
        // buffer them until the file sketch is sealed.  This gives callers a
        // stable result order independent of parsing details:
        // file-level results, record-level results, then collection results.
        std::vector<BuildResult> record_results;
        std::unique_ptr<MultiSketchBuilder> file;
        if (!file_configs.empty())
            file.reset(new MultiSketchBuilder(file_configs));
        IO::FastxReader reader(path);
        IO::FastxRecord record;
        while (reader.next(record)) {
            if (file) file->update(record);
            if (collection) collection->update(record);
            if (!record_configs.empty()) {
                MultiSketchBuilder per_record(record_configs);
                per_record.update(record);
                std::vector<BuildResult> built = per_record.finish(
                    baseName(path) + ':' +
                        cleanLabel(record.name.empty()
                            ? "record-" + std::to_string(record.index)
                            : record.name) + '-' + std::to_string(record.index),
                    path, record.name);
                record_results.insert(record_results.end(),
                    std::make_move_iterator(built.begin()),
                    std::make_move_iterator(built.end()));
            }
        }
        if (file) {
            std::vector<BuildResult> built = file->finish(
                baseName(path), path, std::string());
            result.insert(result.end(),
                std::make_move_iterator(built.begin()),
                std::make_move_iterator(built.end()));
        }
        result.insert(result.end(),
            std::make_move_iterator(record_results.begin()),
            std::make_move_iterator(record_results.end()));
    }
    if (collection) {
        std::vector<BuildResult> built = collection->finish(
            "collection", std::string(), std::string());
        result.insert(result.end(),
            std::make_move_iterator(built.begin()),
            std::make_move_iterator(built.end()));
    }
    return result;
}

} // namespace API
} // namespace Sketch
