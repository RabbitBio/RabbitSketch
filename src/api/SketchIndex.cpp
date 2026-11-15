#include "api/SketchIndex.h"
#include "InvertedIndex.h"

#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>

namespace Sketch { namespace API {
namespace {

struct Counter {
    explicit Counter(size_t size) : counts(size, 0), stamps(size, 0) {}
    std::vector<uint32_t> counts, stamps, touched;
    uint32_t epoch = 0;
    void reset() {
        touched.clear();
        if (++epoch == 0) {
            std::fill(stamps.begin(), stamps.end(), 0);
            epoch = 1;
        }
    }
    void add(uint32_t id) {
        if (stamps[id] != epoch) {
            stamps[id] = epoch;
            counts[id] = 1;
            touched.push_back(id);
        } else ++counts[id];
    }
};

void validateThreshold(double value) {
    if (!std::isfinite(value) || value < 0 || value > 1)
        throw std::invalid_argument("index similarity threshold must be finite and in [0, 1]");
}

struct Bucket {
    size_t operator()(const IndexKey& key) const noexcept {
        return static_cast<size_t>(Rank::RankStream::fmix64(key.value, key.position) >> 54);
    }
};

} // namespace

struct SketchIndex::Impl {
    std::vector<BuiltSketch> sketches;
    std::vector<std::vector<IndexKey>> keys;
    InvertedIndex<IndexKey> index;
    Query::Request request;
    bool approximate = false;

    void collect(uint32_t sample, Counter& counter, bool upper_only) const {
        if (sample >= sketches.size()) throw std::out_of_range("index sample id is out of range");
        counter.reset();
        for (const IndexKey& key : keys[sample]) {
            const auto found = index.postIdx.find(key);
            if (found == index.postIdx.end()) continue;
            const auto range = found->second;
            for (size_t k = range.off; k < range.off + range.cnt; ++k) {
                const uint32_t other = index.csrPosts[k];
                if (other == sample || (upper_only && other < sample)) continue;
                counter.add(other);
            }
        }
        std::sort(counter.touched.begin(), counter.touched.end());
    }

    std::vector<IndexMatch> verify(uint32_t sample, const Counter& counter,
                                   double minimum) const {
        std::vector<IndexMatch> matches;
        for (const uint32_t other : counter.touched) {
            auto result = sketches[sample].query(sketches[other], request);
            if (!result.estimate_available)
                throw std::logic_error("indexed query failed: " + result.plan.reason);
            if (result.value >= minimum)
                matches.push_back({sample, other, std::move(result)});
        }
        return matches;
    }
};

SketchIndex::SketchIndex(std::vector<BuiltSketch> sketches, int num_threads) {
    if (num_threads < 1) throw std::invalid_argument("index thread count must be positive");
    if (sketches.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("index collection exceeds supported sample count");
    auto state = std::make_shared<Impl>();
    state->sketches = std::move(sketches);
    if (!state->sketches.empty()) {
        const auto& reference = state->sketches.front();
        state->approximate = reference.algorithm() == Algorithm::HyperLogLog ||
                             reference.algorithm() == Algorithm::SetSketch;
        if (reference.algorithm() == Algorithm::ProbMinHash &&
            reference.metadata().weight_semantics != Rank::WeightSemantics::UnweightedSet)
            state->request.metric = Query::Metric::WeightedJaccard;
        state->keys.reserve(state->sketches.size());
        for (const auto& sketch : state->sketches) {
            reference.metadata().requireSameFamily(sketch.metadata(), "index_dict");
            const auto plan = Query::planQuery(reference.descriptor(), sketch.descriptor(), state->request);
            if (!plan.executable()) throw std::invalid_argument("index_dict: " + plan.reason);
            state->keys.push_back(sketch.indexKeys());
        }
        state->index = buildCSRIndexFromKeys<IndexKey>(
            state->keys, static_cast<int>(state->keys.size()), num_threads, Bucket{});
    }
    impl_ = std::move(state);
}

size_t SketchIndex::size() const noexcept { return impl_->sketches.size(); }
size_t SketchIndex::keyCount() const noexcept { return impl_->index.postIdx.size(); }
size_t SketchIndex::postingCount() const noexcept { return impl_->index.csrPosts.size(); }
bool SketchIndex::approximateCandidates() const noexcept { return impl_->approximate; }

std::vector<uint32_t> SketchIndex::postings(uint64_t value, uint32_t position) const {
    const auto found = impl_->index.postIdx.find({value, position});
    if (found == impl_->index.postIdx.end()) return {};
    const auto range = found->second;
    auto begin = impl_->index.csrPosts.begin() + range.off;
    std::vector<uint32_t> result(begin, begin + range.cnt);
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<IndexCandidate> SketchIndex::candidates(uint32_t sample) const {
    Counter counter(size());
    impl_->collect(sample, counter, false);
    std::vector<IndexCandidate> result;
    result.reserve(counter.touched.size());
    for (const uint32_t other : counter.touched)
        result.push_back({other, counter.counts[other]});
    return result;
}

std::vector<IndexMatch> SketchIndex::search(uint32_t sample, double min_similarity) const {
    validateThreshold(min_similarity);
    Counter counter(size());
    impl_->collect(sample, counter, false);
    return impl_->verify(sample, counter, min_similarity);
}

std::vector<IndexMatch> SketchIndex::pairs(double min_similarity, int num_threads) const {
    validateThreshold(min_similarity);
    if (num_threads < 1) throw std::invalid_argument("index thread count must be positive");
    if (size() == 0) return {};
    const int workers = std::min(num_threads, static_cast<int>(size()));
    std::vector<std::vector<IndexMatch>> rows(size());
    std::exception_ptr error;
    #pragma omp parallel num_threads(workers)
    {
        Counter counter(size());
        #pragma omp for schedule(dynamic, 1)
        for (size_t sample = 0; sample < size(); ++sample) {
            try {
                impl_->collect(static_cast<uint32_t>(sample), counter, true);
                rows[sample] = impl_->verify(static_cast<uint32_t>(sample), counter, min_similarity);
            } catch (...) {
                #pragma omp critical(rabbitsketch_index_error)
                { if (!error) error = std::current_exception(); }
            }
        }
    }
    if (error) std::rethrow_exception(error);
    std::vector<IndexMatch> result;
    size_t count = 0;
    for (const auto& row : rows) count += row.size();
    result.reserve(count);
    for (auto& row : rows)
        for (auto& match : row) result.push_back(std::move(match));
    return result;
}

SketchIndex indexDict(const std::vector<BuiltSketch>& sketches, int num_threads) {
    return SketchIndex(sketches, num_threads);
}

SketchIndex BuiltSketch::indexDict(const std::vector<BuiltSketch>& sketches, int num_threads) {
    return API::indexDict(sketches, num_threads);
}

}} // namespace Sketch::API
