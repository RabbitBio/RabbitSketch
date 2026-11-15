#ifndef RABBITSKETCH_API_SKETCH_INDEX_H
#define RABBITSKETCH_API_SKETCH_INDEX_H

#include "api/SketchBuilder.h"

namespace Sketch { namespace API {

struct IndexCandidate {
    uint32_t sample = 0;
    uint32_t shared_keys = 0;
};

struct IndexMatch {
    uint32_t left = 0;
    uint32_t right = 0;
    Query::Result estimate;
};

/** CSR index for comparisons within one compatible collection.
 * Singletons are discarded. HLL/SetSketch witness candidate generation is
 * approximate; every returned match is verified by the native estimator.
 */
class SketchIndex {
public:
    explicit SketchIndex(std::vector<BuiltSketch> sketches, int num_threads = 1);
    size_t size() const noexcept;
    size_t keyCount() const noexcept;
    size_t postingCount() const noexcept;
    bool approximateCandidates() const noexcept;
    std::vector<uint32_t> postings(uint64_t value, uint32_t position = 0) const;
    std::vector<IndexCandidate> candidates(uint32_t sample) const;
    std::vector<IndexMatch> search(uint32_t sample, double min_similarity = 0.0) const;
    std::vector<IndexMatch> pairs(double min_similarity = 0.0, int num_threads = 1) const;
private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};

SketchIndex indexDict(const std::vector<BuiltSketch>& sketches, int num_threads = 1);

}} // namespace Sketch::API
#endif
