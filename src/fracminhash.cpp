#include "fracminhash.h"

#include "rank/CanonicalKmer.h"
#include "rank/RankStream.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>

namespace Sketch {
namespace {

uint64_t thresholdFor(uint32_t scaled) noexcept {
    return scaled == 1 ? UINT64_MAX : UINT64_MAX / scaled;
}

struct IntersectionStats {
    uint64_t common = 0;
    uint64_t left = 0;
    uint64_t right = 0;
};

IntersectionStats compareAtThreshold(const std::vector<uint64_t>& left,
                                     const std::vector<uint64_t>& right,
                                     uint64_t threshold) {
    IntersectionStats result;
    result.left = static_cast<uint64_t>(
        std::upper_bound(left.begin(), left.end(), threshold) - left.begin());
    result.right = static_cast<uint64_t>(
        std::upper_bound(right.begin(), right.end(), threshold) - right.begin());
    size_t left_index = 0;
    size_t right_index = 0;
    while (left_index < result.left && right_index < result.right) {
        if (left[left_index] < right[right_index]) ++left_index;
        else if (right[right_index] < left[left_index]) ++right_index;
        else {
            ++result.common;
            ++left_index;
            ++right_index;
        }
    }
    return result;
}

} // namespace

FracMinHash::FracMinHash(uint32_t scaled, int kmer_size, uint64_t seed)
    : scaled_(scaled), kmer_size_(kmer_size), seed_(seed) {
    if (scaled_ == 0)
        throw std::invalid_argument("FracMinHash scaled must be positive");
    if (kmer_size_ < 1 || kmer_size_ > 32)
        throw std::invalid_argument(
            "FracMinHash k-mer size must be in [1, 32]");
}

FracMinHash FracMinHash::fromHashes(
    uint32_t scaled, int kmer_size, uint64_t seed,
    const std::vector<uint64_t>& hashes) {
    FracMinHash result(scaled, kmer_size, seed);
    if (!std::is_sorted(hashes.begin(), hashes.end()) ||
        std::adjacent_find(hashes.begin(), hashes.end()) != hashes.end())
        throw std::invalid_argument(
            "FracMinHash restored hashes must be sorted and distinct");
    if (!hashes.empty() && hashes.back() > result.threshold())
        throw std::invalid_argument(
            "FracMinHash restored hash exceeds the scaled threshold");
    result.hashes_ = hashes;
    result.sealed_ = true;
    return result;
}

void FracMinHash::requireMutable(const char* operation) const {
    if (sealed_)
        throw std::logic_error(std::string(operation) +
                               " called after finalize");
}

void FracMinHash::update(const char* sequence, uint64_t length) {
    requireMutable("FracMinHash::update");
    if (sequence == nullptr && length != 0)
        throw std::invalid_argument(
            "FracMinHash sequence is null with nonzero length");
    Rank::CanonicalKmerIterator iterator(
        sequence, static_cast<size_t>(length),
        static_cast<uint8_t>(kmer_size_));
    const Rank::RankStream ranks(seed_);
    uint64_t canonical = 0;
    while (iterator.next(canonical))
        addHash(ranks.fingerprint(canonical));
}

void FracMinHash::addHash(uint64_t coordinated_fingerprint) {
    requireMutable("FracMinHash::addHash");
    if (coordinated_fingerprint <= threshold())
        build_hashes_.insert(coordinated_fingerprint);
}

void FracMinHash::finalize() {
    if (sealed_) return;
    hashes_.assign(build_hashes_.begin(), build_hashes_.end());
    std::sort(hashes_.begin(), hashes_.end());
    phmap::flat_hash_set<uint64_t>().swap(build_hashes_);
    hashes_.shrink_to_fit();
    sealed_ = true;
}

void FracMinHash::ensureFinalized() const {
    const_cast<FracMinHash*>(this)->finalize();
}

void FracMinHash::requireCompatible(const FracMinHash& other,
                                    const char* operation) const {
    metadata().requireSameFamily(other.metadata(), operation, true);
}

uint64_t FracMinHash::threshold() const noexcept {
    return thresholdFor(scaled_);
}

double FracMinHash::samplingProbability() const noexcept {
    if (scaled_ == 1) return 1.0;
    return std::ldexp(static_cast<double>(threshold()) + 1.0, -64);
}

const std::vector<uint64_t>& FracMinHash::hashes() const {
    ensureFinalized();
    return hashes_;
}

size_t FracMinHash::size() const {
    return hashes().size();
}

size_t FracMinHash::retainedAt(uint32_t target_scaled) const {
    if (target_scaled < scaled_)
        throw std::invalid_argument(
            "FracMinHash cannot project to a denser scaled value");
    const std::vector<uint64_t>& values = hashes();
    return static_cast<size_t>(std::upper_bound(
        values.begin(), values.end(), thresholdFor(target_scaled)) -
        values.begin());
}

FracMinHash FracMinHash::project(uint32_t target_scaled) const {
    if (target_scaled < scaled_)
        throw std::invalid_argument(
            "FracMinHash cannot project to a denser scaled value");
    const std::vector<uint64_t>& values = hashes();
    const auto end = std::upper_bound(
        values.begin(), values.end(), thresholdFor(target_scaled));
    return fromHashes(target_scaled, kmer_size_, seed_,
                      std::vector<uint64_t>(values.begin(), end));
}

FracMinHash FracMinHash::merge(const FracMinHash& other) const {
    requireCompatible(other, "FracMinHash::merge");
    const uint32_t common_scaled = std::max(scaled_, other.scaled_);
    const FracMinHash left = project(common_scaled);
    const FracMinHash right = other.project(common_scaled);
    std::vector<uint64_t> merged;
    merged.reserve(left.hashes_.size() + right.hashes_.size());
    std::set_union(left.hashes_.begin(), left.hashes_.end(),
                   right.hashes_.begin(), right.hashes_.end(),
                   std::back_inserter(merged));
    return fromHashes(common_scaled, kmer_size_, seed_, merged);
}

double FracMinHash::jaccard(const FracMinHash& other) const {
    requireCompatible(other, "FracMinHash::jaccard");
    const uint64_t common_threshold = thresholdFor(
        std::max(scaled_, other.scaled_));
    const IntersectionStats stats = compareAtThreshold(
        hashes(), other.hashes(), common_threshold);
    const uint64_t union_size = stats.left + stats.right - stats.common;
    return union_size == 0
        ? 1.0 : static_cast<double>(stats.common) / union_size;
}

double FracMinHash::containment(const FracMinHash& other) const {
    requireCompatible(other, "FracMinHash::containment");
    const uint64_t common_threshold = thresholdFor(
        std::max(scaled_, other.scaled_));
    const IntersectionStats stats = compareAtThreshold(
        hashes(), other.hashes(), common_threshold);
    return stats.left == 0
        ? 1.0 : static_cast<double>(stats.common) / stats.left;
}

double FracMinHash::cardinality() const {
    return static_cast<double>(size()) / samplingProbability();
}

double FracMinHash::distance(const FracMinHash& other) const {
    const double similarity = jaccard(other);
    if (!(similarity > 0.0))
        return std::numeric_limits<double>::infinity();
    if (similarity >= 1.0) return 0.0;
    return -std::log(2.0 * similarity / (1.0 + similarity)) /
        static_cast<double>(kmer_size_);
}

double FracMinHash::ani(const FracMinHash& other) const {
    const double mash_distance = distance(other);
    if (!std::isfinite(mash_distance)) return 0.0;
    return std::exp(-mash_distance);
}

size_t FracMinHash::memoryBytes() const {
    ensureFinalized();
    return hashes_.capacity() * sizeof(uint64_t);
}

Rank::RankMetadata FracMinHash::metadata() const {
    Rank::RankMetadata meta;
    meta.hash_profile = Rank::HashProfile::FracRankStreamV1;
    meta.backend = Rank::Backend::FracMinHash;
    meta.weight_semantics = Rank::WeightSemantics::UnweightedSet;
    meta.resolution_kind = Rank::ResolutionKind::Scaled;
    meta.fingerprint_bits = 64;
    meta.kmer_size = static_cast<uint16_t>(kmer_size_);
    meta.resolution = scaled_;
    meta.seed = seed_;
    meta.parameter_fingerprint = Rank::RankStream::fmix64(
        1, UINT64_C(0x667261636d696e68));
    return meta;
}

} // namespace Sketch
