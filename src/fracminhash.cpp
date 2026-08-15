#include "fracminhash.h"

#include "api/RuntimeInfo.h"
#include "rank/CanonicalKmer.h"
#include "rank/RankStream.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace Sketch {
namespace {

uint64_t thresholdFor(uint32_t scaled) noexcept {
    return scaled == 1 ? UINT64_MAX : UINT64_MAX / scaled;
}

double samplingProbabilityFor(uint32_t scaled) noexcept {
    if (scaled == 1) return 1.0;
    return std::ldexp(static_cast<double>(thresholdFor(scaled)) + 1.0, -64);
}

constexpr uint32_t BUFFERED_BUILD_MAX_SCALED = 64;
constexpr uint64_t MAX_SPECULATIVE_RESERVE = UINT64_C(8) * 1024 * 1024;

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
    result.common = Runtime::countSortedIntersectionU64(
        left.data(), static_cast<size_t>(result.left),
        right.data(), static_cast<size_t>(result.right));
    return result;
}

} // namespace

FracMinHash::FracMinHash(uint32_t scaled, int kmer_size, uint64_t seed)
    : scaled_(scaled), kmer_size_(kmer_size), seed_(seed),
      buffered_build_(scaled <= BUFFERED_BUILD_MAX_SCALED) {
    if (scaled_ == 0)
        throw std::invalid_argument("FracMinHash scaled must be positive");
    if (kmer_size_ < 1 || kmer_size_ > 32)
        throw std::invalid_argument(
            "FracMinHash k-mer size must be in [1, 32]");
}

FracMinHash FracMinHash::fromHashes(
    uint32_t scaled, int kmer_size, uint64_t seed,
    std::vector<uint64_t> hashes) {
    FracMinHash result(scaled, kmer_size, seed);
    if (!std::is_sorted(hashes.begin(), hashes.end()) ||
        std::adjacent_find(hashes.begin(), hashes.end()) != hashes.end())
        throw std::invalid_argument(
            "FracMinHash restored hashes must be sorted and distinct");
    if (!hashes.empty() && hashes.back() > result.threshold())
        throw std::invalid_argument(
            "FracMinHash restored hash exceeds the scaled threshold");
    result.hashes_ = std::move(hashes);
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
    if (length < static_cast<uint64_t>(kmer_size_)) return;
    reserveForUpdate(length);

    Rank::CanonicalKmerIterator iterator(
        sequence, static_cast<size_t>(length),
        static_cast<uint8_t>(kmer_size_));
    const uint64_t admission_threshold = threshold();
    uint64_t canonical[8] = {};
    uint64_t fingerprints[8] = {};
    bool have_last = false;
    uint64_t last = 0;

    for (;;) {
        unsigned count = 0;
        while (count < 8 && iterator.next(canonical[count])) ++count;
        if (count == 0) break;
        const uint8_t valid_mask = count == 8
            ? UINT8_MAX
            : static_cast<uint8_t>((UINT16_C(1) << count) - 1);
        uint8_t accepted = Runtime::filterFmix64x8(
            canonical, valid_mask, seed_, admission_threshold, fingerprints);
        while (accepted != 0) {
            unsigned lane = 0;
            while ((accepted & static_cast<uint8_t>(UINT8_C(1) << lane)) == 0)
                ++lane;
            accepted &= static_cast<uint8_t>(accepted - 1);
            const uint64_t fingerprint = fingerprints[lane];
            if (!have_last || fingerprint != last)
                retainFingerprint(fingerprint);
            last = fingerprint;
            have_last = true;
        }
        if (count != 8) break;
    }
}

void FracMinHash::addHash(uint64_t coordinated_fingerprint) {
    requireMutable("FracMinHash::addHash");
    if (coordinated_fingerprint <= threshold())
        retainFingerprint(coordinated_fingerprint);
}

void FracMinHash::reserveForUpdate(uint64_t length) {
    const uint64_t windows = length - static_cast<uint64_t>(kmer_size_) + 1;
    uint64_t expected = windows / scaled_ + (windows % scaled_ != 0);
    const uint64_t headroom = expected / 16 + 64;
    if (expected <= UINT64_MAX - headroom) expected += headroom;
    expected = std::min(expected, MAX_SPECULATIVE_RESERVE);

    const size_t current = buffered_build_
        ? build_buffer_.size() : build_hashes_.size();
    const size_t increment = static_cast<size_t>(std::min<uint64_t>(
        expected, static_cast<uint64_t>(
            std::numeric_limits<size_t>::max() - current)));
    const size_t target = current + increment;
    if (buffered_build_) {
        if (build_buffer_.capacity() < target) build_buffer_.reserve(target);
    } else if (build_hashes_.capacity() < target) {
        build_hashes_.reserve(target);
    }
}

void FracMinHash::retainFingerprint(uint64_t fingerprint) {
    if (buffered_build_) build_buffer_.push_back(fingerprint);
    else build_hashes_.insert(fingerprint);
}

void FracMinHash::finalize() {
    if (sealed_) return;
    if (buffered_build_) {
        std::sort(build_buffer_.begin(), build_buffer_.end());
        build_buffer_.erase(
            std::unique(build_buffer_.begin(), build_buffer_.end()),
            build_buffer_.end());
        hashes_ = std::move(build_buffer_);
    } else {
        hashes_.assign(build_hashes_.begin(), build_hashes_.end());
        std::sort(hashes_.begin(), hashes_.end());
    }
    phmap::flat_hash_set<uint64_t>().swap(build_hashes_);
    std::vector<uint64_t>().swap(build_buffer_);
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
    return samplingProbabilityFor(scaled_);
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
    const std::vector<uint64_t>& left = hashes();
    const std::vector<uint64_t>& right = other.hashes();
    const uint64_t common_threshold = thresholdFor(common_scaled);
    const auto left_end = std::upper_bound(
        left.begin(), left.end(), common_threshold);
    const auto right_end = std::upper_bound(
        right.begin(), right.end(), common_threshold);
    std::vector<uint64_t> merged;
    merged.reserve(static_cast<size_t>(left_end - left.begin()) +
                   static_cast<size_t>(right_end - right.begin()));
    std::set_union(left.begin(), left_end, right.begin(), right_end,
                   std::back_inserter(merged));
    return fromHashes(common_scaled, kmer_size_, seed_, std::move(merged));
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
    return cardinalityAt(scaled_);
}

double FracMinHash::cardinalityAt(uint32_t target_scaled) const {
    const size_t retained = target_scaled == scaled_
        ? size() : retainedAt(target_scaled);
    return static_cast<double>(retained) /
        samplingProbabilityFor(target_scaled);
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
