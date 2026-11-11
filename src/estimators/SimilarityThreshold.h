#ifndef RABBITSKETCH_SIMILARITY_THRESHOLD_H
#define RABBITSKETCH_SIMILARITY_THRESHOLD_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace Sketch { namespace detail {

inline void validateJaccardThreshold(double threshold) {
    if (!std::isfinite(threshold) || threshold < 0.0 || threshold > 1.0)
        throw std::invalid_argument("min_jaccard must be finite and in [0, 1]");
}

inline void validateDistanceThreshold(double threshold) {
    if (std::isnan(threshold) || threshold < 0.0)
        throw std::invalid_argument("max_distance must be nonnegative or +infinity");
}

inline double minimumJaccardForMashDistance(double maximum, int kmer_size) {
    validateDistanceThreshold(maximum);
    if (!std::isfinite(maximum)) return 0.0;
    const double kd = maximum * static_cast<double>(kmer_size);
    const double p = std::exp(-kd);
    if (p == 0.0) return 0.0;
    const double threshold = p / (2.0 - p);
    // Conservatively round the inverse transform: a pair at the distance
    // boundary must survive pruning. The caller still filters the fully
    // evaluated distance, so this small margin cannot admit extra output rows.
    const double margin = 32.0 * std::numeric_limits<double>::epsilon() * (1.0 + kd);
    return threshold * (1.0 - margin);
}

inline uint64_t minimumMatchesFixedDenominator(double threshold, uint64_t denominator) {
    if (threshold <= 0.0 || denominator == 0) return 0;
    uint64_t count = static_cast<uint64_t>(threshold * static_cast<double>(denominator));
    const auto qualifies = [&](uint64_t n) {
        return static_cast<double>(n) / static_cast<double>(denominator) >= threshold;
    };
    // ceil(threshold * denominator) can round an exact boundary up by one.
    // Adjust using the same division as the unfiltered estimator instead.
    while (count > 0 && qualifies(count - 1)) --count;
    while (count < denominator && !qualifies(count)) ++count;
    return count;
}

inline uint64_t minimumMatchesSetJaccard(double threshold, uint64_t size_sum) {
    if (threshold <= 0.0 || size_sum == 0) return 0;
    uint64_t count = static_cast<uint64_t>(
        (threshold / (1.0 + threshold)) * static_cast<double>(size_sum));
    const auto qualifies = [&](uint64_t n) {
        return n == size_sum ||
            static_cast<double>(n) / static_cast<double>(size_sum - n) >= threshold;
    };
    while (count > 0 && qualifies(count - 1)) --count;
    while (count < size_sum && !qualifies(count)) ++count;
    return count;
}

}} // namespace Sketch::detail
#endif
