#ifndef RABBITSKETCH_PARTITIONED_SETSKETCH_H
#define RABBITSKETCH_PARTITIONED_SETSKETCH_H

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace Sketch { namespace detail {

// The native update assigns each item to ONE register. For Poisson load l,
// P(R <= k) = exp(-l * C[k]), C[k] = 1-exp(-a*b^-k).
// Rank zero includes empty registers; rank q+1 is the censored upper tail.
struct PartitionedSetSketchModel {
    std::array<double,64> c{}, delta{};
    uint32_t q;
    PartitionedSetSketchModel(double base, double a, uint32_t finite_q) : q(finite_q) {
        for (uint32_t k=0; k<=q; ++k)
            c[k] = -std::expm1(-a * std::pow(base, -static_cast<double>(k)));
        c[q+1] = 0.0;
        for (uint32_t k=1; k<=q+1; ++k) delta[k] = c[k-1]-c[k];
    }
    double logRank(uint32_t k, double load) const {
        if (k==0) return -load*c[0];
        return -load*c[k] + std::log(-std::expm1(-load*delta[k]));
    }
    double logDiagonal(uint32_t k, double load_a, double load_b, double shared) const {
        const double union_load = load_a + load_b - shared;
        if (k==0) return -union_load*c[0];
        const double d = delta[k];
        // Either a shared item reaches this rank, or both exclusive parts do.
        // This avoids cancellation in four nearly equal joint-CDF terms.
        const double shared_hit = -std::expm1(-shared*d);
        const double both_exclusive = std::exp(-shared*d) *
            (-std::expm1(-(load_a-shared)*d)) *
            (-std::expm1(-(load_b-shared)*d));
        return -union_load*c[k] + std::log(shared_hit+both_exclusive);
    }
    double logJoint(uint32_t rank_a, uint32_t rank_b,
                    double load_a, double load_b, double shared) const {
        if (rank_a==rank_b) return logDiagonal(rank_a,load_a,load_b,shared);
        if (rank_a<rank_b) return logRank(rank_a,load_a) + logRank(rank_b,load_b-shared);
        return logRank(rank_a,load_a-shared) + logRank(rank_b,load_b);
    }
};

}} // namespace Sketch::detail
#endif
