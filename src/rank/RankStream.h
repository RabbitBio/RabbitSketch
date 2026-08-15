#ifndef RABBITSKETCH_RANK_STREAM_H
#define RABBITSKETCH_RANK_STREAM_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace Sketch {
namespace Rank {

/**
 * Versioned integer random stream used by RabbitSketch2 backends.
 *
 * The element fingerprint deliberately matches the historical RabbitSketch
 * single-round murmur3_fmix(canonical_2bit_kmer, seed) pipeline.  Keeping this
 * operation separate from derive() is important: fingerprints are coordinated
 * across backends, while admission/register/race variables use explicit domain
 * separation and can therefore evolve without silently changing identity.
 */
class RankStream {
public:
    static constexpr uint16_t VERSION = 1;

    enum class Domain : uint64_t {
        Admission    = UINT64_C(0x243f6a8885a308d3),
        Bucket       = UINT64_C(0x13198a2e03707344),
        Register     = UINT64_C(0xa4093822299f31d0),
        WeightedRace = UINT64_C(0x082efa98ec4e6c89)
    };

    explicit constexpr RankStream(uint64_t seed = 42) noexcept : seed_(seed) {}

    constexpr uint64_t seed() const noexcept { return seed_; }

    static inline uint64_t fmix64(uint64_t value, uint64_t seed) noexcept {
        value ^= seed;
        value ^= value >> 33;
        value *= UINT64_C(0xff51afd7ed558ccd);
        value ^= value >> 33;
        value *= UINT64_C(0xc4ceb9fe1a85ec53);
        value ^= value >> 33;
        return value;
    }

    /** Stable, backend-independent identity for a canonical k-mer. */
    inline uint64_t fingerprint(uint64_t canonical_code) const noexcept {
        return fmix64(canonical_code, seed_);
    }

    /**
     * Derive a deterministic substream value from an element fingerprint.
     * A counter identifies a register, retry, or independent lane.
     */
    inline uint64_t derive(uint64_t fingerprint_value,
                           Domain domain,
                           uint64_t counter = 0) const noexcept {
        const uint64_t step = UINT64_C(0x9e3779b97f4a7c15) * (counter + 1);
        const uint64_t domain_key = static_cast<uint64_t>(domain);
        return fmix64(fingerprint_value ^ domain_key ^ step,
                      seed_ ^ UINT64_C(0xd1b54a32d192ed03));
    }

    /** Convert an integer rank to the open interval (0, 1). */
    static inline double uniformOpen01(uint64_t rank) noexcept {
        constexpr double INV_TWO_POW_53 = 1.0 / 9007199254740992.0;
        return (static_cast<double>(rank >> 11) + 0.5) * INV_TWO_POW_53;
    }

    inline double exponentialRank(uint64_t fingerprint_value,
                                  uint64_t counter,
                                  double weight) const {
        if (!(weight > 0.0) || !std::isfinite(weight)) {
            throw std::invalid_argument(
                "RankStream::exponentialRank requires a finite positive weight");
        }
        const double u = uniformOpen01(
            derive(fingerprint_value, Domain::WeightedRace, counter));
        return -std::log(u) / weight;
    }

    /** Top precision bits used by a partitioned register sketch. */
    static inline uint32_t bucket(uint64_t fingerprint_value,
                                  uint8_t precision) {
        if (precision == 0 || precision >= 32) {
            throw std::invalid_argument("bucket precision must be in [1, 31]");
        }
        return static_cast<uint32_t>(fingerprint_value >> (64 - precision));
    }

    /** HLL rho value after removing the top precision bucket bits. */
    static inline uint8_t hllRank(uint64_t fingerprint_value,
                                  uint8_t precision) {
        if (precision == 0 || precision >= 64) {
            throw std::invalid_argument("HLL precision must be in [1, 63]");
        }
        const uint64_t remainder = fingerprint_value << precision;
        if (remainder == 0) {
            return static_cast<uint8_t>(65 - precision);
        }
#if defined(__GNUC__) || defined(__clang__)
        return static_cast<uint8_t>(__builtin_clzll(remainder) + 1);
#else
        uint8_t zeros = 0;
        uint64_t mask = UINT64_C(1) << 63;
        while ((remainder & mask) == 0) {
            ++zeros;
            mask >>= 1;
        }
        return static_cast<uint8_t>(zeros + 1);
#endif
    }

private:
    uint64_t seed_;
};

} // namespace Rank
} // namespace Sketch

#endif // RABBITSKETCH_RANK_STREAM_H
