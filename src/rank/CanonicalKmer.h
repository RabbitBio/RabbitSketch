#ifndef RABBITSKETCH_CANONICAL_KMER_H
#define RABBITSKETCH_CANONICAL_KMER_H

#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace Sketch {
namespace Rank {

enum class Canonicalization : uint16_t {
    Lexicographic2BitReverseComplementV1 = 1,
    ForwardOnlyV1 = 2
};

inline uint8_t encodeBase(unsigned char base) noexcept {
    static constexpr uint8_t TABLE[256] = {
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,  0,255,  1,255,255,255,  2,255,255,255,255,255,255,255,255,
        255,255,255,255,  3,255,255,255,255,255,255,255,255,255,255,255,
        255,  0,255,  1,255,255,255,  2,255,255,255,255,255,255,255,255,
        255,255,255,255,  3,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255
    };
    return TABLE[base];
}

inline bool validBase(uint8_t encoded) noexcept { return encoded <= 3; }
inline uint8_t complementBase(uint8_t encoded) noexcept {
    return validBase(encoded) ? static_cast<uint8_t>(3 - encoded) : 255;
}

inline uint64_t kmerMask(uint8_t kmer_size) {
    if (kmer_size == 0 || kmer_size > 32) {
        throw std::invalid_argument("k-mer size must be in [1, 32]");
    }
    return kmer_size == 32
        ? UINT64_MAX
        : ((UINT64_C(1) << (2 * kmer_size)) - 1);
}

inline uint64_t canonicalCode(uint64_t forward,
                              uint64_t reverse_complement) noexcept {
    return forward <= reverse_complement ? forward : reverse_complement;
}

/**
 * Correctness-first iterator for the shared canonical k-mer contract.
 *
 * Optimized backends may retain vectorized rolling loops, but their golden
 * tests are checked against this iterator.  Invalid-base windows are skipped.
 */
class CanonicalKmerIterator {
public:
    CanonicalKmerIterator(const char* sequence, size_t length, uint8_t kmer_size)
        : sequence_(sequence), length_(length), kmer_size_(kmer_size),
          mask_(kmerMask(kmer_size)), total_(0), position_(0),
          forward_(0), reverse_(0), invalid_(0) {
        if (sequence_ == nullptr && length_ != 0) {
            throw std::invalid_argument("sequence must not be null");
        }
        if (length_ < kmer_size_) return;
        total_ = length_ - kmer_size_ + 1;
        for (uint8_t i = 0; i < kmer_size_; ++i) {
            const uint8_t encoded = encodeBase(
                static_cast<unsigned char>(sequence_[i]));
            if (!validBase(encoded)) ++invalid_;
            const uint8_t safe = validBase(encoded) ? encoded : 0;
            forward_ = ((forward_ << 2) | safe) & mask_;
            reverse_ = (reverse_ >> 2)
                | (static_cast<uint64_t>(validBase(encoded)
                      ? complementBase(encoded) : 0)
                   << (2 * (kmer_size_ - 1)));
        }
    }

    bool next(uint64_t& canonical_code, size_t* position = nullptr) {
        while (position_ < total_) {
            const size_t current = position_;
            const bool valid = invalid_ == 0;
            if (valid) canonical_code = canonicalCode(forward_, reverse_);
            advance();
            if (valid) {
                if (position != nullptr) *position = current;
                return true;
            }
        }
        return false;
    }

private:
    void advance() noexcept {
        if (position_ + 1 < total_) {
            const uint8_t outgoing = encodeBase(
                static_cast<unsigned char>(sequence_[position_]));
            const uint8_t incoming = encodeBase(
                static_cast<unsigned char>(sequence_[position_ + kmer_size_]));
            if (!validBase(outgoing)) --invalid_;
            if (!validBase(incoming)) ++invalid_;
            const uint8_t safe = validBase(incoming) ? incoming : 0;
            forward_ = ((forward_ << 2) | safe) & mask_;
            reverse_ = (reverse_ >> 2)
                | (static_cast<uint64_t>(validBase(incoming)
                      ? complementBase(incoming) : 0)
                   << (2 * (kmer_size_ - 1)));
        }
        ++position_;
    }

    const char* sequence_;
    size_t length_;
    uint8_t kmer_size_;
    uint64_t mask_;
    size_t total_;
    size_t position_;
    uint64_t forward_;
    uint64_t reverse_;
    int invalid_;
};

} // namespace Rank
} // namespace Sketch

#endif // RABBITSKETCH_CANONICAL_KMER_H
