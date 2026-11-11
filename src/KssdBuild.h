#ifndef RABBITSKETCH_KSSD_BUILD_H
#define RABBITSKETCH_KSSD_BUILD_H

#include "Sketch.h"
#include "rank/CanonicalKmer.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#if defined(__AVX2__) || defined(__AVX512BW__)
#include <immintrin.h>
#endif

namespace Sketch {
namespace KssdBuild {

// Selection happens once per span, never once per k-mer. All three paths use
// the same permutation and produce exactly the same reduced tuples.
enum class Lookup { Automatic, Bitmap, Map, Direct };

class Stream {
public:
    explicit Stream(const kssd_parameter_t& p, Lookup lookup = Lookup::Automatic,
                    int drlevel = -1)
        : p_(p), dr_(drlevel < 0 ? p.drlevel : drlevel),
          dim_end_(1u << (4 * (p.half_subk - dr_))),
          lookup_(lookup == Lookup::Automatic
              ? (p.admitted_dimensions.empty() ? Lookup::Direct : Lookup::Bitmap)
              : lookup) {}

    void reset() noexcept { forward_ = reverse_ = 0; valid_ = 0; }

    // Continue a record across spans. Reset at record boundaries. FASTA line
    // breaks are ignored, whereas every other non-ACGT byte breaks the window.
    template<bool FastaLines = false, typename Emit>
    void append(const char* sequence, size_t length, Emit&& emit) {
        switch (lookup_) {
        case Lookup::Map: appendImpl<Lookup::Map, FastaLines>(sequence, length, emit); break;
        case Lookup::Direct: appendImpl<Lookup::Direct, FastaLines>(sequence, length, emit); break;
        default: appendImpl<Lookup::Bitmap, FastaLines>(sequence, length, emit); break;
        }
    }

private:
    template<Lookup Method, bool FastaLines, typename Emit>
    void appendImpl(const char* sequence, size_t length, Emit& emit) {
        alignas(64) uint8_t encoded[4096];
        const int dimension_shift = p_.half_outctx_len * 2;
        const int context_shift = p_.kmer_size * 2 - p_.half_outctx_len * 4;
        for (size_t start = 0; start < length; start += sizeof(encoded)) {
            const size_t count = std::min(sizeof(encoded), length - start);
            const char* source = sequence + start;
            size_t pos = 0;
#if defined(__AVX512BW__)
            const __m512i lut = _mm512_broadcast_i32x4(_mm_setr_epi8(
                -1,0,-1,1,3,-1,-1,2,-1,-1,-1,-1,-1,-1,-1,-1));
            for (; pos + 64 <= count; pos += 64) {
                const __m512i raw = _mm512_loadu_si512(source + pos);
                const __m512i upper = _mm512_and_si512(raw, _mm512_set1_epi8(0xDF));
                const __mmask64 valid =
                    _mm512_cmpeq_epi8_mask(upper, _mm512_set1_epi8('A')) |
                    _mm512_cmpeq_epi8_mask(upper, _mm512_set1_epi8('C')) |
                    _mm512_cmpeq_epi8_mask(upper, _mm512_set1_epi8('G')) |
                    _mm512_cmpeq_epi8_mask(upper, _mm512_set1_epi8('T'));
                const __m512i codes = _mm512_shuffle_epi8(lut,
                    _mm512_and_si512(raw, _mm512_set1_epi8(15)));
                _mm512_storeu_si512(encoded + pos,
                    _mm512_mask_mov_epi8(_mm512_set1_epi8(-1), valid, codes));
            }
#elif defined(__AVX2__)
            const __m256i lut = _mm256_broadcastsi128_si256(_mm_setr_epi8(
                -1,0,-1,1,3,-1,-1,2,-1,-1,-1,-1,-1,-1,-1,-1));
            for (; pos + 32 <= count; pos += 32) {
                const __m256i raw = _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(source + pos));
                const __m256i upper = _mm256_and_si256(raw, _mm256_set1_epi8(0xDF));
                const __m256i valid = _mm256_or_si256(
                    _mm256_or_si256(_mm256_cmpeq_epi8(upper, _mm256_set1_epi8('A')),
                                    _mm256_cmpeq_epi8(upper, _mm256_set1_epi8('C'))),
                    _mm256_or_si256(_mm256_cmpeq_epi8(upper, _mm256_set1_epi8('G')),
                                    _mm256_cmpeq_epi8(upper, _mm256_set1_epi8('T'))));
                const __m256i codes = _mm256_shuffle_epi8(lut,
                    _mm256_and_si256(raw, _mm256_set1_epi8(15)));
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(encoded + pos),
                    _mm256_blendv_epi8(_mm256_set1_epi8(-1), codes, valid));
            }
#endif
            for (; pos < count; ++pos)
                encoded[pos] = Rank::encodeBase(static_cast<unsigned char>(source[pos]));
            for (size_t i = 0; i < count; ++i) {
                const uint8_t base = encoded[i];
                if (base > 3) {
                    if constexpr (FastaLines) {
                        if (source[i] == '\n' || source[i] == '\r') continue;
                    }
                    reset();
                    continue;
                }
                forward_ = ((forward_ << 2) | base) & p_.tupmask;
                reverse_ = (reverse_ >> 2) |
                    (static_cast<uint64_t>(base ^ 3) << p_.rev_add_move);
                if (valid_ < p_.kmer_size && ++valid_ < p_.kmer_size) continue;
                const uint64_t canonical = std::min(forward_, reverse_);
                const uint32_t dimension = (canonical & p_.domask) >> dimension_shift;
                int rank;
                if constexpr (Method == Lookup::Direct) {
                    rank = p_.shuffled_dim[dimension];
                } else {
                    if constexpr (Method == Lookup::Bitmap) {
                        if (!p_.admitted_dimensions.empty() &&
                            !(p_.admitted_dimensions[dimension >> 6] &
                              (UINT64_C(1) << (dimension & 63)))) continue;
                    }
                    const auto found = p_.shuffled_map.find(dimension);
                    if (found == p_.shuffled_map.end()) continue;
                    rank = found->second;
                }
                if (rank < p_.dim_start || static_cast<uint32_t>(rank) >= dim_end_) continue;
                const uint64_t reduced =
                    (((canonical & p_.undomask0) |
                      ((canonical & p_.undomask1) << context_shift)) >> (dr_ * 4)) |
                    static_cast<uint64_t>(rank - p_.dim_start);
                emit(reduced);
            }
        }
    }

    const kssd_parameter_t& p_;
    int dr_;
    uint32_t dim_end_;
    Lookup lookup_;
    uint64_t forward_ = 0, reverse_ = 0;
    unsigned valid_ = 0;
};

// Small sampled sets only need sequential appends followed by sort/unique.
// Cap the duplicate buffer: long or repetitive inputs switch to a hash set,
// whose memory grows with distinct tuples, not read coverage.
template<typename Key>
class Accumulator {
public:
    explicit Accumulator(bool buffered = true) : buffered_(buffered) {}
    void reset(size_t estimated = 1024) {
        values_.clear(); set_.clear(); hashed_ = !buffered_;
        if (hashed_) set_.reserve(std::min(estimated, BUFFER_LIMIT));
        else values_.reserve(std::min(estimated, BUFFER_LIMIT));
    }
    void insert(uint64_t value) {
        const Key key = static_cast<Key>(value);
        if (hashed_) { set_.insert(key); return; }
        values_.push_back(key);
        if (values_.size() == BUFFER_LIMIT) {
            set_.reserve(BUFFER_LIMIT);
            set_.insert(values_.begin(), values_.end());
            values_.clear(); hashed_ = true;
        }
    }
    void finish(std::vector<Key>& output) {
        if (hashed_) output.assign(set_.begin(), set_.end());
        else output.assign(values_.begin(), values_.end());
        std::sort(output.begin(), output.end());
        output.erase(std::unique(output.begin(), output.end()), output.end());
    }
    bool usingHashSet() const noexcept { return hashed_; }
private:
    static constexpr size_t BUFFER_LIMIT = 16384;
    bool buffered_, hashed_ = false;
    std::vector<Key> values_;
    phmap::flat_hash_set<Key> set_;
};

} // namespace KssdBuild
} // namespace Sketch
#endif
