/**
 * BinDash – b-bit One Permutation Hashing implementation.
 *
 * Algorithm:
 *   1. update(): canonical k-mer sliding window → murmur3_fmix hash →
 *      partition into nbins bins, keep minimum hash per bin.
 *   2. densify(): fill empty bins using optimal densification
 *      (Shrivastava, ICML 2017).
 *   3. packBits(): retain only the lowest bbits of each bin minimum,
 *      bit-pack into sketchsize64 * bbits uint64_t words.
 *   4. jaccard(): XNOR + AND across bbits layers, popcount matching bits.
 */

#include "BinDash.h"
#include "hash_int.h"

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstring>
#include <immintrin.h>
#include <limits>
#include <stdexcept>

using namespace Sketch;

// ═══════════════════════════════════════════════════════════════════════════
// DNA encoding (same LUT as ProbMinHash / FastKMV)
// ═══════════════════════════════════════════════════════════════════════════

static const uint8_t BD_ENC_LUT[256] = {
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
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
};

#define BD_ENC(c)    (BD_ENC_LUT[(uint8_t)(c)])
#define BD_COMP(e)   ((uint8_t)((e) <= 3 ? 3 - (e) : 255))
#define BD_VALID(e)  ((e) <= 3)

// ═══════════════════════════════════════════════════════════════════════════
// Constructor
// ═══════════════════════════════════════════════════════════════════════════

BinDash::BinDash(uint32_t sketchsize64, int kmer_size,
                 uint32_t bbits, uint64_t seed)
    : sketchsize64_(sketchsize64),
      kmer_size_(kmer_size),
      bbits_(bbits),
      seed_(seed),
      nbins_(sketchsize64 * 64),
      signs_(sketchsize64 * 64, UINT64_MAX),
      usigs_(),
      finalized_(false),
      raw_nonempty_(0)
{
    if (sketchsize64 == 0)
        throw std::invalid_argument("BinDash sketchsize64 must be positive");
    if (kmer_size < 1 || kmer_size > 32)
        throw std::invalid_argument("BinDash k-mer size must be in [1, 32]");
    if (bbits < 1 || bbits > 64)
        throw std::invalid_argument("BinDash b-bit width must be in [1, 64]");
}

BinDash BinDash::fromPacked(uint32_t sketchsize64,
                            int kmer_size,
                            uint32_t bbits,
                            uint64_t seed,
                            uint32_t raw_nonempty,
                            const std::vector<uint64_t>& signatures) {
    BinDash result(sketchsize64, kmer_size, bbits, seed);
    const size_t expected = static_cast<size_t>(sketchsize64) * bbits;
    if (signatures.size() != expected)
        throw std::invalid_argument(
            "BinDash packed signature word count is inconsistent");
    if (raw_nonempty > result.nbins_)
        throw std::invalid_argument(
            "BinDash raw nonempty count exceeds the number of bins");
    result.usigs_ = signatures;
    result.raw_nonempty_ = raw_nonempty;
    std::vector<uint64_t>().swap(result.signs_);
    result.finalized_ = true;
    return result;
}

Rank::RankMetadata BinDash::metadata() const {
    Rank::RankMetadata meta;
    meta.backend = Rank::Backend::BinDash;
    meta.hash_profile = Rank::HashProfile::BinDashDoubleFmixV1;
    meta.weight_semantics = Rank::WeightSemantics::UnweightedSet;
    meta.resolution_kind = Rank::ResolutionKind::RegisterCount;
    meta.fingerprint_bits = static_cast<uint8_t>(bbits_);
    meta.kmer_size = static_cast<uint16_t>(kmer_size_);
    meta.resolution = nbins_;
    meta.seed = seed_;
    meta.parameter_fingerprint = Rank::RankStream::fmix64(
        bbits_, UINT64_C(0x62696e646173682d));
    return meta;
}

// ═══════════════════════════════════════════════════════════════════════════
// update – canonical k-mer sliding window + bin-minimum tracking
// ═══════════════════════════════════════════════════════════════════════════

void BinDash::update(const char* seq) {
    update(seq, static_cast<uint64_t>(std::strlen(seq)));
}

void BinDash::update(const char* seq, uint64_t length) {
    if (finalized_)
        throw std::logic_error("BinDash::update called after finalize");
    if (seq == nullptr && length != 0)
        throw std::invalid_argument("BinDash::update sequence must not be null");

    const int K = kmer_size_;
    if (length < static_cast<uint64_t>(K)) return;

    const uint64_t kmer_mask = (K == 32) ? ~0ULL : ((1ULL << (2 * K)) - 1);
    const uint64_t loc_seed  = seed_;
    const uint32_t loc_nbins = nbins_;

    uint64_t fwd = 0, rev = 0;
    int inv = 0;

    for (int k = 0; k < K; ++k) {
        uint8_t ef = BD_ENC(seq[k]);
        if (!BD_VALID(ef)) inv++;
        fwd = (fwd << 2) | (BD_VALID(ef) ? (ef & 3u) : 0u);
        uint8_t er = BD_VALID(ef) ? BD_COMP(ef) : 0u;
        rev = (rev >> 2) | (static_cast<uint64_t>(er & 3u) << (2 * (K - 1)));
    }

    const uint64_t N_body = length - static_cast<uint64_t>(K);

    for (uint64_t i = 0; i <= N_body; ++i) {
        if (inv == 0) {
            const uint64_t canonical = (fwd <= rev) ? fwd : rev;
            uint64_t h = mc::murmur3_fmix(canonical, loc_seed);
            uint32_t binidx = static_cast<uint32_t>((h >> 32) % loc_nbins);
            uint64_t value  = mc::murmur3_fmix(h, 0ULL);
            if (value < signs_[binidx]) {
                if (signs_[binidx] == UINT64_MAX) ++raw_nonempty_;
                signs_[binidx] = value;
            }
        }

        if (i < N_body) {
            uint8_t ef_out = BD_ENC(seq[i]);
            uint8_t ef_in  = BD_ENC(seq[i + K]);
            if (!BD_VALID(ef_out)) inv--;
            if (!BD_VALID(ef_in))  inv++;
            fwd = ((fwd << 2) | (BD_VALID(ef_in) ? (ef_in & 3u) : 0u)) & kmer_mask;
            uint8_t er_in = BD_VALID(ef_in) ? (BD_COMP(ef_in) & 3u) : 0u;
            rev = (rev >> 2) | (static_cast<uint64_t>(er_in) << (2 * (K - 1)));
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// densify – optimal densification (Shrivastava, ICML 2017)
// ═══════════════════════════════════════════════════════════════════════════

static inline uint64_t univhash2(uint64_t s, uint64_t t) {
    uint64_t x = 1009ULL * s + 1000003ULL * t;
    return (48271ULL * x + 11ULL) % ((1ULL << 31) - 1);
}

void BinDash::densify() {
    uint64_t minval = UINT64_MAX;
    uint64_t maxval = 0;
    for (auto sign : signs_) {
        if (sign < minval) minval = sign;
        if (sign > maxval) maxval = sign;
    }
    if (maxval != UINT64_MAX) return;
    if (minval == UINT64_MAX) return;

    const uint64_t n = signs_.size();
    for (uint64_t i = 0; i < n; i++) {
        uint64_t j = i;
        uint64_t nattempts = 0;
        while (signs_[j] == UINT64_MAX) {
            j = univhash2(i, nattempts) % n;
            nattempts++;
        }
        signs_[i] = signs_[j];
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// packBits – b-bit truncation + column-major bit packing
// ═══════════════════════════════════════════════════════════════════════════

void BinDash::packBits() {
    if (usigs_.size() != static_cast<size_t>(sketchsize64_) * bbits_)
        usigs_.assign(static_cast<size_t>(sketchsize64_) * bbits_, 0ULL);
    else
        std::fill(usigs_.begin(), usigs_.end(), 0ULL);
    const uint32_t bb = bbits_;
    for (uint32_t signidx = 0; signidx < nbins_; signidx++) {
        uint64_t sign = signs_[signidx];
        uint32_t leftshift = signidx % 64;
        for (uint32_t b = 0; b < bb; b++) {
            uint64_t bit = (sign >> b) & 1ULL;
            usigs_[signidx / 64 * bb + b] |= (bit << leftshift);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// ensureFinalized – lazy densify + pack + free temporary signs_
// ═══════════════════════════════════════════════════════════════════════════

void BinDash::finalize() {
    if (finalized_) return;
    densify();
    packBits();
    { std::vector<uint64_t>().swap(signs_); }
    finalized_ = true;
}

void BinDash::ensureFinalized() const {
    if (finalized_) return;
    const_cast<BinDash*>(this)->finalize();
}

static uint64_t bd_count_same_packed(const uint64_t* a, const uint64_t* b,
                                     uint32_t sketchsize64, uint32_t bbits) {
    uint64_t same = 0;
    for (uint32_t g = 0; g < sketchsize64; ++g) {
        uint64_t eqmask = ~0ULL;
        const uint32_t base = g * bbits;
        uint32_t bit = 0;

#if defined(__AVX512F__)
        const __m512i all1_512 = _mm512_set1_epi64(-1LL);
        __m512i eqv512 = all1_512;
        for (; bit + 8 <= bbits; bit += 8) {
            const __m512i va = _mm512_loadu_si512((const void*)(a + base + bit));
            const __m512i vb = _mm512_loadu_si512((const void*)(b + base + bit));
            const __m512i x  = _mm512_xor_si512(va, vb);
            const __m512i xn = _mm512_andnot_si512(x, all1_512); // XNOR
            eqv512 = _mm512_and_si512(eqv512, xn);
        }
        alignas(64) uint64_t lanes512[8];
        _mm512_store_si512((void*)lanes512, eqv512);
        eqmask &= lanes512[0] & lanes512[1] & lanes512[2] & lanes512[3]
               & lanes512[4] & lanes512[5] & lanes512[6] & lanes512[7];
#elif defined(__AVX2__)
        const __m256i all1_256 = _mm256_set1_epi64x(-1LL);
        __m256i eqv256 = all1_256;
        for (; bit + 4 <= bbits; bit += 4) {
            const __m256i va = _mm256_loadu_si256((const __m256i*)(a + base + bit));
            const __m256i vb = _mm256_loadu_si256((const __m256i*)(b + base + bit));
            const __m256i x  = _mm256_xor_si256(va, vb);
            const __m256i xn = _mm256_andnot_si256(x, all1_256); // XNOR
            eqv256 = _mm256_and_si256(eqv256, xn);
        }
        alignas(32) uint64_t lanes256[4];
        _mm256_store_si256((__m256i*)lanes256, eqv256);
        eqmask &= lanes256[0] & lanes256[1] & lanes256[2] & lanes256[3];
#endif

        for (; bit < bbits; ++bit)
            eqmask &= ~(a[base + bit] ^ b[base + bit]); // XNOR
        same += static_cast<uint64_t>(__builtin_popcountll(eqmask));
    }
    return same;
}

// ═══════════════════════════════════════════════════════════════════════════
// Public static wrapper — enables flat-array distance computation without
// constructing BinDash objects.
// ═══════════════════════════════════════════════════════════════════════════

uint64_t BinDash::countSameBits(const uint64_t* a, const uint64_t* b,
                                uint32_t sk64, uint32_t bb) {
    return bd_count_same_packed(a, b, sk64, bb);
}

// ═══════════════════════════════════════════════════════════════════════════
// jaccard – XNOR + AND across bbits layers, popcount
// ═══════════════════════════════════════════════════════════════════════════

double BinDash::jaccardPacked(const BinDash& other) const {
    metadata().requireSameFamily(other.metadata(), "BinDash::jaccard");

    ensureFinalized();
    other.ensureFinalized();

    const uint64_t samebits = bd_count_same_packed(usigs_.data(), other.usigs_.data(),
                                                   sketchsize64_, bbits_);
    const double maxnbits = static_cast<double>(nbins_);
    const double p_match  = static_cast<double>(samebits) / maxnbits;
    const double p_random = std::ldexp(1.0, -static_cast<int>(bbits_));
    double j = (p_match - p_random) / (1.0 - p_random);
    if (j < 0.0) j = 0.0;
    if (j > 1.0) j = 1.0;
    return j;
}

double BinDash::jaccard(const BinDash& other) const {
    return jaccardPacked(other);
}

double BinDash::cardinalityEstimate() const {
    const double m = static_cast<double>(nbins_);
    const double occupied = std::min(static_cast<double>(raw_nonempty_), m);
    if (occupied <= 0.0) return 0.0;
    if (occupied >= m) return std::numeric_limits<double>::infinity();
    return -m * std::log1p(-occupied / m);
}

double BinDash::containment(const BinDash& other) const {
    const double card_a = cardinalityEstimate();
    if (card_a <= 0.0 || !std::isfinite(card_a)) return 0.0;
    const double j = jaccard(other);
    if (j <= 0.0) return 0.0;
    const double card_b = other.cardinalityEstimate();
    if (card_b <= 0.0) return 0.0;
    return j * (card_a + card_b) / (card_a * (1.0 + j));
}

double BinDash::ani(const BinDash& other) const {
    const double j = jaccard(other);
    if (j <= 0.0) return 0.0;
    if (j >= 1.0) return 1.0;
    return std::pow(2.0 * j / (1.0 + j),
                    1.0 / static_cast<double>(kmer_size_));
}

double BinDash::distance(const BinDash& other) const {
    const double j = jaccard(other);
    if (j <= 0.0) return std::numeric_limits<double>::infinity();
    if (j >= 1.0) return 0.0;
    const double ratio = 2.0 * j / (1.0 + j);
    return -std::log(ratio) / static_cast<double>(kmer_size_);
}

#undef BD_ENC
#undef BD_COMP
#undef BD_VALID
