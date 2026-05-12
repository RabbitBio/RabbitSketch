/**
 * FastKMV – K Minimum Values sketch with ntHash rolling + fmix finalizer.
 *
 * Lazy-sort warmup: the first 2k keys are appended unsorted in O(1).
 * When the buffer fills, compactify() sorts, deduplicates, and truncates
 * to the bottom-k in one pass, producing a sorted array with exact
 * threshold — identical to the original design from that point on.
 *
 * Warmup (sorted_==false): O(1) append into 2k buffer; at 2k entries,
 * compactify() sorts, dedupes, truncates to bottom-k (may repeat if still <k).
 * Steady state: sorted bottom-k, O(log k) lower_bound + O(k) memmove;
 * threshold_ == v[k-1] is exact.
 *
 * Rolling hash: ntHash (Mohamadi et al. 2016), canonical = min(fwd, rc).
 * Finalizer: murmur3 fmix (default 1 round; compile with -DFASTKMV_DOUBLE_FMUX
 * for 2 rounds), 8-wide AVX-512 when available.
 *
 * Ablation: -DFASTKMV_NO_FMUX skips the finalizer and uses (canonical_ntHash >> 11)
 * as the sketch key (testing only; expect different Jaccard accuracy).
 *
 * Jaccard estimator (standard KMV):
 *   merge two sorted sketches, take k smallest distinct values, count
 *   how many appear in both:  J ≈ common / distinct.
 */

#include "fastkmv.h"
#include "hash_int.h"
#include <cmath>

// Default: single fmix (faster). Define FASTKMV_DOUBLE_FMUX for two rounds (legacy).
// Define FASTKMV_NO_FMUX to disable fmix entirely (raw ntHash keys; ablation test).
#ifdef FASTKMV_DOUBLE_FMUX
#define FASTKMV_FMUX_ROUNDS 2
#else
#define FASTKMV_FMUX_ROUNDS 1
#endif

#include <immintrin.h>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cassert>
#include <climits>
#include <limits>
#include <memory>

using namespace Sketch;

// ═══════════════════════════════════════════════════════════════════════════
// DNA encoding via bit operations (no 256-byte LUT)
// ═══════════════════════════════════════════════════════════════════════════

static inline uint8_t fkmv_enc_or_invalid(uint8_t c) {
    const uint8_t cu = c & 0xDF;
    if (cu == 'A' || cu == 'C' || cu == 'G' || cu == 'T')
        return (c >> 1) & 3;
    return 0xFF;
}
#define FKMV_VALID(e) ((e) <= 3)

// ═══════════════════════════════════════════════════════════════════════════
// ntHash seed tables  (Mohamadi et al., Bioinformatics 2016)
// ═══════════════════════════════════════════════════════════════════════════

static const uint64_t NT_SEED_FWD[4] = {
    0x3c8bfbb395c60474ULL,
    0x3193c18562a02b4cULL,
    0x295549f54be24456ULL,
    0x20323ed082572324ULL,
};
static const uint64_t NT_SEED_RC[4] = {
    NT_SEED_FWD[2], NT_SEED_FWD[3], NT_SEED_FWD[0], NT_SEED_FWD[1],
};

static inline uint64_t rol64(uint64_t x, unsigned n) {
    n &= 63;
    return (x << n) | (x >> ((64 - n) & 63));
}

#ifdef __AVX2__
static inline __m256i avx2_mullo_epi64(__m256i a, __m256i b) {
    __m256i a_hi  = _mm256_srli_epi64(a, 32);
    __m256i b_hi  = _mm256_srli_epi64(b, 32);
    __m256i lo_lo = _mm256_mul_epu32(a, b);
    __m256i cross = _mm256_add_epi64(_mm256_mul_epu32(a_hi, b),
                                      _mm256_mul_epu32(a, b_hi));
    return _mm256_add_epi64(lo_lo, _mm256_slli_epi64(cross, 32));
}
#endif

// ═══════════════════════════════════════════════════════════════════════════
// FastKMV  –  constructor / copy / assign
// ═══════════════════════════════════════════════════════════════════════════

FastKMV::FastKMV(uint32_t k, int kmer_size, uint64_t seed)
    : k_(k), kmer_size_(kmer_size), seed_(seed),
      buf_cap_(k * 2),
      vals_(new uint64_t[k * 2]),
      enc_buf_(nullptr),
      enc_cap_(0),
      size_(0),
      threshold_(UINT64_MAX),
      sorted_(false)
{
    assert(k > 1);
    assert(kmer_size >= 1 && kmer_size <= 32);
    // Warmup: append-only until buf full → compactify (see insertKey).
    // Sentinels UINT64_MAX are applied in compactify / merge, not here.
}

FastKMV::FastKMV(const FastKMV& o)
    : k_(o.k_), kmer_size_(o.kmer_size_), seed_(o.seed_),
      buf_cap_(o.buf_cap_),
      vals_(new uint64_t[o.buf_cap_]),
      enc_buf_(nullptr),
      enc_cap_(0),
      size_(o.size_),
      threshold_(o.threshold_),
      sorted_(o.sorted_)
{
    std::copy(o.vals_.get(), o.vals_.get() + size_, vals_.get());
    if (sorted_ && size_ < k_)
        std::fill(vals_.get() + size_, vals_.get() + k_, UINT64_MAX);
}

FastKMV& FastKMV::operator=(FastKMV other) {
    std::swap(k_,         other.k_);
    std::swap(kmer_size_, other.kmer_size_);
    std::swap(seed_,      other.seed_);
    std::swap(buf_cap_,   other.buf_cap_);
    std::swap(vals_,      other.vals_);
    std::swap(enc_buf_,   other.enc_buf_);
    std::swap(enc_cap_,   other.enc_cap_);
    std::swap(size_,      other.size_);
    std::swap(threshold_, other.threshold_);
    std::swap(sorted_,    other.sorted_);
    return *this;
}

// ═══════════════════════════════════════════════════════════════════════════
// insertKey
//
//   Phase A — warmup (!sorted_): O(1) append; threshold_ is UINT64_MAX.
//   Phase B — sorted but unfilled (sorted_ && size_<k): lower_bound + memmove.
//   Phase C — full KMV (sorted_ && size_==k): threshold + lower_bound + memmove.
// ═══════════════════════════════════════════════════════════════════════════

void FastKMV::insertKey(uint64_t key) {
    if (key >= threshold_) return;

    if (!sorted_) {
        vals_[size_++] = key;
        if (size_ >= buf_cap_)
            compactify();
        return;
    }

    uint64_t* v = vals_.get();

    if (size_ < k_) {
        uint64_t* pos = std::lower_bound(v, v + size_, key);
        uint32_t idx = static_cast<uint32_t>(pos - v);
        if (idx < size_ && v[idx] == key) return;
        std::memmove(v + idx + 1, v + idx, (size_ - idx) * sizeof(uint64_t));
        v[idx] = key;
        ++size_;
        if (size_ == k_)
            threshold_ = v[k_ - 1];
        return;
    }

    uint64_t* pos = std::lower_bound(v, v + k_, key);
    uint32_t idx = static_cast<uint32_t>(pos - v);
    if (idx < k_ && v[idx] == key) return;
    std::memmove(v + idx + 1, v + idx, (k_ - 1 - idx) * sizeof(uint64_t));
    v[idx] = key;
    threshold_ = v[k_ - 1];
}

// ═══════════════════════════════════════════════════════════════════════════
// compactify  –  sort, deduplicate, truncate to k
// ═══════════════════════════════════════════════════════════════════════════

void FastKMV::compactify() const {
    uint64_t* v = vals_.get();
    std::sort(v, v + size_);

    uint32_t w = 0;
    for (uint32_t r = 0; r < size_; ++r) {
        if (r == 0 || v[r] != v[r - 1]) {
            v[w++] = v[r];
            if (w == k_) break;
        }
    }
    size_ = w;
    threshold_ = (size_ >= k_) ? v[k_ - 1] : UINT64_MAX;
    std::fill(v + size_, v + k_, UINT64_MAX);
    sorted_ = true;
}

void FastKMV::ensureSorted() const {
    if (!sorted_)
        compactify();
}

// ═══════════════════════════════════════════════════════════════════════════
// addHash  –  hash → 53-bit integer key → insert into bottom-k
// ═══════════════════════════════════════════════════════════════════════════

void FastKMV::addHash(uint64_t h) {
#ifdef FASTKMV_NO_FMUX
    insertKey(h >> 11);
#else
    uint64_t h1 = mc::murmur3_fmix(h, seed_);
    insertKey(h1 >> 11);
#endif
}

// ═══════════════════════════════════════════════════════════════════════════
// update  –  ntHash rolling + SIMD fmix finalizer + threshold filter
// ═══════════════════════════════════════════════════════════════════════════

void FastKMV::finalize() {
    enc_buf_.reset();
    enc_cap_ = 0;
}

void FastKMV::update(const char* seq, uint64_t length) {
    const int K = kmer_size_;
    if (length < static_cast<uint64_t>(K)) return;

    const uint64_t loc_seed = seed_;

    // ── Phase 0: SIMD bulk sequence encoding ────────────────────────────
    if (enc_cap_ < length) {
        enc_buf_.reset(new uint8_t[length]);
        enc_cap_ = length;
    }
    uint8_t* enc = enc_buf_.get();

    uint64_t p = 0;
#if defined(__AVX512BW__)
    {
        const __m512i vlut = _mm512_broadcast_i32x4(_mm_setr_epi8(
            -1,  0, -1,  1,  2, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1));
        const __m512i vmask_lo = _mm512_set1_epi8(0x0F);
        for (; p + 64 <= length; p += 64) {
            __m512i vraw = _mm512_loadu_si512(seq + p);
            __m512i vnib = _mm512_and_si512(vraw, vmask_lo);
            __m512i venc = _mm512_shuffle_epi8(vlut, vnib);
            _mm512_storeu_si512(enc + p, venc);
        }
    }
#elif defined(__AVX2__)
    {
        const __m256i vlut = _mm256_broadcastsi128_si256(_mm_setr_epi8(
            -1,  0, -1,  1,  2, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1));
        const __m256i vmask_lo = _mm256_set1_epi8(0x0F);
        for (; p + 32 <= length; p += 32) {
            __m256i vraw = _mm256_loadu_si256((const __m256i*)(seq + p));
            __m256i vnib = _mm256_and_si256(vraw, vmask_lo);
            __m256i venc = _mm256_shuffle_epi8(vlut, vnib);
            _mm256_storeu_si256((__m256i*)(enc + p), venc);
        }
    }
#endif
    for (; p < length; ++p)
        enc[p] = fkmv_enc_or_invalid((uint8_t)seq[p]);

    // ── Phase 1: ntHash rolling hash ────────────────────────────────────
    uint64_t sf_k[4], sc_ror1[4], sc_km1[4];
    const unsigned km1 = static_cast<unsigned>(K > 1 ? K - 1 : 63);
    for (int i = 0; i < 4; ++i) {
        sf_k[i]    = rol64(NT_SEED_FWD[i], static_cast<unsigned>(K));
        sc_ror1[i] = rol64(NT_SEED_RC[i],  63u);
        sc_km1[i]  = rol64(NT_SEED_RC[i],  km1);
    }

    uint64_t fwd_h = 0, rc_h = 0;
    int inv = 0;
    for (int k = 0; k < K; ++k) {
        uint8_t e = enc[k];
        if (e > 3) { ++inv; e = 0; }
        fwd_h ^= rol64(NT_SEED_FWD[e], static_cast<unsigned>(K - 1 - k));
        rc_h  ^= rol64(NT_SEED_RC[e],  static_cast<unsigned>(k));
    }

    const int      lanes   = 8;
    const uint64_t N_body  = length - static_cast<uint64_t>(K);
    const uint64_t N_batch = (N_body >= 1) ? (N_body / lanes) * lanes : 0;

#if defined(__AVX512F__) && defined(__AVX512DQ__) && !defined(FASTKMV_NO_FMUX)
    const __m512i vs  = _mm512_set1_epi64((int64_t)loc_seed);
    const __m512i vc1 = _mm512_set1_epi64(0xff51afd7ed558ccdLL);
    const __m512i vc2 = _mm512_set1_epi64(0xc4ceb9fe1a85ec53LL);
#endif
#if defined(__AVX512F__) && defined(__AVX512DQ__)
    __m512i vthresh = _mm512_set1_epi64((int64_t)threshold_);
#elif defined(__AVX2__)
#ifndef FASTKMV_NO_FMUX
    const __m256i vs  = _mm256_set1_epi64x((long long)loc_seed);
    const __m256i vc1 = _mm256_set1_epi64x(0xff51afd7ed558ccdLL);
    const __m256i vc2 = _mm256_set1_epi64x(0xc4ceb9fe1a85ec53LL);
#endif
    const __m256i vsign = _mm256_set1_epi64x((long long)0x8000000000000000ULL);
    __m256i vthresh = _mm256_set1_epi64x((long long)threshold_);
#endif

    for (uint64_t i = 0; i < N_batch; i += lanes) {
        uint64_t resv[8];
        bool     lane_valid[8];

        for (int j = 0; j < lanes; ++j) {
            const uint64_t pos = i + j;

            lane_valid[j] = (inv == 0);
            resv[j] = lane_valid[j] ? (fwd_h < rc_h ? fwd_h : rc_h) : 0;

            const uint8_t e_out = enc[pos];
            const uint8_t e_in  = enc[pos + K];
            if (e_out > 3) --inv;
            if (e_in  > 3) ++inv;
            const uint8_t oi = (e_out <= 3) ? e_out : 0;
            const uint8_t ii = (e_in  <= 3) ? e_in  : 0;

            fwd_h = rol64(fwd_h, 1)  ^ sf_k[oi]   ^ NT_SEED_FWD[ii];
            rc_h  = rol64(rc_h, 63u) ^ sc_ror1[oi] ^ sc_km1[ii];
        }

#if defined(__AVX512F__) && defined(__AVX512DQ__)
        __mmask8 valid_mask = 0;
        for (int j = 0; j < lanes; ++j)
            if (lane_valid[j]) valid_mask |= (1u << j);
        if (!valid_mask) continue;

#ifdef FASTKMV_NO_FMUX
        __m512i vkeys = _mm512_srli_epi64(_mm512_loadu_si512((const void*)resv), 11);
#else
        __m512i vb = _mm512_loadu_si512((const void*)resv);
        __m512i va = _mm512_xor_epi64(vb, vs);
        __m512i vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);
        va = _mm512_mullo_epi64(vb, vc1);
        vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);
        va = _mm512_mullo_epi64(vb, vc2);
        vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);

#if FASTKMV_FMUX_ROUNDS >= 2
        va = _mm512_xor_epi64(vb, vs);
        vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);
        va = _mm512_mullo_epi64(vb, vc1);
        vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);
        va = _mm512_mullo_epi64(vb, vc2);
        vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);
#endif

        __m512i vkeys = _mm512_srli_epi64(vb, 11);
#endif
        __mmask8 pass_mask = _mm512_mask_cmp_epu64_mask(
            valid_mask, vkeys, vthresh, _MM_CMPINT_LT);

        if (pass_mask) {
            uint64_t keyv[8];
            _mm512_storeu_si512(keyv, vkeys);
            while (pass_mask) {
                int j = __builtin_ctz(pass_mask);
                pass_mask &= pass_mask - 1;
                insertKey(keyv[j]);
            }
            vthresh = _mm512_set1_epi64((int64_t)threshold_);
        }
#elif defined(__AVX2__)
        for (int half = 0; half < 2; ++half) {
            const int base = half * 4;
            if (!(lane_valid[base]|lane_valid[base+1]|lane_valid[base+2]|lane_valid[base+3]))
                continue;

#ifdef FASTKMV_NO_FMUX
            __m256i vkeys = _mm256_srli_epi64(
                _mm256_loadu_si256((const __m256i*)(resv + base)), 11);
#else
            __m256i vb = _mm256_loadu_si256((const __m256i*)(resv + base));
            __m256i va = _mm256_xor_si256(vb, vs);
            __m256i vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);
            va = avx2_mullo_epi64(vb, vc1);
            vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);
            va = avx2_mullo_epi64(vb, vc2);
            vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);

#ifndef PMH_FAST_HASH
            va = _mm256_xor_si256(vb, vs);
            vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);
            va = avx2_mullo_epi64(vb, vc1);
            vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);
            va = avx2_mullo_epi64(vb, vc2);
            vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);
#endif

            __m256i vkeys = _mm256_srli_epi64(vb, 11);
#endif
            __m256i cmp = _mm256_cmpgt_epi64(
                _mm256_xor_si256(vthresh, vsign),
                _mm256_xor_si256(vkeys,   vsign));
            int pass_bits = _mm256_movemask_pd(_mm256_castsi256_pd(cmp));

            if (pass_bits) {
                uint64_t keyv[4];
                _mm256_storeu_si256((__m256i*)keyv, vkeys);
                for (int j = 0; j < 4; ++j) {
                    if ((pass_bits & (1 << j)) && lane_valid[base + j])
                        insertKey(keyv[j]);
                }
                vthresh = _mm256_set1_epi64x((long long)threshold_);
            }
        }
#else
        if (!(lane_valid[0]|lane_valid[1]|lane_valid[2]|lane_valid[3]|
              lane_valid[4]|lane_valid[5]|lane_valid[6]|lane_valid[7])) continue;

        for (int j = 0; j < lanes; ++j) {
            if (!lane_valid[j]) continue;
#ifdef FASTKMV_NO_FMUX
            uint64_t key = resv[j] >> 11;
#else
            uint64_t h0 = mc::murmur3_fmix(resv[j], loc_seed);
#if FASTKMV_FMUX_ROUNDS >= 2
            uint64_t h1 = mc::murmur3_fmix(h0, loc_seed);
#else
            uint64_t h1 = h0;
#endif
            uint64_t key = h1 >> 11;
#endif
            if (key >= threshold_) continue;
            insertKey(key);
        }
#endif
    }

    // ── Remainder loop ──────────────────────────────────────────────────
    for (uint64_t i = N_batch; i <= N_body; ++i) {
        if (inv == 0) {
            uint64_t canon = (fwd_h < rc_h) ? fwd_h : rc_h;
#ifdef FASTKMV_NO_FMUX
            uint64_t key = canon >> 11;
#else
            uint64_t h0 = mc::murmur3_fmix(canon, loc_seed);
#if FASTKMV_FMUX_ROUNDS >= 2
            uint64_t h1 = mc::murmur3_fmix(h0, loc_seed);
#else
            uint64_t h1 = h0;
#endif
            uint64_t key = h1 >> 11;
#endif
            if (key < threshold_)
                insertKey(key);
        }
        if (i < N_body) {
            const uint8_t e_out = enc[i];
            const uint8_t e_in  = enc[i + K];
            if (e_out > 3) --inv;
            if (e_in  > 3) ++inv;
            const uint8_t oi = (e_out <= 3) ? e_out : 0;
            const uint8_t ii = (e_in  <= 3) ? e_in  : 0;
            fwd_h = rol64(fwd_h, 1)  ^ sf_k[oi]   ^ NT_SEED_FWD[ii];
            rc_h  = rol64(rc_h, 63u) ^ sc_ror1[oi] ^ sc_km1[ii];
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// jaccard  –  SIMD-accelerated KMV bottom-k intersection estimator
//
// Sorted-set intersection via rotational comparison:
//   AVX-512: 8-wide (8×8 all-pairs per iteration)
//   AVX2:    4-wide (4×4 all-pairs per iteration)
//   Scalar fallback for remainder / small sketches.
// ═══════════════════════════════════════════════════════════════════════════

static uint64_t fkmv_scalar_intersect(const uint64_t* list1, uint64_t size1,
                                      const uint64_t* list2, uint64_t size2,
                                      uint64_t budget,
                                      uint64_t* i_a, uint64_t* i_b)
{
    uint64_t counter = 0;
    *i_a = 0;
    *i_b = 0;
    while (*i_a < size1 && *i_b < size2 && budget > 0) {
        if (list1[*i_a] < list2[*i_b])      { ++(*i_a); }
        else if (list1[*i_a] > list2[*i_b]) { ++(*i_b); }
        else { ++counter; ++(*i_a); ++(*i_b); }
        --budget;
    }
    return counter;
}

double FastKMV::jaccard(const FastKMV& other, double min_jaccard) const {
    assert(k_ == other.k_);
    ensureSorted();
    other.ensureSorted();

    const uint64_t* __restrict__ a = vals_.get();
    const uint64_t* __restrict__ b = other.vals_.get();
    const uint64_t sa = size_;
    const uint64_t sb = other.size_;
    const uint64_t k  = k_;

    if (sa == 0 || sb == 0) return 0.0;

    // Jaccard-threshold prefilter (active iff min_jaccard > 0):
    // KMV denominator d ≤ k, so to reach J ≥ min_jaccard the number of
    // matching hashes must satisfy common ≥ ceil(min_jaccard * k).  When
    // the achievable common falls below this bound we abort and return
    // 0.0, which is strictly below any positive threshold.
    const uint64_t need = (min_jaccard > 0.0)
        ? static_cast<uint64_t>(std::ceil(min_jaccard * static_cast<double>(k)))
        : 0ULL;
    if (need > 0 && std::min(sa, sb) < need) return 0.0;

    uint64_t ia = 0, ib = 0;
    uint64_t common = 0;

#if defined(__AVX512F__)
    {
        const uint64_t st_a = (sa / 8) * 8;
        const uint64_t st_b = (sb / 8) * 8;

        if (k > 8 && st_a > 0 && st_b > 0) {
            const uint64_t stop = k - 8;

            __m512i sv0 = _mm512_set_epi64(0,7,6,5,4,3,2,1);
            __m512i sv1 = _mm512_set_epi64(1,0,7,6,5,4,3,2);
            __m512i sv2 = _mm512_set_epi64(2,1,0,7,6,5,4,3);
            __m512i sv3 = _mm512_set_epi64(3,2,1,0,7,6,5,4);
            __m512i sv4 = _mm512_set_epi64(4,3,2,1,0,7,6,5);
            __m512i sv5 = _mm512_set_epi64(5,4,3,2,1,0,7,6);
            __m512i sv6 = _mm512_set_epi64(6,5,4,3,2,1,0,7);

            while (ia < st_a && ib < st_b) {
                __m512i va = _mm512_loadu_si512(&a[ia]);
                __m512i vb = _mm512_loadu_si512(&b[ib]);

                uint64_t a_max = a[ia + 7];
                uint64_t b_max = b[ib + 7];

                ia += (a_max <= b_max) * 8;
                ib += (a_max >= b_max) * 8;

                __mmask8 cmp0 = _mm512_cmpeq_epu64_mask(va, vb);
                __mmask8 cmp1 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv0, vb));
                __mmask8 cmp2 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv1, vb));
                __mmask8 cmp3 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv2, vb));
                __mmask8 cmpA = cmp0 | cmp1 | cmp2 | cmp3;

                __mmask8 cmp4 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv3, vb));
                __mmask8 cmp5 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv4, vb));
                __mmask8 cmp6 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv5, vb));
                __mmask8 cmp7 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv6, vb));
                __mmask8 cmpB = cmp4 | cmp5 | cmp6 | cmp7;

                __mmask8 hits = cmpA | cmpB;
                common += _mm_popcnt_u64(hits);

                if (need > 0 &&
                    common + std::min(sa - ia, sb - ib) < need) {
                    return 0.0;
                }

                if (ia + ib - common >= stop) {
                    common -= _mm_popcnt_u64(hits);
                    ia -= (a_max <= b_max) * 8;
                    ib -= (a_max >= b_max) * 8;
                    break;
                }
            }
        }
    }
#elif defined(__AVX2__)
    {
        const uint64_t st_a = (sa / 4) * 4;
        const uint64_t st_b = (sb / 4) * 4;

        if (k > 8 && st_a > 0 && st_b > 0) {
            const uint64_t stop = k - 8;

            while (ia < st_a && ib < st_b) {
                __m256i va = _mm256_loadu_si256((const __m256i*)&a[ia]);
                __m256i vb = _mm256_loadu_si256((const __m256i*)&b[ib]);

                uint64_t a_max = a[ia + 3];
                uint64_t b_max = b[ib + 3];

                ia += (a_max <= b_max) * 4;
                ib += (a_max >= b_max) * 4;

                __m256i cmp1 = _mm256_cmpeq_epi64(va, vb);
                __m256i rot1 = _mm256_permute4x64_epi64(vb, 0x39);
                __m256i cmp2 = _mm256_cmpeq_epi64(va, rot1);
                __m256i rot2 = _mm256_permute4x64_epi64(vb, 0x4E);
                __m256i cmp3 = _mm256_cmpeq_epi64(va, rot2);
                __m256i rot3 = _mm256_permute4x64_epi64(vb, 0x93);
                __m256i cmp4 = _mm256_cmpeq_epi64(va, rot3);

                __m256i combined = _mm256_or_si256(
                    _mm256_or_si256(cmp1, cmp2),
                    _mm256_or_si256(cmp3, cmp4));
                int mask = _mm256_movemask_pd(_mm256_castsi256_pd(combined));
                common += _mm_popcnt_u64(static_cast<unsigned>(mask));

                if (need > 0 &&
                    common + std::min(sa - ia, sb - ib) < need) {
                    return 0.0;
                }

                if (ia + ib - common >= stop) break;
            }
        }
    }
#endif

    // scalar remainder
    uint64_t remaining = k - (ia + ib - common);
    uint64_t ia_s, ib_s;
    common += fkmv_scalar_intersect(a + ia, sa - ia, b + ib, sb - ib,
                                    remaining, &ia_s, &ib_s);
    ia += ia_s;
    ib += ib_s;

    if (need > 0 && common < need) return 0.0;

    uint64_t distinct = ia + ib - common;
    while (distinct < k && ia < sa) { ++distinct; ++ia; }
    while (distinct < k && ib < sb) { ++distinct; ++ib; }
    if (distinct > k) distinct = k;

    return (distinct > 0) ? static_cast<double>(common) / static_cast<double>(distinct)
                          : 0.0;
}

// ═══════════════════════════════════════════════════════════════════════════
// cardinality  –  KMV estimator: (valid_k - 1) * KEY_MAX / tau_k
// ═══════════════════════════════════════════════════════════════════════════

double FastKMV::cardinality() const {
    ensureSorted();
    if (size_ == 0) return 0.0;
    if (size_ == 1) return 1.0;
    // vals_[size_-1] is the largest (= k-th minimum when sketch is full).
    // For a partial sketch (size_ < k_) we know the exact count.
    if (size_ < k_) return static_cast<double>(size_);
    const uint64_t tau = vals_[k_ - 1];
    if (tau == 0) return static_cast<double>(k_);
    return static_cast<double>(k_ - 1) * static_cast<double>(KEY_MAX)
           / static_cast<double>(tau);
}

// ═══════════════════════════════════════════════════════════════════════════
// containment  –  C(this ⊆ other) = |A∩B| / |A|
//   Uses cardinality-based formula:
//     C = J * (cardA + cardB) / (cardA * (1 + J))
// ═══════════════════════════════════════════════════════════════════════════

double FastKMV::containment(const FastKMV& other) const {
    const double card_a = cardinality();
    if (card_a <= 0.0) return 0.0;
    const double j = jaccard(other);
    if (j <= 0.0) return 0.0;
    const double card_b = other.cardinality();
    // |AUB| = (cardA + cardB) / (1 + J)  =>  |A∩B| = J * |AUB|
    return j * (card_a + card_b) / (card_a * (1.0 + j));
}

// ═══════════════════════════════════════════════════════════════════════════
// distance  –  Mash distance from Jaccard
//   D = -ln(2J/(1+J)) / k   (Ondov et al. 2016)
// ═══════════════════════════════════════════════════════════════════════════

double FastKMV::distance(const FastKMV& other, double max_distance) const {
    // Derive the minimum Jaccard that still meets max_distance via the
    // inverse Mash formula: D = -ln(2J/(1+J))/k  ⟹  J = 1/(2·exp(k·D) - 1).
    // When max_distance = +∞ (default) we get min_jaccard = 0 ⟹ no pruning.
    double min_jaccard = 0.0;
    if (std::isfinite(max_distance) && max_distance >= 0.0) {
        const double kd    = static_cast<double>(kmer_size_) * max_distance;
        const double denom = 2.0 * std::exp(kd) - 1.0;
        min_jaccard = (denom > 0.0) ? (1.0 / denom) : 1.0;
        if (min_jaccard < 0.0) min_jaccard = 0.0;
        if (min_jaccard > 1.0) min_jaccard = 1.0;
    }

    const double j = jaccard(other, min_jaccard);
    if (j <= 0.0) return std::numeric_limits<double>::infinity();
    if (j >= 1.0) return 0.0;
    const double ratio = 2.0 * j / (1.0 + j);
    return -std::log(ratio) / static_cast<double>(kmer_size_);
}

// ═══════════════════════════════════════════════════════════════════════════
// ani  –  Average Nucleotide Identity from Jaccard
//   ANI = (2J / (1+J))^(1/k)   (Ondov et al. 2016 / Mash)
// ═══════════════════════════════════════════════════════════════════════════

double FastKMV::ani(const FastKMV& other) const {
    const double j = jaccard(other);
    if (j <= 0.0) return 0.0;
    if (j >= 1.0) return 1.0;
    return std::pow(2.0 * j / (1.0 + j),
                    1.0 / static_cast<double>(kmer_size_));
}

// ═══════════════════════════════════════════════════════════════════════════
// merge  –  sorted merge of two bottom-k lists, take k smallest distinct
// ═══════════════════════════════════════════════════════════════════════════

FastKMV FastKMV::merge(const FastKMV& other) const {
    assert(k_ == other.k_ && kmer_size_ == other.kmer_size_);
    ensureSorted();
    other.ensureSorted();

    FastKMV ret(k_, kmer_size_, seed_);

    const uint64_t* a = vals_.get();
    const uint64_t* b = other.vals_.get();
    uint32_t ia = 0, ib = 0;

    while (ret.size_ < k_ && ia < size_ && ib < other.size_) {
        if (a[ia] == UINT64_MAX && b[ib] == UINT64_MAX) break;
        uint64_t next;
        if (a[ia] == b[ib]) {
            next = a[ia];
            ++ia;
            ++ib;
        } else if (a[ia] < b[ib]) {
            next = a[ia++];
        } else {
            next = b[ib++];
        }
        ret.vals_[ret.size_++] = next;
    }
    while (ret.size_ < k_ && ia < size_ && a[ia] != UINT64_MAX) {
        ret.vals_[ret.size_++] = a[ia++];
    }
    while (ret.size_ < k_ && ib < other.size_ && b[ib] != UINT64_MAX) {
        ret.vals_[ret.size_++] = b[ib++];
    }

    ret.threshold_ = (ret.size_ >= k_) ? ret.vals_[k_ - 1] : UINT64_MAX;
    ret.sorted_ = true;
    return ret;
}

// ═══════════════════════════════════════════════════════════════════════════
// printSketch
// ═══════════════════════════════════════════════════════════════════════════

void FastKMV::printSketch() const {
    ensureSorted();
    std::fprintf(stdout,
                 "FastKMV k=%u kmer=%d fill=%u/%u vals[0..19]: ",
                 k_, kmer_size_, size_, k_);
    for (uint32_t i = 0; i < k_ && i < 20; ++i)
        std::fprintf(stdout, "%016lx ", (unsigned long)vals_[i]);
    if (k_ > 20) std::fprintf(stdout, "...");
    std::fprintf(stdout, "\n");
}

#undef FKMV_VALID
