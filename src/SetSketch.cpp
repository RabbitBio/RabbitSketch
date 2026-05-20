/**
 * SetSketch extreme-optimized implementation.
 *
 * add hot path: 1 integer compare with cached global_thresh_ → skip 95%+ hashes.
 * update():     SIMD hash + inline add with global + per-register early exit.
 * cardinality:  SIMD gather from precomputed base_inv_pow_ table.
 * union_size:   inline max + SIMD gather, zero allocation.
 */
#include "SetSketch.h"
#include "Sketch.h"
#include "MurmurHash3.h"
#include "hash_int.h"
#include <immintrin.h>
#include <algorithm>
#include <cstring>
#include <cstdio>

using namespace Sketch;

// ── AVX2 64-bit lane multiply helper (no AVX-512DQ needed) ───────────────────
#ifdef __AVX2__
static inline __m256i ss_avx2_mullo_epi64(__m256i a, __m256i b) {
    __m256i hi_a = _mm256_srli_epi64(a, 32);
    __m256i hi_b = _mm256_srli_epi64(b, 32);
    __m256i lo   = _mm256_mul_epu32(a, b);
    __m256i mid  = _mm256_add_epi64(_mm256_mul_epu32(hi_a, b),
                                     _mm256_mul_epu32(a, hi_b));
    return _mm256_add_epi64(lo, _mm256_slli_epi64(mid, 32));
}
#endif

// ── SIMD equal-register counter ───────────────────────────────────────────────
namespace {
static int setsketch_count_equal_regs(const uint8_t* __restrict__ a,
                                       const uint8_t* __restrict__ b,
                                       int n) {
  int count = 0, i = 0;
#if defined(__AVX512BW__)
  for (; i + 64 <= n; i += 64) {
    __m512i va = _mm512_loadu_si512((const void*)(a + i));
    __m512i vb = _mm512_loadu_si512((const void*)(b + i));
    count += (int)__builtin_popcountll(
        (uint64_t)_mm512_cmpeq_epi8_mask(va, vb));
  }
#endif
#if defined(__AVX2__)
  for (; i + 32 <= n; i += 32) {
    __m256i va = _mm256_loadu_si256((const __m256i*)(a + i));
    __m256i vb = _mm256_loadu_si256((const __m256i*)(b + i));
    __m256i eq = _mm256_cmpeq_epi8(va, vb);
    count += (int)__builtin_popcount((uint32_t)_mm256_movemask_epi8(eq));
  }
#endif
  for (; i < n; i++)
    count += (int)(a[i] == b[i]);
  return count;
}

// Count registers where c1[i] > c2[i] and c1[i] < c2[i] in one pass.
// Returns (n_greater, n_less). n_equal = m - n_greater - n_less.
static inline void setsketch_count_greater_less(const uint8_t* __restrict__ a,
                                                 const uint8_t* __restrict__ b,
                                                 int n,
                                                 int& n_greater,
                                                 int& n_less) {
  int g = 0, l = 0, i = 0;
#if defined(__AVX512BW__)
  for (; i + 64 <= n; i += 64) {
    __m512i va = _mm512_loadu_si512((const void*)(a + i));
    __m512i vb = _mm512_loadu_si512((const void*)(b + i));
    // unsigned a > b  iff  max(a,b) == a  AND  a != b
    uint64_t mask_eq = (uint64_t)_mm512_cmpeq_epi8_mask(va, vb);
    __m512i vmax    = _mm512_max_epu8(va, vb);
    uint64_t mask_amax = (uint64_t)_mm512_cmpeq_epi8_mask(va, vmax);
    uint64_t mask_g = mask_amax & ~mask_eq;
    uint64_t mask_l = (~mask_amax) & ~mask_eq;
    g += __builtin_popcountll(mask_g);
    l += __builtin_popcountll(mask_l);
  }
#endif
#if defined(__AVX2__)
  for (; i + 32 <= n; i += 32) {
    __m256i va = _mm256_loadu_si256((const __m256i*)(a + i));
    __m256i vb = _mm256_loadu_si256((const __m256i*)(b + i));
    __m256i veq = _mm256_cmpeq_epi8(va, vb);
    __m256i vmax = _mm256_max_epu8(va, vb);
    __m256i va_is_max = _mm256_cmpeq_epi8(va, vmax);
    __m256i vg = _mm256_andnot_si256(veq, va_is_max);
    __m256i vl = _mm256_andnot_si256(veq, _mm256_andnot_si256(va_is_max, _mm256_set1_epi8(-1)));
    g += __builtin_popcount((uint32_t)_mm256_movemask_epi8(vg));
    l += __builtin_popcount((uint32_t)_mm256_movemask_epi8(vl));
  }
#endif
  for (; i < n; i++) {
    if (a[i] > b[i]) g++;
    else if (a[i] < b[i]) l++;
  }
  n_greater = g;
  n_less    = l;
}

static inline double setsketch_sum_max_registers(const uint8_t* __restrict__ c1,
                                                 const uint8_t* __restrict__ c2,
                                                 int m,
                                                 const double* __restrict__ baseInvPow) {
  double sum = 0.0;
  int i = 0;
#if defined(__AVX512BW__) && defined(__AVX512F__)
  __m512d vsum0 = _mm512_setzero_pd();
  __m512d vsum1 = _mm512_setzero_pd();
  for (; i + 16 <= m; i += 16) {
    __m128i va = _mm_loadu_si128((const __m128i*)(c1 + i));
    __m128i vb = _mm_loadu_si128((const __m128i*)(c2 + i));
    __m128i vmax = _mm_max_epu8(va, vb);
    __m256i vidx0 = _mm256_cvtepu8_epi32(vmax);
    __m128i vmax_hi = _mm_srli_si128(vmax, 8);
    __m256i vidx1 = _mm256_cvtepu8_epi32(vmax_hi);
    vsum0 = _mm512_add_pd(vsum0, _mm512_i32gather_pd(vidx0, baseInvPow, 8));
    vsum1 = _mm512_add_pd(vsum1, _mm512_i32gather_pd(vidx1, baseInvPow, 8));
  }
  sum += _mm512_reduce_add_pd(_mm512_add_pd(vsum0, vsum1));
#elif defined(__AVX2__)
  __m256d vacc0 = _mm256_setzero_pd();
  __m256d vacc1 = _mm256_setzero_pd();
  __m256d vacc2 = _mm256_setzero_pd();
  __m256d vacc3 = _mm256_setzero_pd();
  for (; i + 16 <= m; i += 16) {
    __m128i va = _mm_loadu_si128((const __m128i*)(c1 + i));
    __m128i vb = _mm_loadu_si128((const __m128i*)(c2 + i));
    __m128i vmax = _mm_max_epu8(va, vb);

    __m128i idx0_8 = vmax;
    __m128i idx1_8 = _mm_srli_si128(vmax, 8);
    __m256i idx0_32 = _mm256_cvtepu8_epi32(idx0_8);
    __m256i idx1_32 = _mm256_cvtepu8_epi32(idx1_8);
    __m128i idx0_lo = _mm256_castsi256_si128(idx0_32);
    __m128i idx0_hi = _mm256_extracti128_si256(idx0_32, 1);
    __m128i idx1_lo = _mm256_castsi256_si128(idx1_32);
    __m128i idx1_hi = _mm256_extracti128_si256(idx1_32, 1);

    vacc0 = _mm256_add_pd(vacc0, _mm256_i32gather_pd(baseInvPow, idx0_lo, 8));
    vacc1 = _mm256_add_pd(vacc1, _mm256_i32gather_pd(baseInvPow, idx0_hi, 8));
    vacc2 = _mm256_add_pd(vacc2, _mm256_i32gather_pd(baseInvPow, idx1_lo, 8));
    vacc3 = _mm256_add_pd(vacc3, _mm256_i32gather_pd(baseInvPow, idx1_hi, 8));
  }
  __m256d vacc01 = _mm256_add_pd(vacc0, vacc1);
  __m256d vacc23 = _mm256_add_pd(vacc2, vacc3);
  __m256d vacc = _mm256_add_pd(vacc01, vacc23);
  alignas(32) double buf[4];
  _mm256_store_pd(buf, vacc);
  sum += buf[0] + buf[1] + buf[2] + buf[3];
#endif
  for (; i < m; i++) {
    uint8_t r = (c1[i] > c2[i]) ? c1[i] : c2[i];
    sum += baseInvPow[r];
  }
  return sum;
}
} // namespace

// ── Constructor: precompute threshold + base_inv_pow tables ───────────────────
SetSketch::SetSketch(int np, double base, double a, int kmerlen, bool track_witnesses)
    : np_(np), q_(62), base_(base), a_(a),
      min_reg_(0), value_(0.0), is_calculated_(0), kmerLen_(kmerlen) {
  assert(np >= 4 && np <= 16);
  assert(base > 1.0);
  assert(a > 0.0);

  const uint64_t m = 1ULL << np;
  core_.resize(m, 0);
  track_witnesses_ = track_witnesses;
  if (track_witnesses_) witnesses_.resize(m, 0);
  shift_ = 64 - np;
  mask_u_ = (shift_ >= 64) ? ~0ULL : (1ULL << shift_) - 1;

  const double log_base = std::log(base);
  factor_ = m * (base - 1.0) / (base * log_base * a);

  // base^(-k) table
  for (int k = 0; k < 64; k++)
    base_inv_pow_[k] = std::pow(base, -(double)k);

  // CDF integer thresholds: threshold[k] = exp(-a * base^(-k)) * 2^shift_
  const double scale = (double)(mask_u_ + 1);
  for (int k = 0; k <= (int)q_; k++) {
    double cdf = std::exp(-a * base_inv_pow_[k]);
    uint64_t t = (uint64_t)(cdf * scale);
    thresholds_[k] = (t > mask_u_) ? mask_u_ : t;
  }
  thresholds_[q_ + 1] = mask_u_ + 1; // sentinel

  // Global lower bound
  count_at_min_ = (uint32_t)m;
  global_thresh_ = thresholds_[0];
}

// ── Recompute min_reg_ by scanning all registers ─────────────────────────────
void SetSketch::recompute_min() {
  uint8_t mn = core_[0];
  for (size_t i = 1; i < core_.size(); i++)
    if (core_[i] < mn) mn = core_[i];
  min_reg_ = mn;
  uint32_t cnt = 0;
  for (size_t i = 0; i < core_.size(); i++)
    cnt += (core_[i] == mn);
  count_at_min_ = cnt;
  global_thresh_ = thresholds_[mn];
}

// ── add_slow: fallback used by remainder loop (not hot path) ──────────────────
void SetSketch::add_slow(uint64_t hashval) {
  const uint64_t rest = hashval & mask_u_;
  if (rest < global_thresh_) return;

  const uint32_t index = (uint32_t)(hashval >> shift_);
  const uint8_t cur = core_[index];
  if (rest < thresholds_[cur]) return;

  uint8_t k = cur + 1;
  while (k < q_ && rest >= thresholds_[k]) ++k;
  core_[index] = k;
  if (track_witnesses_) witnesses_[index] = hashval;
  is_calculated_ = 0;

  if (cur == min_reg_) {
    if (--count_at_min_ == 0)
      recompute_min();
  }
}

// ── Cardinality: SIMD gather from base_inv_pow_ table ─────────────────────────
void SetSketch::ensure_cardinality() const {
  if (is_calculated_) return;
  const size_t sz = core_.size();
  const double* __restrict__ tbl = base_inv_pow_;
  const uint8_t* __restrict__ c = core_.data();
  double sum = 0.0;

#if defined(__AVX512F__)
  size_t i = 0;
  __m512d vsum0 = _mm512_setzero_pd();
  __m512d vsum1 = _mm512_setzero_pd();
  for (; i + 16 <= sz; i += 16) {
    __m128i vidx8_0 = _mm_loadl_epi64((__m128i*)(c + i));
    __m256i vidx32_0 = _mm256_cvtepu8_epi32(vidx8_0);
    vsum0 = _mm512_add_pd(vsum0, _mm512_i32gather_pd(vidx32_0, tbl, 8));

    __m128i vidx8_1 = _mm_loadl_epi64((__m128i*)(c + i + 8));
    __m256i vidx32_1 = _mm256_cvtepu8_epi32(vidx8_1);
    vsum1 = _mm512_add_pd(vsum1, _mm512_i32gather_pd(vidx32_1, tbl, 8));
  }
  sum = _mm512_reduce_add_pd(_mm512_add_pd(vsum0, vsum1));
  for (; i < sz; i++) sum += tbl[c[i]];
#elif defined(__AVX2__)
  // 8 elements per iteration: zero-extend bytes → int32 indices → gather doubles
  size_t i = 0;
  __m256d vacc0 = _mm256_setzero_pd();
  __m256d vacc1 = _mm256_setzero_pd();
  for (; i + 8 <= sz; i += 8) {
    __m128i vidx8   = _mm_loadl_epi64((const __m128i*)(c + i));
    __m256i vidx32  = _mm256_cvtepu8_epi32(vidx8);
    __m128i vlo     = _mm256_castsi256_si128(vidx32);
    __m128i vhi     = _mm256_extracti128_si256(vidx32, 1);
    vacc0 = _mm256_add_pd(vacc0, _mm256_i32gather_pd(tbl, vlo, 8));
    vacc1 = _mm256_add_pd(vacc1, _mm256_i32gather_pd(tbl, vhi, 8));
  }
  alignas(32) double buf[4];
  _mm256_store_pd(buf, _mm256_add_pd(vacc0, vacc1));
  sum += buf[0] + buf[1] + buf[2] + buf[3];
  for (; i < sz; i++) sum += tbl[c[i]];
#else
  for (size_t i = 0; i < sz; i++) sum += tbl[c[i]];
#endif

  value_ = (sum > 1e-300) ? factor_ / sum : 0.0;
  is_calculated_ = 1;
}

double SetSketch::cardinality() const {
  ensure_cardinality();
  return value_;
}

// ── union_size: inline max + SIMD gather, no allocation ───────────────────────
double SetSketch::union_size(const SetSketch& other) const {
  const int sz = static_cast<int>(core_.size());
  const double* __restrict__ tbl = base_inv_pow_;
  const uint8_t* __restrict__ c1 = core_.data();
  const uint8_t* __restrict__ c2 = other.core_.data();
  const double sum = setsketch_sum_max_registers(c1, c2, sz, tbl);

  return (sum > 1e-300) ? factor_ / sum : 0.0;
}

double SetSketch::jaccard_index_inclexcl(const SetSketch& other) const {
  double us = union_size(other);
  if (us <= 0.0) return 0.0;
  double c1 = cardinality();
  double c2 = other.cardinality();
  double inter = c1 + c2 - us;
  return (inter > 0.0) ? inter / us : 0.0;
}

// Default = joint-MLE estimator. Definition is below; we just forward to it.
double SetSketch::jaccard_index(const SetSketch& other) const {
  return jaccard_index_mle(other);
}

// ── Joint MLE Jaccard estimator (Ertl 2021, "estimateJointNew") ──────────────
//
// Given two SetSketch states with cardinality estimates c1, c2 and counts of
// equal / A>B / A<B registers, the log-likelihood of true Jaccard J is:
//
//   logL(J) = N_eq * log1p( log_b(1+(c2*J-c1)*z) + log_b(1+(c1*J-c2)*z) )
//           + N_gt * log( -log_b(1+(c2*J-c1)*z) )
//           + N_lt * log( -log_b(1+(c1*J-c2)*z) )
//   z = (1 - 1/b) / (c1 + c2)
//
// We maximize logL on J ∈ [0, min(c1/c2, c2/c1)] using Brent's algorithm.
// If both cardinalities are 0 or any register-extreme corner case appears,
// fall back to inclusion-exclusion (which is what the reference also does).
//
// References:
//   Ertl O. (2021). SetSketch: Filling the Gap between MinHash and HLL.
//   set-sketch-paper/c++/sketch.hpp :: estimateJointNew
namespace {

// Brent's 1-D minimizer (no-derivative, golden + parabolic interpolation).
// Searches min of f on [ax, cx]. Returns x* with |x* - true min| < tol*(|x*|+ZEPS).
template <typename F>
static double brent_minimize(F f, double ax, double cx,
                             double tol = 1e-10, int max_iter = 200) {
  const double GOLD = 0.3819660112501051;  // (3 - sqrt(5)) / 2
  const double ZEPS = 1e-12;

  double bx = ax + GOLD * (cx - ax);
  double a = ax, b = cx;
  double x = bx, w = bx, v = bx;
  double fx = f(x), fw = fx, fv = fx;
  double d = 0.0, e = 0.0;

  for (int iter = 0; iter < max_iter; ++iter) {
    double xm = 0.5 * (a + b);
    double tol1 = tol * std::fabs(x) + ZEPS;
    double tol2 = 2.0 * tol1;
    if (std::fabs(x - xm) <= tol2 - 0.5 * (b - a)) {
      return x;
    }

    bool use_golden = true;
    if (std::fabs(e) > tol1) {
      double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) p = -p;
      q = std::fabs(q);
      double etemp = e;
      e = d;
      if (std::fabs(p) < std::fabs(0.5 * q * etemp) &&
          p > q * (a - x) && p < q * (b - x)) {
        d = p / q;
        double u_test = x + d;
        if (u_test - a < tol2 || b - u_test < tol2)
          d = (xm - x >= 0.0 ? tol1 : -tol1);
        use_golden = false;
      }
    }

    if (use_golden) {
      e = (x >= xm) ? (a - x) : (b - x);
      d = GOLD * e;
    }

    double u = (std::fabs(d) >= tol1)
                   ? (x + d)
                   : (x + (d >= 0.0 ? tol1 : -tol1));
    double fu = f(u);

    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }
  return x;
}

} // namespace

// ── Shared MLE kernel (used by both the OO and from-cores entry points) ─────
//
// Inputs:
//   c1, c2       : two SetSketch register arrays of length m
//   m            : # registers
//   baseInvPow   : pre-computed table base^(-k), used by incl-excl fallback
//   factor       : SetSketch's `numRegisters / (base * logBaseDivBaseMinus1 * a)`
//   card1, card2 : SetSketch cardinality estimates of the two sketches
//   base         : SetSketch base parameter
//   q            : SetSketch q parameter (max non-saturated register value)
//
// Returns Ĵ in [0, 1]. When the joint MLE objective is degenerate (sketch
// extremely under-filled), falls back to inclusion-exclusion via the static
// `jaccardFromCores` helper to remain numerically safe.
static double setsketch_mle_kernel(const uint8_t* __restrict__ c1,
                                   const uint8_t* __restrict__ c2,
                                   int m,
                                   const double* __restrict__ baseInvPow,
                                   double factor,
                                   double card1, double card2,
                                   double base, uint32_t q) {
  if (m <= 0) return 0.0;
  if (card1 <= 0.0 || card2 <= 0.0) return 0.0;

  // 1) Joint register counts (SIMD).
  int n_gt = 0, n_lt = 0;
  setsketch_count_greater_less(c1, c2, m, n_gt, n_lt);
  const int n_eq = m - n_gt - n_lt;

  // 2) Boundary cases.
  if (n_gt + n_lt == 0) return 1.0;  // all equal → J ≈ 1
  if (n_eq == m) {
    // unreachable here (handled above), but for completeness
    return SetSketch::jaccardFromCores(c1, c2, m, baseInvPow, factor, card1, card2);
  }

  // 3) Range correction: fall back when too many registers carry no joint info
  //    (both empty or both saturated). Threshold 30%.
  const uint8_t Q = static_cast<uint8_t>(q + 1);
  int n_both_zero = 0;
  int n_both_max  = 0;
  for (int i = 0; i < m; ++i) {
    if (c1[i] == 0 && c2[i] == 0)       n_both_zero++;
    else if (c1[i] == Q && c2[i] == Q)  n_both_max++;
  }
  if ((n_both_zero + n_both_max) * 10 > m * 3) {
    return SetSketch::jaccardFromCores(c1, c2, m, baseInvPow, factor, card1, card2);
  }

  // 4) MLE objective: minimize −logL(J).
  const double inv_log_base = 1.0 / std::log(base);
  const double z = (1.0 - 1.0 / base) / (card1 + card2);
  const double j_max = (card1 >= card2) ? (card2 / card1) : (card1 / card2);
  if (j_max <= 1e-12) return 0.0;

  auto neg_logL = [&](double j) -> double {
    double log1px1 = 0.0, log1px2 = 0.0;
    bool need1 = (n_eq > 0) || (n_gt > 0);
    bool need2 = (n_eq > 0) || (n_lt > 0);
    if (need1) {
      double arg1 = (card2 * j - card1) * z;
      if (arg1 <= -1.0) return std::numeric_limits<double>::infinity();
      log1px1 = std::log1p(arg1) * inv_log_base;
    }
    if (need2) {
      double arg2 = (card1 * j - card2) * z;
      if (arg2 <= -1.0) return std::numeric_limits<double>::infinity();
      log1px2 = std::log1p(arg2) * inv_log_base;
    }

    double ret = 0.0;
    if (n_eq > 0) {
      double inner = log1px1 + log1px2;
      if (inner <= -1.0) return std::numeric_limits<double>::infinity();
      ret += n_eq * std::log1p(inner);
    }
    if (n_gt > 0) {
      if (log1px1 >= 0.0) return std::numeric_limits<double>::infinity();
      ret += n_gt * std::log(-log1px1);
    }
    if (n_lt > 0) {
      if (log1px2 >= 0.0) return std::numeric_limits<double>::infinity();
      ret += n_lt * std::log(-log1px2);
    }
    if (std::isnan(ret)) return std::numeric_limits<double>::infinity();
    return -ret;
  };

  const double EPS = 1e-9;
  double lo = EPS;
  double hi = std::max(EPS * 2.0, j_max - EPS);
  if (hi <= lo) return 0.0;

  double j_hat = brent_minimize(neg_logL, lo, hi, 1e-10, 200);
  if (std::isnan(j_hat) || j_hat < 0.0) j_hat = 0.0;
  if (j_hat > j_max)                    j_hat = j_max;
  return j_hat;
}

double SetSketch::jaccard_index_mle(const SetSketch& other) const {
  const int m = static_cast<int>(core_.size());
  if (m <= 0 || (int)other.core_.size() != m) return 0.0;
  if (np_ != other.np_ || base_ != other.base_ || a_ != other.a_) {
    // Incompatible sketches → fall back to inclusion-exclusion.
    return jaccard_index_inclexcl(other);
  }

  const double c1 = cardinality();
  const double c2 = other.cardinality();
  return setsketch_mle_kernel(core_.data(), other.core_.data(), m,
                              base_inv_pow_, factor_,
                              c1, c2, base_, q_);
}

// ── Static MLE entry point for flat-cores all-pairs hot paths ────────────────
double SetSketch::jaccardFromCoresMLE(
    const uint8_t* c1, const uint8_t* c2, int m,
    const double* baseInvPow, double factor,
    double card1, double card2,
    double base, uint32_t q)
{
  return setsketch_mle_kernel(c1, c2, m, baseInvPow, factor,
                              card1, card2, base, q);
}

SetSketch SetSketch::merge(const SetSketch& other) const {
  assert(np_ == other.np_ && base_ == other.base_ && a_ == other.a_);
  SetSketch ret(np_, base_, a_);
  for (size_t i = 0; i < core_.size(); ++i)
    ret.core_[i] = std::max(core_[i], other.core_[i]);
  ret.is_calculated_ = 0;
  ret.recompute_min();
  return ret;
}

double SetSketch::equalRegisterFraction(const SetSketch& other) const {
  const int n = static_cast<int>(core_.size());
  if (n == 0 || n != static_cast<int>(other.core_.size())) return 0.0;
  return (double)setsketch_count_equal_regs(core_.data(), other.core_.data(), n) / n;
}

double SetSketch::distanceFiltered(const SetSketch& other,
                                   double min_jaccard,
                                   double prefilter_factor) const {
  if (equalRegisterFraction(other) < min_jaccard * prefilter_factor)
    return -1.0;
  return distance(other);
}

void SetSketch::printSketch() {
  fprintf(stdout, "SetSketch core[%zu]: ", core_.size());
  for (size_t i = 0; i < core_.size() && i < 20; i++)
    fprintf(stdout, "%u ", core_[i]);
  if (core_.size() > 20) fprintf(stdout, "...");
  fprintf(stdout, "\n");
}

// ── update(seq): rolling k-mer + SIMD hash + INLINE add with global filter ───
static const uint8_t ENCODE_LUT[256] = {
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
#define ENC(c)   (ENCODE_LUT[(uint8_t)(c)])
#define COMP(e)  ((uint8_t)((e) <= 3 ? 3 - (e) : 255))
#define VALID(e) ((e) <= 3)

void SetSketch::update(char* seq) {
  update(seq, strlen(seq));
}

void SetSketch::update(char* seq, size_t len) {
  const uint64_t LENGTH = static_cast<uint64_t>(len);
  const int KMERLEN = kmerLen_;
  const uint64_t fwd_mask = (KMERLEN >= 32) ? UINT64_MAX : ((1ULL << (2 * KMERLEN)) - 1);
  if (LENGTH < (uint64_t)KMERLEN) return;

  // Cache class members in locals for the hot loop
  const uint32_t shift = shift_;
  const uint64_t mask = mask_u_;
  const uint32_t qmax = q_;
  uint8_t*       core = core_.data();
  uint64_t*      wit  = track_witnesses_ ? witnesses_.data() : nullptr;
  const uint64_t* thresh = thresholds_;
  uint8_t  loc_min_reg     = min_reg_;
  uint32_t loc_count_at_min = count_at_min_;
  uint64_t loc_global_thresh = global_thresh_;
  uint8_t  loc_is_calc = is_calculated_;
  bool     loc_min_dirty = false;

  uint64_t fwd_enc = 0, rev_enc = 0;
  int invalid_count = 0;
  for (int k = 0; k < KMERLEN; k++) {
    uint8_t ef = ENC(seq[k]);
    if (!VALID(ef)) invalid_count++;
    fwd_enc = ((fwd_enc << 2) | (VALID(ef) ? (ef & 3u) : 0u)) & fwd_mask;
    uint8_t er = VALID(ef) ? (COMP(ef) & 3u) : 0u;
    rev_enc = (rev_enc >> 2) | (static_cast<uint64_t>(er) << (2 * (KMERLEN - 1)));
  }

  const int lanes = 8;
  const uint64_t total_kmers = LENGTH - KMERLEN + 1;
  const uint64_t N = (total_kmers / lanes) * lanes;

  for (uint64_t i = 0; i < N; i += lanes) {
    uint64_t resv[8];
    bool lane_valid[8];
    for (int j = 0; j < lanes; j++) {
      lane_valid[j] = (invalid_count == 0);
      resv[j] = lane_valid[j] ? ((fwd_enc <= rev_enc) ? fwd_enc : rev_enc) : 0;

      uint8_t ef_out = ENC(seq[i + j]);
      uint8_t ef_in  = ENC(seq[i + j + KMERLEN]);
      if (!VALID(ef_out)) invalid_count--;
      if (!VALID(ef_in))  invalid_count++;
      fwd_enc = ((fwd_enc << 2) | (VALID(ef_in) ? (ef_in & 3u) : 0u)) & fwd_mask;
      uint8_t er_in = VALID(ef_in) ? (COMP(ef_in) & 3u) : 0u;
      rev_enc = (rev_enc >> 2) | ((uint64_t)er_in << (2 * (KMERLEN - 1)));
    }

    // ── SIMD hash (identical to HLL) ──────────────────────────────────────
    uint64_t hashvalv[8];
#if defined(__AVX512F__) && defined(__AVX512DQ__)
    {
      __m512i vb = _mm512_loadu_si512((void*)resv);
      __m512i vseed = _mm512_set1_epi64(42);
      __m512i va = _mm512_xor_epi64(vb, vseed);
      __m512i vtmp = _mm512_srli_epi64(va, 33);
      vb = _mm512_xor_epi64(va, vtmp);
      va = _mm512_mullo_epi64(vb, _mm512_set1_epi64(0xff51afd7ed558ccdULL));
      vtmp = _mm512_srli_epi64(va, 33);
      vb = _mm512_xor_epi64(va, vtmp);
      va = _mm512_mullo_epi64(vb, _mm512_set1_epi64(0xc4ceb9fe1a85ec53ULL));
      vtmp = _mm512_srli_epi64(va, 33);
      vb = _mm512_xor_epi64(va, vtmp);
      _mm512_storeu_si512(hashvalv, vb);
    }
#elif defined(__AVX2__)
    {
      // murmur3_fmix across 8 lanes using two 256-bit registers (4 lanes each)
      const __m256i C1   = _mm256_set1_epi64x(0xff51afd7ed558ccdLL);
      const __m256i C2   = _mm256_set1_epi64x(0xc4ceb9fe1a85ec53LL);
      const __m256i SEED = _mm256_set1_epi64x(42LL);
      __m256i vb0 = _mm256_loadu_si256((const __m256i*)resv);
      __m256i vb1 = _mm256_loadu_si256((const __m256i*)(resv + 4));
      // lanes 0-3
      __m256i va0 = _mm256_xor_si256(vb0, SEED);
      __m256i vt0 = _mm256_srli_epi64(va0, 33);
      vb0 = _mm256_xor_si256(va0, vt0);
      va0 = ss_avx2_mullo_epi64(vb0, C1);
      vt0 = _mm256_srli_epi64(va0, 33);
      vb0 = _mm256_xor_si256(va0, vt0);
      va0 = ss_avx2_mullo_epi64(vb0, C2);
      vt0 = _mm256_srli_epi64(va0, 33);
      vb0 = _mm256_xor_si256(va0, vt0);
      // lanes 4-7
      __m256i va1 = _mm256_xor_si256(vb1, SEED);
      __m256i vt1 = _mm256_srli_epi64(va1, 33);
      vb1 = _mm256_xor_si256(va1, vt1);
      va1 = ss_avx2_mullo_epi64(vb1, C1);
      vt1 = _mm256_srli_epi64(va1, 33);
      vb1 = _mm256_xor_si256(va1, vt1);
      va1 = ss_avx2_mullo_epi64(vb1, C2);
      vt1 = _mm256_srli_epi64(va1, 33);
      vb1 = _mm256_xor_si256(va1, vt1);
      _mm256_storeu_si256((__m256i*)hashvalv,       vb0);
      _mm256_storeu_si256((__m256i*)(hashvalv + 4), vb1);
    }
#else
    for (int j = 0; j < lanes; j++)
      hashvalv[j] = mc::murmur3_fmix(resv[j], 42);
#endif

    // ── SIMD index/rest extraction + global filter ────────────────────────
#if defined(__AVX512F__)
    {
      __m512i vhash = _mm512_loadu_si512((void*)hashvalv);
      __m512i vmask = _mm512_set1_epi64(mask);
      __m512i vrest = _mm512_and_epi64(vhash, vmask);

      __m512i vgt = _mm512_set1_epi64(loc_global_thresh);
      __mmask8 pass_mask = _mm512_cmpge_epu64_mask(vrest, vgt);

      if (pass_mask != 0) {
        uint64_t restv[8];
        _mm512_storeu_si512(restv, vrest);

        for (int j = 0; j < lanes; j++) {
          if (!lane_valid[j]) continue;
          if (!(pass_mask & (1 << j))) continue;

          uint32_t idx = (uint32_t)(hashvalv[j] >> shift);
          uint64_t rest = restv[j];
          uint8_t cur = core[idx];
          if (rest < thresh[cur]) continue;

          uint8_t k = cur + 1;
          while (k < qmax && rest >= thresh[k]) ++k;
          core[idx] = k;
          if (wit) wit[idx] = hashvalv[j];
          loc_is_calc = 0;

          if (cur == loc_min_reg) {
            if (!loc_min_dirty && loc_count_at_min > 0 && --loc_count_at_min == 0)
              loc_min_dirty = true;
          }
        }
      }
    }
#elif defined(__AVX2__)
    {
      const __m256i vmask2 = _mm256_set1_epi64x((int64_t)mask);
      __m256i vh0 = _mm256_loadu_si256((const __m256i*)hashvalv);
      __m256i vh1 = _mm256_loadu_si256((const __m256i*)(hashvalv + 4));
      __m256i vr0 = _mm256_and_si256(vh0, vmask2);
      __m256i vr1 = _mm256_and_si256(vh1, vmask2);

      // Unsigned >=: rest >= thresh  ⟺  !(thresh > rest)  via sign-flip trick
      const __m256i sign2 = _mm256_set1_epi64x((int64_t)0x8000000000000000ULL);
      __m256i vgt_adj = _mm256_xor_si256(_mm256_set1_epi64x((int64_t)loc_global_thresh), sign2);
      // fail bit i = 1 when thresh > rest[i] (unsigned)
      int f0 = _mm256_movemask_pd(_mm256_castsi256_pd(
                   _mm256_cmpgt_epi64(vgt_adj, _mm256_xor_si256(vr0, sign2))));
      int f1 = _mm256_movemask_pd(_mm256_castsi256_pd(
                   _mm256_cmpgt_epi64(vgt_adj, _mm256_xor_si256(vr1, sign2))));
      uint8_t pass_mask = (uint8_t)~((f1 << 4) | f0);

      if (pass_mask != 0) {
        uint64_t restv[8];
        _mm256_storeu_si256((__m256i*)restv,       vr0);
        _mm256_storeu_si256((__m256i*)(restv + 4), vr1);

        for (int j = 0; j < lanes; j++) {
          if (!lane_valid[j]) continue;
          if (!(pass_mask & (1u << j))) continue;

          uint32_t idx = (uint32_t)(hashvalv[j] >> shift);
          uint64_t rest = restv[j];
          uint8_t cur = core[idx];
          if (rest < thresh[cur]) continue;

          uint8_t k = cur + 1;
          while (k < qmax && rest >= thresh[k]) ++k;
          core[idx] = k;
          if (wit) wit[idx] = hashvalv[j];
          loc_is_calc = 0;

          if (cur == loc_min_reg) {
            if (!loc_min_dirty && loc_count_at_min > 0 && --loc_count_at_min == 0)
              loc_min_dirty = true;
          }
        }
      }
    }
#else
    // Scalar path with global filter
    for (int j = 0; j < lanes; j++) {
      if (!lane_valid[j]) continue;
      uint64_t rest = hashvalv[j] & mask;
      if (rest < loc_global_thresh) continue;

      uint32_t idx = (uint32_t)(hashvalv[j] >> shift);
      uint8_t cur = core[idx];
      if (rest < thresh[cur]) continue;

      uint8_t k = cur + 1;
      while (k < qmax && rest >= thresh[k]) ++k;
      core[idx] = k;
      if (wit) wit[idx] = hashvalv[j];
      loc_is_calc = 0;

      if (cur == loc_min_reg) {
        if (!loc_min_dirty && loc_count_at_min > 0 && --loc_count_at_min == 0)
          loc_min_dirty = true;
      }
    }
#endif
  }

  // ── Remainder loop ──────────────────────────────────────────────────────
  for (uint64_t i = N; i < total_kmers; ++i) {
    if (invalid_count == 0) {
      uint64_t res = (fwd_enc <= rev_enc) ? fwd_enc : rev_enc;
      uint64_t hashval = mc::murmur3_fmix(res, 42);
      uint64_t rest = hashval & mask;
      if (rest >= loc_global_thresh) {
        uint32_t idx = (uint32_t)(hashval >> shift);
        uint8_t cur = core[idx];
        if (rest >= thresh[cur]) {
          uint8_t k = cur + 1;
          while (k < qmax && rest >= thresh[k]) ++k;
          core[idx] = k;
          if (wit) wit[idx] = hashval;
          loc_is_calc = 0;
          if (cur == loc_min_reg) {
            if (!loc_min_dirty && loc_count_at_min > 0 && --loc_count_at_min == 0)
              loc_min_dirty = true;
          }
        }
      }
    }
    uint8_t ef_out = ENC(seq[i]);
    uint8_t ef_in  = ENC(seq[i + KMERLEN]);
    if (!VALID(ef_out)) invalid_count--;
    if (!VALID(ef_in))  invalid_count++;
    fwd_enc = ((fwd_enc << 2) | (VALID(ef_in) ? (ef_in & 3u) : 0u)) & fwd_mask;
    uint8_t er_in = VALID(ef_in) ? (COMP(ef_in) & 3u) : 0u;
    rev_enc = (rev_enc >> 2) | ((uint64_t)er_in << (2 * (KMERLEN - 1)));
  }

  // Write back locals. If min tracking was exhausted, rebuild once at end.
  if (loc_min_dirty) {
    is_calculated_ = loc_is_calc;
    recompute_min();
  } else {
    min_reg_ = loc_min_reg;
    count_at_min_ = loc_count_at_min;
    global_thresh_ = loc_global_thresh;
    is_calculated_ = loc_is_calc;
  }
}

// ── containment ──────────────────────────────────────────────────────────────
// C(this ⊆ other) = (|A| + |B| - |AUB|) / |A|
double SetSketch::containment(const SetSketch& other) const
{
  const double card_a = cardinality();
  if (card_a <= 0.0) return 0.0;
  const double card_b = other.cardinality();
  const double us     = union_size(other);
  const double inter  = card_a + card_b - us;
  return (inter > 0.0) ? inter / card_a : 0.0;
}

// ── ani ──────────────────────────────────────────────────────────────────────
// ANI = (2J / (1+J))^(1/kmer_size)   (Mash / Ondov et al. 2016)
double SetSketch::ani(const SetSketch& other, int kmer_size) const
{
  const double j = jaccard_index(other);
  if (j <= 0.0) return 0.0;
  if (j >= 1.0) return 1.0;
  return std::pow(2.0 * j / (1.0 + j), 1.0 / static_cast<double>(kmer_size));
}

// ── inverted index: block key extraction ─────────────────────────────────────
void SetSketch::getBlockKeys(std::vector<uint32_t>& keys, int blockSize) const {
    const int m = getM();
    const int nBlocks = m / blockSize;
    keys.resize(nBlocks);
    for (int b = 0; b < nBlocks; b++) {
        keys[b] = blockHash(static_cast<uint32_t>(b),
                            core_[b * 3],
                            core_[b * 3 + 1],
                            core_[b * 3 + 2]);
    }
}

// ── inverted index: exact Jaccard from flat core arrays (SIMD) ──────────────
double SetSketch::jaccardFromCores(
    const uint8_t* __restrict__ c1,
    const uint8_t* __restrict__ c2,
    int m,
    const double* __restrict__ baseInvPow,
    double factor,
    double card1, double card2)
{
    const double sum = setsketch_sum_max_registers(c1, c2, m, baseInvPow);
    if (sum <= 1e-300) return 0.0;
    double us    = factor / sum;
    double inter = card1 + card2 - us;
    return (inter > 0.0) ? inter / us : 0.0;
}

double SetSketch::jaccardFromCoresEarlyAbort(
    const uint8_t* __restrict__ c1,
    const uint8_t* __restrict__ c2,
    int m,
    const double* __restrict__ baseInvPow,
    double factor,
    double card1, double card2,
    double minJaccard)
{
    if (m <= 0) return 0.0;
    if (minJaccard <= 0.0) {
        return jaccardFromCores(c1, c2, m, baseInvPow, factor, card1, card2);
    }
    const double cardSum = card1 + card2;
    const double maxTerm = baseInvPow[0];
    double sum = 0.0;
    for (int i = 0; i < m; ++i) {
        uint8_t r = (c1[i] > c2[i]) ? c1[i] : c2[i];
        sum += baseInvPow[r];

        const int remain = m - i - 1;
        const double maxPossibleSum = sum + maxTerm * static_cast<double>(remain);
        const double maxPossibleJ = (cardSum * maxPossibleSum / factor) - 1.0;
        if (maxPossibleJ < minJaccard) return -1.0;
    }
    if (sum <= 1e-300) return 0.0;
    const double us = factor / sum;
    const double inter = cardSum - us;
    return (inter > 0.0) ? inter / us : 0.0;
}

double SetSketch::jaccardFromCoresBatch(
    const uint8_t* __restrict__ c1,
    const uint8_t* __restrict__ c2,
    int m,
    const double* __restrict__ baseInvPow,
    double factor,
    double card1, double card2,
    double minJaccard,
    const double* __restrict__ tailSum1,
    const double* __restrict__ tailSum2,
    int tailStep)
{
    if (m <= 0) return 0.0;
    if (minJaccard <= 0.0)
        return jaccardFromCores(c1, c2, m, baseInvPow, factor, card1, card2);

    const double cardSum = card1 + card2;
    const double targetSum = (1.0 + minJaccard) * factor / cardSum;
    double sum = 0.0;
    int i = 0, cp = 0;

    for (; i + tailStep <= m; i += tailStep, ++cp) {
        sum += setsketch_sum_max_registers(c1 + i, c2 + i, tailStep, baseInvPow);
        double tail = tailSum1[cp + 1];
        double t2   = tailSum2[cp + 1];
        if (t2 < tail) tail = t2;
        if (sum + tail < targetSum) return -1.0;
    }
    for (; i < m; ++i) {
        uint8_t r = (c1[i] > c2[i]) ? c1[i] : c2[i];
        sum += baseInvPow[r];
    }
    if (sum <= 1e-300) return 0.0;
    double us = factor / sum;
    double inter = cardSum - us;
    return (inter > 0.0) ? inter / us : 0.0;
}

#undef ENC
#undef COMP
#undef VALID
