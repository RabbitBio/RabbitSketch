/**
 * SetSketch - Extreme-optimized: O(1) integer-only add with global lower-bound filter,
 * SIMD-accelerated update(), and SIMD-gather cardinality/distance.
 *
 * Hot path: 1 comparison against cached global threshold → skip >95% of hashes.
 * Register distribution: P(K ≤ k) = exp(-a * base^(-k)).
 */
#ifndef _SETSKETCH_H_
#define _SETSKETCH_H_

#include <cstdint>
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

namespace Sketch {

class SetSketch {
public:
  SetSketch(int np = 14, double base = 2.0, double a = 5.0, int kmerlen = 32,
            bool track_witnesses = false);
  ~SetSketch() = default;

  void update(char* seq);
  void update(char* seq, size_t len);
  SetSketch merge(const SetSketch& other) const;
  double cardinality() const;

  /// Default Jaccard estimator: joint MLE from Ertl 2021 ("estimateJointNew").
  /// In our benchmarks this yields 2–10× lower MAE than inclusion-exclusion
  /// across all practical sequence lengths (≥10 kb) and is the SetSketch
  /// paper's main contribution. Internally falls back to inclusion-exclusion
  /// for degenerate corner cases (extremely short input vs. large sketch).
  double jaccard_index(const SetSketch& other) const;

  /// Joint MLE estimator. Identical to `jaccard_index()` — kept for explicit
  /// callers that want to make the estimator choice unambiguous.
  double jaccard_index_mle(const SetSketch& other) const;

  /// Classical inclusion-exclusion Jaccard estimator
  /// (HLL-style: J = (|A|+|B|-|A∪B|)/|A∪B|). Fast but inherits HLL's variance
  /// amplification at high J. Available for backward compatibility and for
  /// head-to-head comparisons; new code should use `jaccard_index()`.
  double jaccard_index_inclexcl(const SetSketch& other) const;

  double distance(const SetSketch& other) const { return 1.0 - jaccard_index(other); }

  /// Containment of *this in other: |A ∩ B| / |A|.
  ///   C(A⊆B) = (|A| + |B| - |AUB|) / |A|
  double containment(const SetSketch& other) const;

  /// Average Nucleotide Identity from Jaccard similarity.
  ///   ANI = (2J / (1+J))^(1/kmer_size)
  /// @param kmer_size  k-mer length used during sketching (default 32).
  double ani(const SetSketch& other, int kmer_size = 32) const;
  const std::vector<uint8_t>& getCore() const { return core_; }
  const std::vector<uint64_t>& getWitnesses() const { return witnesses_; }
  bool tracksWitnesses() const { return track_witnesses_; }
  double equalRegisterFraction(const SetSketch& other) const;
  double distanceFiltered(const SetSketch& other,
                          double min_jaccard,
                          double prefilter_factor = 0.3) const;
  void printSketch();

  // Expose internals for inline distance computation in test harness
  double getFactor() const { return factor_; }
  const double* getBaseInvPow() const { return base_inv_pow_; }
  int getM() const { return (int)(1ULL << np_); }
  double getBase() const { return base_; }
  uint32_t getQ() const { return q_; }

  // ── Inverted index support (block-of-3 registers) ─────────────────────────
  // Individual 8-bit registers have only 256 values → too low entropy for
  // effective inverted-index filtering.  Grouping 3 adjacent registers into
  // one block key drops random collision from 1/256 to ~(1/256)^3 ≈ 6e-8.
  // Candidates passing the block-match threshold are verified with exact
  // SetSketch Jaccard — zero accuracy loss.

  /** FNV-1a hash of a block of 3 register values. */
  static uint32_t blockHash(uint32_t blockIdx,
                            uint8_t v1, uint8_t v2, uint8_t v3) {
      uint32_t h = 2166136261u;
      h ^= blockIdx; h *= 16777619u;
      h ^= v1;       h *= 16777619u;
      h ^= v2;       h *= 16777619u;
      h ^= v3;       h *= 16777619u;
      return h;
  }

  /** Fill @p keys with one uint32_t key per block of 3 registers. */
  void getBlockKeys(std::vector<uint32_t>& keys, int blockSize = 3) const;

  /** Number of full blocks for a given register count and block size. */
  static int numBlocks(int m, int blockSize = 3) { return m / blockSize; }

  /**
   * Conservative min matching-block threshold for direct-distance filter.
   * P(block match | J) ≈ J^3.  6-sigma safety margin → essentially zero
   * false negatives.
   */
  static int minMatchBlocksForDist(double maxDist, int nBlocks) {
      double minJac  = 1.0 - maxDist;
      double pBlock  = minJac * minJac * minJac;
      double expected = nBlocks * pBlock;
      double sd       = std::sqrt(expected * (1.0 - pBlock));
      return std::max(1, static_cast<int>(std::floor(expected - 6.0 * sd)));
  }

  /**
   * Exact Jaccard from two flat core arrays (SIMD accelerated).
   * Used by the inverted-index Phase 3 for candidate verification.
   */
  /// Static batch Jaccard functions (used by the all-pairs / list-allpairs
  /// hot path in `rabbitsketch`).
  ///
  /// `jaccardFromCores*` implement the **inclusion-exclusion** estimator —
  /// fast SIMD reduce + tail-sum early-abort. Use them for *screening*: which
  /// pair passes a `minJaccard` threshold.
  ///
  /// `jaccardFromCoresMLE` implements the **joint MLE** (Ertl 2021) on top of
  /// the same flat-cores layout, with the same `setsketch_mle_kernel` used by
  /// `jaccard_index_mle()`. Use it for *refinement*: re-estimate the J of
  /// pairs that survived the incl-excl screen.
  ///
  /// Typical pattern in the all-pairs hot loop:
  ///   double j = jaccardFromCoresBatch(...);          // SIMD + early-abort
  ///   if (j < minJaccard) continue;                   // 99% of pairs cut here
  ///   j = jaccardFromCoresMLE(c1, c2, m, baseInvPow,  // MLE refines the rest
  ///                            factor, c1_card, c2_card, base, q);
  static double jaccardFromCores(const uint8_t* c1, const uint8_t* c2, int m,
                                 const double* baseInvPow, double factor,
                                 double card1, double card2);
  static double jaccardFromCoresEarlyAbort(const uint8_t* c1, const uint8_t* c2, int m,
                                           const double* baseInvPow, double factor,
                                           double card1, double card2, double minJaccard);

  /// Joint MLE Jaccard estimator on a flat cores layout. Used for refinement
  /// after `jaccardFromCoresBatch` clears the screening threshold.
  /// Internally falls back to incl-excl if the joint state is too degenerate
  /// (e.g. >30% registers both empty).
  static double jaccardFromCoresMLE(const uint8_t* c1, const uint8_t* c2, int m,
                                    const double* baseInvPow, double factor,
                                    double card1, double card2,
                                    double base, uint32_t q);

  /**
   * SIMD-batched Jaccard with early abort using precomputed suffix sums.
   * Processes registers in SIMD chunks of @p tailStep, checking a tight
   * upper bound at each checkpoint.  Mathematically exact: never misses
   * a pair above threshold.
   *
   * @param tailSum1  Suffix sums for c1, length (m/tailStep + 1).
   *                  tailSum1[cp] = Σ_{k≥cp*tailStep} baseInvPow[c1[k]].
   * @param tailSum2  Same for c2.
   * @param tailStep  Registers per SIMD chunk (multiple of 16, must divide m).
   */
  static double jaccardFromCoresBatch(
      const uint8_t* c1, const uint8_t* c2, int m,
      const double* baseInvPow, double factor,
      double card1, double card2, double minJaccard,
      const double* tailSum1, const double* tailSum2,
      int tailStep);

private:
  void add_slow(uint64_t hashval);
  void recompute_min();
  double union_size(const SetSketch& other) const;
  void ensure_cardinality() const;

  std::vector<uint8_t> core_;
  std::vector<uint64_t> witnesses_;  // hash that "won" each register
  bool track_witnesses_;
  uint32_t np_;
  uint32_t q_;
  double base_;
  double a_;
  double factor_;

  // Precomputed lookup tables
  uint64_t thresholds_[64];            // CDF thresholds (stack-allocated, always in cache)
  double   base_inv_pow_[64];          // base^(-k) for k=0..63
  uint64_t mask_u_;
  uint32_t shift_;

  // Global lower-bound tracking (like HLL's implicit filter via clz)
  uint8_t  min_reg_;                   // current min register value
  uint32_t count_at_min_;              // registers still at min_reg_
  uint64_t global_thresh_;             // = thresholds_[min_reg_], cached in register

  mutable double value_;
  mutable uint8_t is_calculated_;
  int      kmerLen_;    // k-mer length used in update() (default 32)
};

} // namespace Sketch

#endif
