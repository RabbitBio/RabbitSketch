/**
 * ProbMinHash4 sketch for probability Jaccard similarity estimation.
 *
 * Based on: Otmar Ertl, "ProbMinHash – A Class of Locality-Sensitive Hash
 * Algorithms for the (Probability) Jaccard Similarity", IEEE TKDE 2020.
 * https://arxiv.org/abs/1911.00675
 *
 * Optimizations vs. baseline:
 *  1. Ziggurat exponential distribution  (~5x faster than -log(u))
 *  2. SIMD batched murmur3_fmix hashing  (AVX-512: 8-lane; AVX2: 8-lane emulated)
 *  3. Early return in addHashFromRng BEFORE perm_.reset()
 *  4. Cached class members in hot loop locals
 *  5. Global-max pre-filter to skip most addHashFromRng calls
 *  6. SIMD Jaccard / merge  (AVX-512 / AVX2 / scalar)
 *  7. PermStream version overflow protection
 *  8. No strlen(): caller passes length directly (eliminates one full scan)
 *  9. TED parameters packed into TedParam[] (better cache locality)
 * 10. tracker_ leaves ARE the sketch registers (no duplicate regs_ storage)
 * 11. O(m) build_from_leaves() for copy / merge (vs. O(m log m) on-line updates)
 * 12. cur_max loaded once per batch (register); refreshed only after real update
 * 13. bool[8] lane_valid; all-N batch skipped via branchless OR check
 * 14. Optional single-fmix fast path: compile with -DPMH_FAST_HASH
 *
 * Interface mirrors HyperLogLog / SetSketch in this library:
 *   update(seq, length) – ingest a DNA/RNA sequence (k-mer rolling hash)
 *   updateWeighted(...) – same with per–k-mer or uniform positive weights (weighted
 *                        multiset / probability-mass extension; see Ertl TKDE 2020)
 *   addHash(h, w)       – insert one hashed element with weight w (default 1)
 *   distance(other)     – probability Jaccard distance  [0, 1]
 *   jaccard(other)      – probability Jaccard similarity [0, 1]
 *   merge(other)        – element-wise min merge (returns new sketch)
 */

#ifndef _PROBMH_H_
#define _PROBMH_H_

#include <cstdint>
#include <cmath>
#include <limits>
#include <memory>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cstring>

namespace Sketch {

// ── Internal: tournament-tree max-value tracker ───────────────────────────
// Maintains the current maximum over m registers; isUpdatePossible() lets
// ProbMinHash4 skip elements that can no longer improve any register.
// Leaf values (indices 0..m-1) double as the sketch registers themselves,
// eliminating the need for a separate regs_ array.
class ProbMHMaxTracker {
public:
    explicit ProbMHMaxTracker(uint32_t m);
    void   reset(double infinity);
    bool   update(uint32_t idx, double value);
    bool   isUpdatePossible(double value) const;
    double getMax() const;

    // Direct access to the m leaf values (= sketch registers).
    double*       leaves()       noexcept { return v_.get(); }
    const double* leaves() const noexcept { return v_.get(); }

    // O(m) bottom-up rebuild of internal nodes from already-set leaves.
    // Call after bulk-writing leaves (copy / merge) instead of m individual
    // update() calls, which would cost O(m log m).
    void build_from_leaves();

    friend void swap(ProbMHMaxTracker& a, ProbMHMaxTracker& b) noexcept;

private:
    uint32_t m_;
    uint32_t lastIdx_;
    std::unique_ptr<double[]> v_;   // [0..m-1] leaves, [m..2m-2] internal nodes
};

void swap(ProbMHMaxTracker& a, ProbMHMaxTracker& b) noexcept;

// ── Internal: lazy Fisher-Yates permutation stream ────────────────────────
// SoA (val_ + ver_arr_) for better cache locality when reading version in next().
class ProbMHPermStream {
public:
    explicit ProbMHPermStream(uint32_t m);
    void     reset();
    uint32_t next(uint64_t& rng_state);   // advances rng_state in-place
    // Release the two internal arrays (val_ / ver_arr_).  Call after all
    // update() calls are done; next()/reset() must not be called afterward.
    void     clear() noexcept { val_.reset(); ver_arr_.reset(); }
    friend void swap(ProbMHPermStream& a, ProbMHPermStream& b) noexcept;

private:
    uint32_t m_;
    uint32_t idx_;
    uint32_t ver_;
    std::unique_ptr<uint32_t[]> val_;     // logical value at position
    std::unique_ptr<uint32_t[]> ver_arr_; // version stamp
};

void swap(ProbMHPermStream& a, ProbMHPermStream& b) noexcept;

// ── ProbMinHash4 sketch ────────────────────────────────────────────────────
class ProbMinHash4 {
public:
    /**
     * @param m         number of sketch registers (sketch size); must be > 1
     * @param kmer_size k-mer length for sequence hashing (default 21)
     * @param seed      64-bit seed for wyhash-based RNG (default 42)
     * @param max_L     max register updates per element (Route C truncation).
     *                  0 = unlimited (original ProbMinHash4 behaviour).
     *                  1/2/4/… = Top-L truncation (faster, slight accuracy loss).
     */
    explicit ProbMinHash4(uint32_t m = 1024,
                          int      kmer_size = 21,
                          uint64_t seed = 42,
                          uint32_t max_L = 0);

    ~ProbMinHash4() = default;
    ProbMinHash4(const ProbMinHash4&);
    ProbMinHash4& operator=(ProbMinHash4 other);
    ProbMinHash4(ProbMinHash4&&) = default;
    ProbMinHash4& operator=(ProbMinHash4&&) = default;

    /**
     * Ingest a DNA/RNA sequence.
     * Canonical k-mers (forward vs reverse-complement minimum) are used.
     * Can be called multiple times to add more sequences to the same sketch.
     * The caller must supply the sequence length (avoids an extra strlen pass).
     */
    void update(const char* seq, uint64_t length);

    /**
     * Release construction-only memory after all update() calls are done.
     * Frees ted_params_ (~40 KB/sketch) and the permutation-stream arrays
     * (~8 KB/sketch).  tracker_ (registers) and winners_ are preserved so
     * jaccard(), distance(), and getWinnerIndexKeys() remain usable.
     * Must not call update() / addHash() after finalize().
     */
    void finalize() noexcept;

    /** Every k-mer occurrence uses the same weight @p weight_each (> 0). */
    void updateWeighted(const char* seq, uint64_t length, double weight_each);

    /**
     * Per starting position: weight for the k-mer seq[pos .. pos+k-1] is
     * weight_per_kmer_start[pos].  Array must cover indices 0 .. length - k;
     * entries <= 0 skip that k-mer.  Requires length >= kmer_size.
     */
    void updateWeighted(const char* seq, uint64_t length,
                        const double* weight_per_kmer_start);

    /**
     * Fused entropy-weighted update: computes per-k-mer Shannon entropy weight
     * and processes hash registration in a single O(L) pass, without allocating
     * an intermediate weight vector.  ~2–3x faster than the separate
     * fill_kmer_entropy_weights() + updateWeighted() combination.
     *
     * @param w_min  minimum weight for degenerate / N-containing k-mers (default 0.1)
     */
    void updateEntropy(const char* seq, uint64_t length, double w_min = 0.1);

    /** Insert one element from raw k-mer hash @p h with weight @p weight (default 1). */
    void addHash(uint64_t h, double weight = 1.0);

    /**
     * Return the probability Jaccard similarity in [0, 1].
     * Compares register VALUES (exact float equality).  Correct for
     * unweighted update() / updateEntropy() where the same element
     * always produces the same priority in any sketch.
     */
    double jaccard(const ProbMinHash4& other) const;

    /**
     * Frequency-weighted Jaccard estimator via winner-hash comparison.
     *
     * Uses the winners_ array (set by addHash / addHashFromRng) to compare
     * WHICH element won each register, not the priority VALUE.  This gives
     * a correct estimate of WJ = Σmin(wA,wB)/Σmax(wA,wB) even when the
     * same k-mer has different counts (weights) in the two sketches.
     *
     * Usage:
     *   for (auto& [h, cnt] : freq_A) pmh_A.addHash(h, cnt);
     *   for (auto& [h, cnt] : freq_B) pmh_B.addHash(h, cnt);
     *   double wj = pmh_A.jaccard_weighted(pmh_B);
     *
     * Both sketches must share the same m and seed.
     */
    double jaccard_weighted(const ProbMinHash4& other) const;

    /**
     * Return Mash distance from (probability) Jaccard:
     *   D = -ln(2J/(1+J)) / kmer_size
     */
    double distance(const ProbMinHash4& other) const;

    /**
     * Weighted containment of *this in other: weighted_intersection / weight_A.
     *   C_w(A⊆B) = J_w * (w_A + w_B) / (w_A * (1 + J_w))
     * Requires total_weight() > 0.
     */
    double containment(const ProbMinHash4& other) const;

    /**
     * Average Nucleotide Identity estimated from weighted Jaccard.
     *   ANI = (2J / (1+J))^(1/kmer_size)
     */
    double ani(const ProbMinHash4& other) const;

    /**
     * Total accumulated weight inserted into this sketch (sum of all element weights).
     * For unweighted update() calls each k-mer contributes weight 1.
     */
    double total_weight() const { return total_weight_; }

    /**
     * Merge two sketches by taking element-wise minimum of hash values.
     * Both sketches must have the same parameters.
     */
    ProbMinHash4 merge(const ProbMinHash4& other) const;

    /**
     * Return pointer to the m raw register values (for serialization / LSH).
     * The pointer remains valid as long as the sketch is alive and unmodified.
     */
    const double* getRegisters() const noexcept { return tracker_.leaves(); }

    int      getKmerSize()  const { return kmer_size_; }
    uint32_t getM()         const { return m_; }
    uint32_t getMaxL()      const { return max_L_; }

    void printSketch() const;

    // ── Inverted index support ───────────────────────────────────────────────
    // Key = raw_bits(register_value) XOR (reg_idx * golden_ratio_constant).
    // Infinity registers are skipped.  Jaccard = matching_registers / M.

    /** Fill @p keys with one uint64_t inverted-index key per non-infinity register.
     *  Suitable for UNWEIGHTED update() / updateEntropy() sketches where two
     *  registers match iff their float priority values are identical. */
    void getInvertedIndexKeys(std::vector<uint64_t>& keys) const;

    /** Fill @p keys with one uint64_t inverted-index key per occupied register,
     *  encoding both the register index and the winner element identity.
     *
     *  Key = winners_[k] XOR (k * golden_ratio_constant).
     *
     *  Suitable for WEIGHTED addHash() sketches (PMH-norm).  Two keys match
     *  iff the SAME element won the SAME register in both sketches — which is
     *  exactly the condition tested by jaccard_weighted().
     *  Jaccard = matching_keys / M. */
    void getWinnerIndexKeys(std::vector<uint64_t>& keys) const;

    /** Jaccard from inverted-index common count: common / m. */
    static double jaccardFromCommon(int common, uint32_t m) {
        return static_cast<double>(common) / m;
    }

    /** Minimum common count for direct-distance threshold: ceil((1-maxDist)*m). */
    static int minCommonForDist(double maxDist, uint32_t m) {
        return std::max(1, static_cast<int>(std::ceil((1.0 - maxDist) * m)));
    }

private:
    void addHashFromRng(uint64_t rng, double weight);
    void addHashFromRngInv(uint64_t rng, double wInv);
    void updateWeightedImpl(const char* seq, uint64_t length,
                            const double* weight_per_kmer_start,
                            double uniform_weight);

    // Packed TED (Truncated-Exponential Distribution) parameters.
    // Replaces five separate arrays (boundaries_, ted_rate_, ted_c1/c2/c3_)
    // for better spatial locality: hot-loop accesses stride over one array.
    struct TedParam {
        double boundary;  // normalised cumulative boundary; [0] = 1.0
        double gap;       // boundary[i] - boundary[i-1]; used for delta (i >= 1)
        double c1, c2, c3;
    };

    // Build (or fetch from a global cache) the TED-parameter table for sketch
    // size @p m and write firstBoundary^-1 into @p out_firstBoundaryInv.  The
    // table depends only on m, so all sketches with the same m share one copy.
    // Thread-safe.  Eliminates ~40 KB malloc and ~4 K transcendentals per
    // sketch in large-N runs.
    static std::shared_ptr<const TedParam[]>
    getOrBuildTedParams(uint32_t m, double& out_firstBoundaryInv);

    uint32_t           m_;
    int                kmer_size_;
    uint64_t           seed_;
    uint32_t           max_L_;    // Route C: max updates per element (m_ = unlimited)
    double             total_weight_;  // accumulated sum of all element weights

    // Shared, immutable.  Constructor copies a pointer; no per-sketch malloc
    // or transcendental recomputation when m is identical across sketches.
    std::shared_ptr<const TedParam[]>  ted_params_;     // [m-1]
    double                             firstBoundaryInv_;

    ProbMHMaxTracker   tracker_;
    ProbMHPermStream   perm_;

    // ── Frequency-weighted Jaccard support ──────────────────────────────────
    // winners_[k] = the rng (= murmur3_fmix(element_hash, seed_)) of the
    // element that currently holds register k.  0 means "not yet set".
    // Two registers from different sketches match iff the SAME element won
    // both — even when that element has a different weight in each sketch.
    // Used by jaccard_weighted() to correctly estimate weighted Jaccard
    // WJ = Σmin(wA,wB) / Σmax(wA,wB) via addHash(hash, count).
    std::unique_ptr<uint64_t[]>  winners_;        // [m]; 0 = unset sentinel
};

// ── One-Permutation ProbMinHash (Route A) ──────────────────────────────────
//
// Each element maps to exactly ONE bucket via b(x) = hash1(x) % m, with a
// single key g(x) = Uniform(0,1) from hash2(x).  Each register keeps the
// minimum key:  S_j = min_{x: b(x)=j} g(x).
//
// Compared to ProbMinHash4:
//   • No PermStream (eliminated perm_.reset() / perm_.next())
//   • No MaxTracker tournament tree (flat array + scalar global_max)
//   • No TED sampling; one Uniform(0,1) per element
//   • O(1) per-element update (vs O(log m) tree + permutation)
//   • global_max pre-filter rejects ~99%+ elements after warm-up
//
// Trade-offs:
//   • Empty buckets possible when #elements < m  (rare for genomic data)
//   • Jaccard estimator uses count(A[j]==B[j] && finite) / m
//   • Statistical distribution differs from original ProbMinHash
//
class ProbMinHash4OP {
public:
    explicit ProbMinHash4OP(uint32_t m = 1024,
                            int      kmer_size = 21,
                            uint64_t seed = 42);

    ~ProbMinHash4OP() = default;
    ProbMinHash4OP(const ProbMinHash4OP&);
    ProbMinHash4OP& operator=(ProbMinHash4OP other);
    ProbMinHash4OP(ProbMinHash4OP&&) = default;
    ProbMinHash4OP& operator=(ProbMinHash4OP&&) = default;

    void update(const char* seq, uint64_t length);

    double jaccard(const ProbMinHash4OP& other) const;

    double distance(const ProbMinHash4OP& other) const;

    ProbMinHash4OP merge(const ProbMinHash4OP& other) const;

    const double* getRegisters() const noexcept { return regs_.get(); }
    uint32_t getM()        const { return m_; }
    int      getKmerSize() const { return kmer_size_; }
    uint32_t numEmpty()    const;

    void printSketch() const;

private:
    void addHash(uint64_t h);
    void refreshMax();

    uint32_t m_;
    int      kmer_size_;
    uint64_t seed_;

    std::unique_ptr<double[]> regs_;
    double   global_max_;
    uint32_t num_nonempty_;
};

} // namespace Sketch

#endif // _PROBMH_H_
