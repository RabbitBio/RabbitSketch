/**
 * FastKMV – fast K Minimum Values sketch for genomic k-mers.
 *
 * Canonical k-mer processing: 2-bit lex rolling hash (A→0, C→1, G→2, T→3),
 * canonical = min(fwd_enc, rev_enc), then murmur3_fmix finaliser.
 * SIMD batched fmix (AVX-512 / AVX2) where available.
 * Compile with -DFASTKMV_NO_FMUX to use raw canonical>>11 keys (ablation).
 *
 * Jaccard: standard KMV two-pointer merge on sorted bottom-k keys.
 */

#ifndef _FASTKMV_H_
#define _FASTKMV_H_

#include <cstdint>
#include <limits>
#include <memory>
#include <cassert>
#include <vector>
#include "rank/RankMetadata.h"

namespace Sketch {

class FastKMV {
public:
    static constexpr uint64_t KEY_MAX = (UINT64_C(1) << 53) - 1;

    explicit FastKMV(uint32_t k = 1024,
                     int      kmer_size = 21,
                     uint64_t seed = 42);

    /** Construct from a validated, sorted bottom-K register state. */
    static FastKMV fromRegisters(uint32_t k,
                                 int kmer_size,
                                 uint64_t seed,
                                 const std::vector<uint64_t>& registers);

    /** Hash pipeline compiled into fastkmv.cpp (including ablation builds). */
    static Rank::HashProfile runtimeHashProfile() noexcept;

    ~FastKMV() = default;
    FastKMV(const FastKMV&);
    FastKMV& operator=(FastKMV other);
    FastKMV(FastKMV&&) = default;
    FastKMV& operator=(FastKMV&&) = default;

    void update(const char* seq, uint64_t length);

    /** Seal the sketch.  Query operations remain valid; later updates fail. */
    void finalize();
    bool isSealed() const noexcept { return sealed_; }

    /**
     * KMV Jaccard on two sorted bottom-k lists: count the overlap in the
     * k smallest distinct values of the union.
     *
     * Sketches may have different K values; both are compared at min(K1,K2),
     * which is the exact bottom-k projection to their common resolution.
     * Seed, k-mer semantics, and hash profile must match.
     *
     * @param other        the other coordinated FastKMV sketch.
     * @param min_jaccard  target Jaccard threshold for early-abort. If
     *                     > 0, the intersection loop is aborted as soon
     *                     as the achievable match count is guaranteed to
     *                     fall below ceil(min_jaccard * k); the function
     *                     then returns 0.0 (a value strictly below any
     *                     positive threshold). With the default value
     *                     0.0 no pruning is applied and the exact KMV
     *                     Jaccard is returned.
     */
    double jaccard(const FastKMV& other, double min_jaccard = 0.0) const;

    /**
     * Mash distance D = -ln(2J/(1+J)) / k, derived from KMV Jaccard.
     *
     * @param other         the other sketch (must share k and kmer_size).
     * @param max_distance  target distance threshold. If finite, it is
     *                      mapped to a minimum Jaccard via the inverse
     *                      Mash formula and handed to jaccard() so that
     *                      pairs that cannot meet the threshold abort
     *                      before the intersection finishes; aborted
     *                      pairs return +infinity (strictly above any
     *                      finite threshold). With the default value
     *                      +infinity no pruning is applied and the
     *                      exact Mash distance is returned.
     */
    double distance(const FastKMV& other,
                    double max_distance =
                        std::numeric_limits<double>::infinity()) const;

    /**
     * KMV cardinality estimate: (k-1) * KEY_MAX / tau_k.
     * Returns the exact count when the sketch is not yet full.
     */
    double cardinality() const;

    /**
     * Containment of *this in other: |A ∩ B| / |A|.
     * Uses the cardinality-based formula:
     *   C(A⊆B) = J * (|A| + |B|) / (|A| * (1 + J))
     */
    double containment(const FastKMV& other) const;

    /**
     * Average Nucleotide Identity estimated from Jaccard similarity.
     * ANI = (2J / (1+J))^(1/k)  where k is the k-mer size stored at build time.
     */
    double ani(const FastKMV& other) const;

    /** Return an exact bottom-k prefix projection without original sequence. */
    FastKMV project(uint32_t target_k) const;

    /** Merge at min(K1,K2); incompatible rank spaces raise invalid_argument. */
    FastKMV merge(const FastKMV& other) const;

    const uint64_t* getRegisters() const { ensureSorted(); return vals_.get(); }
    uint32_t getK()        const { return k_; }
    uint32_t getM()        const { return k_; }
    int      getKmerSize() const { return kmer_size_; }
    uint64_t getSeed()     const { return seed_; }
    uint32_t size()        const { ensureSorted(); return size_; }
    Rank::RankMetadata metadata() const;

    void printSketch() const;

private:
    void addHash(uint64_t h);
    void insertKey(uint64_t key);
    void compactify() const;
    void ensureSorted() const;

    uint32_t k_;
    int      kmer_size_;
    uint64_t seed_;
    uint32_t buf_cap_;

    std::unique_ptr<uint64_t[]> vals_;
    mutable uint32_t size_;
    mutable uint64_t threshold_;
    mutable bool     sorted_;
    bool             sealed_;
};

} // namespace Sketch

#endif // _FASTKMV_H_
