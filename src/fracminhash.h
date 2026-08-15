#ifndef RABBITSKETCH_FRACMINHASH_H
#define RABBITSKETCH_FRACMINHASH_H

#include "phmap.h"
#include "rank/RankMetadata.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Sketch {

/**
 * Coordinated threshold (scaled) MinHash for genomic set operations.
 *
 * A fingerprint is retained when it is no greater than floor(UINT64_MAX /
 * scaled).  Increasing scaled is therefore an exact, deterministic sparse
 * projection; decreasing scaled would require information that was not kept.
 */
class FracMinHash {
public:
    explicit FracMinHash(uint32_t scaled = 1000,
                         int kmer_size = 21,
                         uint64_t seed = 42);

    static FracMinHash fromHashes(uint32_t scaled,
                                  int kmer_size,
                                  uint64_t seed,
                                  std::vector<uint64_t> hashes);

    void update(const char* sequence, uint64_t length);
    void addHash(uint64_t coordinated_fingerprint);
    void finalize();
    bool isSealed() const noexcept { return sealed_; }

    FracMinHash project(uint32_t target_scaled) const;
    FracMinHash merge(const FracMinHash& other) const;

    double jaccard(const FracMinHash& other) const;
    double containment(const FracMinHash& other) const;
    double cardinality() const;
    double cardinalityAt(uint32_t target_scaled) const;
    double distance(const FracMinHash& other) const;
    double ani(const FracMinHash& other) const;

    uint32_t getScaled() const noexcept { return scaled_; }
    int getKmerSize() const noexcept { return kmer_size_; }
    uint64_t getSeed() const noexcept { return seed_; }
    uint64_t threshold() const noexcept;
    double samplingProbability() const noexcept;
    size_t size() const;
    size_t retainedAt(uint32_t target_scaled) const;
    const std::vector<uint64_t>& hashes() const;
    size_t memoryBytes() const;
    Rank::RankMetadata metadata() const;

private:
    uint32_t scaled_;
    int kmer_size_;
    uint64_t seed_;
    bool buffered_build_;
    mutable phmap::flat_hash_set<uint64_t> build_hashes_;
    mutable std::vector<uint64_t> build_buffer_;
    mutable std::vector<uint64_t> hashes_;
    mutable bool sealed_ = false;

    void requireMutable(const char* operation) const;
    void reserveForUpdate(uint64_t length);
    void retainFingerprint(uint64_t fingerprint);
    void ensureFinalized() const;
    void requireCompatible(const FracMinHash& other,
                           const char* operation) const;
};

} // namespace Sketch

#endif // RABBITSKETCH_FRACMINHASH_H
