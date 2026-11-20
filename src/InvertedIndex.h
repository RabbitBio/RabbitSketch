/**
 * InvertedIndex.h – generic inverted-index module for sketch-based all-to-all
 * distance computation.
 *
 * Header-only, templated on integer or exact composite keys.
 *
 * Workflow:
 *   1. Caller extracts each sketch's native index keys.
 *   2. buildCSRIndexFromKeys() – bounded partitioned sorting, singleton removal,
 *      and exact delta-compressed posting lists.
 *   3. computeDistances() – posting-list traversal with stamp/epoch counting,
 *      caller-provided Jaccard lambda, minCommon pruning, and Mash distance output.
 */

#ifndef _INVERTED_INDEX_H_
#define _INVERTED_INDEX_H_

#include "PartitionedIndex.h"
#include "common.h"

#include "api/OpenMPCompat.h"
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <climits>
#include <limits>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <sys/stat.h>
#include <unistd.h>

namespace Sketch {

// Synthesized keys (e.g. BinDash position + signature) are returned by value;
// only prefetch an actual referenced element, never a temporary key.
template<class Keys>
inline void prefetchIndexKey(const Keys& keys, size_t position) {
    if constexpr (std::is_lvalue_reference<decltype(keys[position])>::value)
        __builtin_prefetch(&keys[position], 0, 1);
}

// Compatibility bucket selector. Partitioning now uses exact key ranges; a
// caller-provided hash bucket function is accepted without changing key values.
struct MixedIndexBucket {
    size_t operator()(uint64_t key) const {
        key += UINT64_C(0x9e3779b97f4a7c15);
        key = (key ^ (key >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        key = (key ^ (key >> 27)) * UINT64_C(0x94d049bb133111eb);
        return static_cast<size_t>((key ^ (key >> 31)) >> 54);
    }
};

// Ownable extracted keys can be sorted in place. Never deduplicate: some
// sketch adapters intentionally carry repeated occurrences.
template<class Key, class BucketFn = MixedIndexBucket>
InvertedIndex<Key> buildCSRIndexFromKeys(
    std::vector<std::vector<Key>>& keys, int n, int threads,
    BucketFn = {}, IndexBuildOptions options = {}, IndexBuildStats* stats = nullptr) {
    if (n < 0 || static_cast<size_t>(n) > keys.size())
        throw std::invalid_argument("invalid index sample count");
    const int workers = std::max(1, std::min(threads, std::max(1, n)));
    #pragma omp parallel for num_threads(workers) schedule(dynamic)
    for (int i = 0; i < n; ++i)
        if (!std::is_sorted(keys[i].begin(), keys[i].end()))
            std::sort(keys[i].begin(), keys[i].end());
    return buildPartitionedIndex<Key>(keys, n, threads, true, options, stats);
}

// Preserve read-only input and the old const overload. Sorted arrays use the
// same fast path; unsorted arrays use bounded rescans, with no full-size copy.
template<class Key, class BucketFn = MixedIndexBucket>
InvertedIndex<Key> buildCSRIndexFromKeys(
    const std::vector<std::vector<Key>>& keys, int n, int threads,
    BucketFn = {}, IndexBuildOptions options = {}, IndexBuildStats* stats = nullptr) {
    if (n < 0 || static_cast<size_t>(n) > keys.size())
        throw std::invalid_argument("invalid index sample count");
    bool sorted = true;
    for (int i = 0; i < n; ++i)
        sorted = sorted && std::is_sorted(keys[i].begin(), keys[i].end());
    return buildPartitionedIndex<Key>(keys, n, threads, sorted, options, stats);
}

// Compatibility for owners of existing thread-local posting maps. Sorting
// uses their existing ID arrays in place; the result uses the same codec and
// visitor as the key-array builder, without a second decoded posting array.
template<class Key>
InvertedIndex<Key> buildCSRIndex(
    std::vector<phmap::flat_hash_map<Key, std::vector<uint32_t>>>& threadIdx,
    int /*threads*/) {
    phmap::flat_hash_map<Key, std::vector<uint32_t>> merged;
    for (auto& local : threadIdx) {
        for (auto& entry : local) {
            auto found = merged.try_emplace(entry.first, std::move(entry.second));
            if (!found.second) {
                auto& target = found.first->second;
                target.insert(target.end(), entry.second.begin(), entry.second.end());
            }
        }
        phmap::flat_hash_map<Key, std::vector<uint32_t>>().swap(local);
    }
    threadIdx.clear();
    InvertedIndex<Key> result;
    size_t encoded = 0, unique = 0;
    for (auto& entry : merged) {
        auto& ids = entry.second;
        if (ids.size() < 2) continue;
        if (ids.size() > UINT32_MAX) throw std::length_error("index posting list too large");
        std::sort(ids.begin(), ids.end());
        uint32_t last = 0;
        for (uint32_t sid : ids) {
            const size_t bytes = PostingCodec::encodedSize(sid - last);
            if (bytes > SIZE_MAX - encoded) throw std::length_error("compressed index too large");
            encoded += bytes; last = sid;
        }
        ++unique; result.totalPostings += ids.size();
    }
    result.encodedPosts.resize(encoded);
    result.postIdx.reserve(unique);
    size_t cursor = 0;
    for (auto& entry : merged) {
        auto& ids = entry.second;
        if (ids.size() < 2) continue;
        result.postIdx.emplace(entry.first,
            typename InvertedIndex<Key>::PostRange{cursor, static_cast<uint32_t>(ids.size())});
        uint32_t last = 0;
        for (uint32_t sid : ids) {
            cursor += PostingCodec::encode(result.encodedPosts.data() + cursor, sid - last);
            last = sid;
        }
        std::vector<uint32_t>().swap(ids);
    }
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase 3: all-to-all distance via CSR posting traversal.
//
// JaccardFn:    double fn(int common, int s_i, int s_j)  → Jaccard similarity
// MinCommonFn:  int fn(int s_i) → minimum common count to pass distance filter
//
// FastKMV/Kssd call with:
//   jaccardFn   = [](int c, int s0, int s1){ return (double)c/(s0+s1-c); }
//   minCommonFn = [minJac](int s0){ return max(1,(int)ceil(minJac*s0)); }
//
// BinDash call with:
//   jaccardFn   = [p,n](int c,int,int){ return (c/n-p)/(1-p); }
//   minCommonFn = [minSame](int){ return minSame; }   // global constant
// ─────────────────────────────────────────────────────────────────────────────
template<typename KeyT, typename JaccardFn, typename MinCommonFn,
         typename CounterT = uint16_t,
         typename KeyRanges = std::vector<std::vector<KeyT>>,
         typename IndexT = InvertedIndex<KeyT>>
void computeDistances(
    const IndexT&                          idx,
    const KeyRanges&                        skKeys,
    const std::vector<int>&                 sketchSizes,
    const std::vector<std::string>&         fileList,
    int                                     N,
    int                                     kmerSize,
    double                                  maxDist,
    JaccardFn                               jaccardFn,
    MinCommonFn                             minCommonFn,
    const std::string&                      outputPath,
    int                                     nThreads,
    size_t                                  outputBufferLimit = 0)
{
    const bool useMashDist = (kmerSize > 0);

    // Size-ratio pruning (only meaningful for Mash distance mode)
    const double radio = useMashDist
        ? 2.0 * std::exp(maxDist * kmerSize) - 1.0
        : 1e18;

    // If outputPath is a directory, append a default filename
    std::string finalPath = outputPath;
    struct stat st;
    if (stat(finalPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
        if (finalPath.back() != '/') finalPath += '/';
        finalPath += "res.dist";
    }

    FILE* fout = std::fopen(finalPath.c_str(), "w");
    if (!fout) {
        std::cerr << "ERROR: cannot open output file: " << finalPath << std::endl;
        return;
    }
    setvbuf(fout, nullptr, _IOFBF, 1 << 24);
    std::cerr << "output: " << finalPath << std::endl;

    int progress = N / 20;
    if (progress < 1) progress = 1;

    #pragma omp parallel num_threads(nThreads)
    {
        auto&& local = localIndex(idx);
        // Fixed-size callers keep compact 16-bit counters. Threshold sketches
        // can exceed 65535 retained hashes and explicitly select uint32_t.
        std::vector<CounterT> isect(N, 0);
        std::vector<int>      stamp(N, 0);
        int ep = 0;
        std::vector<int> cand;
        cand.reserve(4096);
        std::string buf;
        buf.reserve(outputBufferLimit ? outputBufferLimit : (1 << 24));

        // dynamic,4: good balance between load-balance quality and scheduling
        // overhead. With 128 HT threads, chunk=1 causes excessive contention
        // on the atomic scheduler counter; chunk=4 cuts that 4x with minimal
        // balance loss (max idle time ≤ 4 iterations of the heaviest sequence).
        #pragma omp for schedule(dynamic, 4)
        for (int i = 0; i < N; i++) {
            const int s0 = sketchSizes[i];
            if (__builtin_expect(s0 == 0, 0)) continue;

            const int minCommon = minCommonFn(s0);

            cand.clear();
            ++ep;
            if (__builtin_expect(ep == INT_MAX, 0)) {
                std::memset(stamp.data(), 0, N * sizeof(int));
                ep = 1;
            }

            const auto& keys = skKeys[i];
            const size_t ksz = keys.size();
            for (size_t ki = 0; ki < ksz; ki++) {
                if (__builtin_expect(ki + 1 < ksz, 1))
                    prefetchIndexKey(keys, ki + 1);
                size_t offset;
                uint32_t plSz;
                if (__builtin_expect(!local.findPostingRange(keys[ki], offset, plSz), 0)) continue;
                local.forEachPosting(offset, plSz, [&](uint32_t sid) {
                    int j = static_cast<int>(sid);
                    if (j <= i) return;
                    if (__builtin_expect(stamp[j] != ep, 1)) {
                        stamp[j] = ep;
                        isect[j] = 1;
                        cand.push_back(j);
                    } else {
                        ++isect[j];
                    }
                });
            }

            for (int j : cand) {
                const int common = isect[j];
                if (common < minCommon) continue;

                const int s1 = sketchSizes[j];
                const int mn = s0 < s1 ? s0 : s1;
                const int mx = s0 > s1 ? s0 : s1;
                if (__builtin_expect(mx > radio * mn, 0)) continue;

                double jac = jaccardFn(common, s0, s1);
                if (jac < 0.0) jac = 0.0;
                if (jac > 1.0) jac = 1.0;

                double dist = 1.0 - jac;
                if (useMashDist) {
                    dist = (jac >= 1.0) ? 0.0 : (jac <= 0.0)
                        ? std::numeric_limits<double>::infinity()
                        : -std::log(2.0 * jac / (1.0 + jac))
                          / static_cast<double>(kmerSize);
                }
                // Check the actual cutoff; inverse-transform roundoff can
                // otherwise change output membership at the boundary.
                if (dist < maxDist) {
                    char line[1024];
                    int n = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                        fileList[i].c_str(), fileList[j].c_str(), dist);
                    // Opt-in bounded buffering for BinDash. Flush before an
                    // append would grow the allocation; with many workers,
                    // doubling each output string wastes several GiB. Other
                    // callers keep their existing buffering policy by default.
                    if (outputBufferLimit && !buf.empty() &&
                        buf.size() + static_cast<size_t>(n) > outputBufferLimit) {
                        #pragma omp critical
                        { std::fwrite(buf.data(), 1, buf.size(), fout); }
                        buf.clear();
                    }
                    buf.append(line, static_cast<size_t>(n));
                }
            }

            if (buf.size() > (1 << 24)) {
                #pragma omp critical
                { std::fwrite(buf.data(), 1, buf.size(), fout); }
                buf.clear();
            }
            if (i % progress == 0)
                std::cerr << "  dist " << i << " / " << N << "\n";
        }
        if (!buf.empty()) {
            #pragma omp critical
            { std::fwrite(buf.data(), 1, buf.size(), fout); }
        }
    }
    std::fclose(fout);
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase 3 (exact variant): posting-list candidate generation + caller-provided
// exact pairwise verification.
//
// Same posting-list traversal as computeDistances(), but instead of deriving
// Jaccard from the common-key count, delegates to an ExactJaccardFn that
// receives the two sketch indices and returns the true Jaccard (or negative
// to reject).  Ideal for SetSketch, where Jaccard must be computed from the
// register arrays rather than estimated from shared key counts.
//
// ExactJaccardFn:  double fn(int i, int j) → exact Jaccard, or < 0 to skip
// MinCommonFn:     int fn(int key_count_i) → minimum common keys to verify
// ─────────────────────────────────────────────────────────────────────────────
template<typename KeyT, typename ExactJaccardFn, typename MinCommonFn,
         typename KeyRanges = std::vector<std::vector<KeyT>>,
         typename IndexT = InvertedIndex<KeyT>>
void computeDistancesExact(
    const IndexT&                          idx,
    const KeyRanges&                        skKeys,
    const std::vector<std::string>&         fileList,
    int                                     N,
    int                                     kmerSize,
    double                                  maxDist,
    ExactJaccardFn                          exactJaccardFn,
    MinCommonFn                             minCommonFn,
    const std::string&                      outputPath,
    int                                     nThreads)
{
    const bool useMashDist = (kmerSize > 0);

    std::string finalPath = outputPath;
    struct stat st;
    if (stat(finalPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
        if (finalPath.back() != '/') finalPath += '/';
        finalPath += "res.dist";
    }

    FILE* fout = std::fopen(finalPath.c_str(), "w");
    if (!fout) {
        std::cerr << "ERROR: cannot open output file: " << finalPath << std::endl;
        return;
    }
    setvbuf(fout, nullptr, _IOFBF, 1 << 24);
    std::cerr << "output: " << finalPath << std::endl;

    int progress = N / 20;
    if (progress < 1) progress = 1;

    #pragma omp parallel num_threads(nThreads)
    {
        auto&& local = localIndex(idx);
        // uint16_t isect: same rationale as in computeDistances() —
        // SetSketch witnessesPerSketch defaults to 1024 (m/4 with ssBits=13),
        // safely fits in uint16_t and halves the per-thread footprint.
        std::vector<uint16_t> isect(N, 0);
        std::vector<int>      stamp(N, 0);
        int ep = 0;
        std::vector<int> cand;
        cand.reserve(4096);
        std::string buf;
        buf.reserve(1 << 24);

        #pragma omp for schedule(dynamic, 4)
        for (int i = 0; i < N; i++) {
            const int nk = static_cast<int>(skKeys[i].size());
            if (__builtin_expect(nk == 0, 0)) continue;

            const int minCommon = minCommonFn(nk);

            cand.clear();
            ++ep;
            if (__builtin_expect(ep == INT_MAX, 0)) {
                std::memset(stamp.data(), 0, N * sizeof(int));
                ep = 1;
            }

            const auto& keys = skKeys[i];
            const size_t ksz = keys.size();
            for (size_t ki = 0; ki < ksz; ki++) {
                if (__builtin_expect(ki + 1 < ksz, 1))
                    prefetchIndexKey(keys, ki + 1);
                size_t offset;
                uint32_t count;
                if (!local.findPostingRange(keys[ki], offset, count)) continue;
                local.forEachPosting(offset, count, [&](uint32_t sid) {
                    const int j = static_cast<int>(sid);
                    if (j <= i) return;
                    if (__builtin_expect(stamp[j] != ep, 1)) {
                        stamp[j] = ep;
                        isect[j] = 1;
                        cand.push_back(j);
                    } else {
                        ++isect[j];
                    }
                });
            }

            for (int j : cand) {
                if (isect[j] < minCommon) continue;

                double jac = exactJaccardFn(i, j);
                if (jac < 0.0) continue;
                if (jac > 1.0) jac = 1.0;

                double dist = 1.0 - jac;
                if (useMashDist) {
                    dist = (jac >= 1.0) ? 0.0
                        : -std::log(2.0 * jac / (1.0 + jac))
                          / static_cast<double>(kmerSize);
                }
                if (dist < maxDist) {
                    char line[1024];
                    int n = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                        fileList[i].c_str(), fileList[j].c_str(), dist);
                    buf.append(line, static_cast<size_t>(n));
                }
            }

            if (buf.size() > (1 << 24)) {
                #pragma omp critical
                { std::fwrite(buf.data(), 1, buf.size(), fout); }
                buf.clear();
            }
            if (i % progress == 0)
                std::cerr << "  dist " << i << " / " << N << "\n";
        }
        if (!buf.empty()) {
            #pragma omp critical
            { std::fwrite(buf.data(), 1, buf.size(), fout); }
        }
    }
    std::fclose(fout);
}

} // namespace Sketch

#endif // _INVERTED_INDEX_H_
