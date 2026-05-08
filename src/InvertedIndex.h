/**
 * InvertedIndex.h – generic inverted-index module for sketch-based all-to-all
 * distance computation.
 *
 * Header-only, templated on key type (uint32_t or uint64_t).
 *
 * Workflow:
 *   1. Caller builds per-thread local inverted indices during sketch construction.
 *   2. buildCSRIndex()  – 64-shard parallel merge + singleton removal + CSR flatten.
 *   3. computeDistances() – posting-list traversal with stamp/epoch counting,
 *      caller-provided Jaccard lambda, minCommon pruning, and Mash distance output.
 */

#ifndef _INVERTED_INDEX_H_
#define _INVERTED_INDEX_H_

#include "phmap.h"
#include "common.h"

#include <omp.h>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <climits>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>

namespace Sketch {

// ─────────────────────────────────────────────────────────────────────────────
// Data structure holding the CSR-packed inverted index.
// ─────────────────────────────────────────────────────────────────────────────
template<typename KeyT>
struct InvertedIndex {
    struct PostRange { size_t off; uint32_t cnt; };

    phmap::flat_hash_map<KeyT, PostRange> postIdx;
    std::vector<uint32_t>                 csrPosts;
    size_t                                totalPostings = 0;
};

// ─────────────────────────────────────────────────────────────────────────────
// buildCSRIndexFromKeys – fast path bypassing per-thread phmap inserts.
//
//   1. Per-thread per-bucket counts of (key, sid) pairs (1024 high-bit buckets)
//   2. Compute per-thread per-bucket write offsets (parallel-safe scatter)
//   3. Scatter pairs into bucket-grouped flat array
//   4. Sort each bucket by key (parallel across buckets)
//   5. Per-bucket scan emits posting lists into csrPosts; build phmap postIdx.
//
// Output is a phmap-backed InvertedIndex equivalent to buildCSRIndex's, so
// the downstream computeDistances / computeDistancesExact code is unchanged.
// ─────────────────────────────────────────────────────────────────────────────
template<typename KeyT>
InvertedIndex<KeyT> buildCSRIndexFromKeys(
    const std::vector<std::vector<KeyT>>& skKeys,
    int N,
    int nThreads)
{
    InvertedIndex<KeyT> idx;
    if (N <= 0) return idx;

    constexpr int BK_BITS  = 10;
    constexpr int BK       = 1 << BK_BITS;
    constexpr int KEY_BITS = static_cast<int>(sizeof(KeyT) * 8);
    constexpr int SHIFT    = KEY_BITS - BK_BITS;

    struct KS { KeyT key; uint32_t sid; };

    const int T = std::min(nThreads, N);
    if (T <= 0) return idx;

    // seqLo[tid] .. seqLo[tid+1] – contiguous sequence range owned by thread tid
    std::vector<int> seqLo(T + 1, 0);
    {
        const int base = N / T;
        const int rem  = N % T;
        int p = 0;
        for (int tid = 0; tid < T; ++tid) {
            seqLo[tid] = p;
            p += base + (tid < rem ? 1 : 0);
        }
        seqLo[T] = N;
    }

    // ── Step 1: per-thread per-bucket counts ────────────────────────────────
    std::vector<size_t> cnt(static_cast<size_t>(T) * BK, 0);
    #pragma omp parallel num_threads(T)
    {
        const int tid = omp_get_thread_num();
        if (tid < T) {
            size_t* lc = &cnt[static_cast<size_t>(tid) * BK];
            const int lo = seqLo[tid], hi = seqLo[tid + 1];
            for (int t = lo; t < hi; ++t)
                for (KeyT k : skKeys[t])
                    ++lc[static_cast<size_t>(k) >> SHIFT];
        }
    }

    // ── Step 2: bucket totals + per-thread per-bucket write offsets ─────────
    std::vector<size_t> bktSize(BK, 0);
    for (int b = 0; b < BK; ++b) {
        size_t s = 0;
        for (int tid = 0; tid < T; ++tid)
            s += cnt[static_cast<size_t>(tid) * BK + b];
        bktSize[b] = s;
    }
    std::vector<size_t> bktOff(BK + 1, 0);
    for (int b = 0; b < BK; ++b) bktOff[b + 1] = bktOff[b] + bktSize[b];
    const size_t totalPairs = bktOff[BK];
    if (totalPairs == 0) return idx;

    std::vector<size_t> off(static_cast<size_t>(T) * BK, 0);
    for (int b = 0; b < BK; ++b) {
        size_t cur = bktOff[b];
        for (int tid = 0; tid < T; ++tid) {
            off[static_cast<size_t>(tid) * BK + b] = cur;
            cur += cnt[static_cast<size_t>(tid) * BK + b];
        }
    }
    std::vector<size_t>().swap(cnt);

    // ── Step 3: scatter (key, sid) pairs into bucket-grouped flat array ─────
    std::vector<KS> pairs(totalPairs);
    #pragma omp parallel num_threads(T)
    {
        const int tid = omp_get_thread_num();
        if (tid < T) {
            std::vector<size_t> cur(BK);
            std::memcpy(cur.data(), &off[static_cast<size_t>(tid) * BK],
                        BK * sizeof(size_t));
            const int lo = seqLo[tid], hi = seqLo[tid + 1];
            for (int t = lo; t < hi; ++t) {
                const uint32_t sid = static_cast<uint32_t>(t);
                for (KeyT k : skKeys[t]) {
                    int b = static_cast<int>(static_cast<size_t>(k) >> SHIFT);
                    pairs[cur[b]++] = KS{k, sid};
                }
            }
        }
    }
    std::vector<size_t>().swap(off);

    // ── Step 4: sort each bucket by key ────────────────────────────────────
    #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
    for (int b = 0; b < BK; ++b) {
        std::sort(pairs.data() + bktOff[b], pairs.data() + bktOff[b + 1],
                  [](const KS& a, const KS& bb) { return a.key < bb.key; });
    }

    // ── Step 5: per-bucket count of unique non-singleton keys + post total ──
    std::vector<size_t> bktUnique(BK, 0);
    std::vector<size_t> bktPosts (BK, 0);
    #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
    for (int b = 0; b < BK; ++b) {
        size_t u = 0, p = 0;
        const size_t lo = bktOff[b], hi = bktOff[b + 1];
        size_t i = lo;
        while (i < hi) {
            size_t j = i + 1;
            while (j < hi && pairs[j].key == pairs[i].key) ++j;
            const size_t run = j - i;
            if (run > 1) { ++u; p += run; }
            i = j;
        }
        bktUnique[b] = u;
        bktPosts [b] = p;
    }
    size_t totalUnique = 0, totalPosts = 0;
    std::vector<size_t> postOff(BK + 1, 0);
    for (int b = 0; b < BK; ++b) {
        totalUnique += bktUnique[b];
        postOff[b]   = totalPosts;
        totalPosts  += bktPosts[b];
    }
    postOff[BK] = totalPosts;
    std::cerr << "merge index: " << totalUnique << " unique keys" << std::endl;
    std::cerr << "singleton removal: -> " << totalUnique << std::endl;

    // ── Step 6: fill csrPosts in parallel; postIdx insertion serial ────────
    idx.totalPostings = totalPosts;
    idx.csrPosts.resize(totalPosts);
    idx.postIdx.reserve(totalUnique);

    #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
    for (int b = 0; b < BK; ++b) {
        const size_t lo = bktOff[b], hi = bktOff[b + 1];
        size_t out = postOff[b];
        size_t i = lo;
        while (i < hi) {
            size_t j = i + 1;
            while (j < hi && pairs[j].key == pairs[i].key) ++j;
            const size_t run = j - i;
            if (run > 1) {
                for (size_t k = i; k < j; ++k)
                    idx.csrPosts[out++] = pairs[k].sid;
            }
            i = j;
        }
    }

    // postIdx is small after singleton removal → serial insert is fast.
    for (int b = 0; b < BK; ++b) {
        const size_t lo = bktOff[b], hi = bktOff[b + 1];
        size_t out = postOff[b];
        size_t i = lo;
        while (i < hi) {
            size_t j = i + 1;
            while (j < hi && pairs[j].key == pairs[i].key) ++j;
            const size_t run = j - i;
            if (run > 1) {
                idx.postIdx[pairs[i].key] =
                    typename InvertedIndex<KeyT>::PostRange{
                        out, static_cast<uint32_t>(run)};
                out += run;
            }
            i = j;
        }
    }

    std::cerr << "CSR flatten: " << totalPosts << " postings" << std::endl;
    return idx;
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase 2: merge thread-local maps → singleton removal → CSR flatten
// ─────────────────────────────────────────────────────────────────────────────
template<typename KeyT>
InvertedIndex<KeyT> buildCSRIndex(
    std::vector<phmap::flat_hash_map<KeyT, std::vector<uint32_t>>>& threadIdx,
    int nThreads)
{
    constexpr int NUM_SHARDS = 64;
    const KeyT SHARD_MASK = static_cast<KeyT>(NUM_SHARDS - 1);
    const int tCount = static_cast<int>(threadIdx.size());

    std::vector<phmap::flat_hash_map<KeyT, std::vector<uint32_t>>> invShards(NUM_SHARDS);

    // Phase 2.1: thread-local re-shard (no locks)
    {
        std::vector<std::vector<phmap::flat_hash_map<KeyT, std::vector<uint32_t>>>> localShards(
            tCount, std::vector<phmap::flat_hash_map<KeyT, std::vector<uint32_t>>>(NUM_SHARDS));

        #pragma omp parallel for num_threads(nThreads) schedule(static)
        for (int tid = 0; tid < tCount; ++tid)
        {
            for (auto& [key, vec] : threadIdx[tid]) {
                int shard = static_cast<int>(key & SHARD_MASK);
                auto& shardMap = localShards[tid][shard];
                auto [it, ins] = shardMap.try_emplace(key, std::move(vec));
                if (!ins) {
                    auto& dst = it->second;
                    dst.reserve(dst.size() + vec.size());
                    dst.insert(dst.end(), vec.begin(), vec.end());
                }
            }
            phmap::flat_hash_map<KeyT, std::vector<uint32_t>>().swap(threadIdx[tid]);
        }
        threadIdx.clear();

        // Phase 2.2: shard merge across threads
        #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
        for (int s = 0; s < NUM_SHARDS; s++) {
            auto& merged = invShards[s];
            for (int tid = 0; tid < tCount; ++tid) {
                auto& src = localShards[tid][s];
                for (auto& [key, vec] : src) {
                    auto [it, ins] = merged.try_emplace(key, std::move(vec));
                    if (!ins) {
                        auto& dst = it->second;
                        dst.reserve(dst.size() + vec.size());
                        dst.insert(dst.end(), vec.begin(), vec.end());
                    }
                }
                phmap::flat_hash_map<KeyT, std::vector<uint32_t>>().swap(src);
            }
        }
    }

    size_t totalUnique = 0;
    for (auto& sh : invShards) totalUnique += sh.size();
    std::cerr << "merge index: " << totalUnique << " unique keys" << std::endl;

    // Singleton removal + posting count in one pass
    size_t totalAfter = 0;
    size_t totalPostings = 0;
    #pragma omp parallel for num_threads(nThreads) reduction(+:totalAfter,totalPostings)
    for (int s = 0; s < NUM_SHARDS; s++) {
        for (auto it = invShards[s].begin(); it != invShards[s].end(); ) {
            if (it->second.size() <= 1) it = invShards[s].erase(it);
            else {
                totalPostings += it->second.size();
                ++it;
                totalAfter++;
            }
        }
    }
    std::cerr << "singleton removal: " << totalUnique << " -> " << totalAfter << std::endl;

    // Flatten into CSR — two passes to parallelise the bulk memcpy:
    //   Pass 1 (serial):   compute per-key offset, fill postIdx (fast hash ops).
    //   Pass 2 (parallel): parallel memcpy per shard into csrPosts.
    // postIdx is written only in Pass 1 (serial) so no data race.
    // csrPosts writes are shard-disjoint so no synchronisation needed in Pass 2.
    InvertedIndex<KeyT> idx;
    idx.totalPostings = totalPostings;
    idx.postIdx.reserve(totalAfter);
    idx.csrPosts.resize(totalPostings);

    // Pass 1: assign offsets + populate postIdx (serial, purely hash-map ops)
    std::vector<size_t> shardOff(NUM_SHARDS, 0);
    {
        size_t cursor = 0;
        for (int s = 0; s < NUM_SHARDS; s++) {
            shardOff[s] = cursor;
            for (auto& [key, vec] : invShards[s]) {
                idx.postIdx[key] = {cursor, static_cast<uint32_t>(vec.size())};
                cursor += vec.size();
            }
        }
    }

    // Pass 2: parallel memcpy — each shard owns a disjoint region of csrPosts
    #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
    for (int s = 0; s < NUM_SHARDS; s++) {
        size_t off = shardOff[s];
        for (auto& [key, vec] : invShards[s]) {
            std::memcpy(&idx.csrPosts[off], vec.data(),
                        vec.size() * sizeof(uint32_t));
            off += vec.size();
        }
        phmap::flat_hash_map<KeyT, std::vector<uint32_t>>().swap(invShards[s]);
    }

    std::cerr << "CSR flatten: " << totalPostings << " postings" << std::endl;
    return idx;
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
template<typename KeyT, typename JaccardFn, typename MinCommonFn>
void computeDistances(
    const InvertedIndex<KeyT>&              idx,
    const std::vector<std::vector<KeyT>>&   skKeys,
    const std::vector<int>&                 sketchSizes,
    const std::vector<std::string>&         fileList,
    int                                     N,
    int                                     kmerSize,
    double                                  maxDist,
    JaccardFn                               jaccardFn,
    MinCommonFn                             minCommonFn,
    const std::string&                      outputPath,
    int                                     nThreads)
{
    const bool useMashDist = (kmerSize > 0);
    const uint32_t* csrPtr = idx.csrPosts.data();
    double minJacByDist = 1.0 - maxDist;
    if (useMashDist) {
        const double p_exp = std::exp(-static_cast<double>(kmerSize) * maxDist);
        minJacByDist = p_exp / (2.0 - p_exp);
    }

    // Size-ratio pruning (only meaningful for Mash distance mode)
    const double radio = useMashDist
        ? 2.0 * std::exp(maxDist * (kmerSize - 1)) - 1.0
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
        // uint16_t isect: max common ≤ min(s_i, s_j) ≤ K (sketch size).
        // FastKMV / ProbMinHash default K = 1024; SetSketch witnessesPerSketch
        // = m/4 = 1024 for ssBits=13.  Halving the per-thread footprint
        // (vs int) doubles L2 cache density for random-write hot path and
        // is the largest single contributor to perf scaling vs N.
        std::vector<uint16_t> isect(N, 0);
        std::vector<int>      stamp(N, 0);
        int ep = 0;
        std::vector<int> cand;
        cand.reserve(4096);
        std::string buf;
        buf.reserve(1 << 24);

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
                    __builtin_prefetch(&keys[ki + 1], 0, 1);
                auto it = idx.postIdx.find(keys[ki]);
                if (__builtin_expect(it == idx.postIdx.end(), 0)) continue;
                const uint32_t* pl   = csrPtr + it->second.off;
                const uint32_t  plSz = it->second.cnt;
                for (uint32_t pi = 0; pi < plSz; pi++) {
                    int j = static_cast<int>(pl[pi]);
                    if (j <= i) continue;
                    if (__builtin_expect(stamp[j] != ep, 1)) {
                        stamp[j] = ep;
                        isect[j] = 1;
                        cand.push_back(j);
                    } else {
                        ++isect[j];
                    }
                }
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

                if (jac > minJacByDist) {
                    double dist = 1.0 - jac;
                    if (useMashDist) {
                        dist = (jac >= 1.0) ? 0.0
                            : -std::log(2.0 * jac / (1.0 + jac))
                              / static_cast<double>(kmerSize);
                    }
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
template<typename KeyT, typename ExactJaccardFn, typename MinCommonFn>
void computeDistancesExact(
    const InvertedIndex<KeyT>&              idx,
    const std::vector<std::vector<KeyT>>&   skKeys,
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
    const uint32_t* csrPtr = idx.csrPosts.data();

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
                    __builtin_prefetch(&keys[ki + 1], 0, 1);
                auto it = idx.postIdx.find(keys[ki]);
                if (__builtin_expect(it == idx.postIdx.end(), 0)) continue;
                const uint32_t* pl   = csrPtr + it->second.off;
                const uint32_t  plSz = it->second.cnt;
                for (uint32_t pi = 0; pi < plSz; pi++) {
                    int j = static_cast<int>(pl[pi]);
                    if (j <= i) continue;
                    if (__builtin_expect(stamp[j] != ep, 1)) {
                        stamp[j] = ep;
                        isect[j] = 1;
                        cand.push_back(j);
                    } else {
                        ++isect[j];
                    }
                }
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
