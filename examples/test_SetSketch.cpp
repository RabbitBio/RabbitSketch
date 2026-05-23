/**
 * test_SetSketch – SetSketch all-to-all with automatic method selection.
 *
 * Two modes, chosen automatically:
 *
 *   1. INVERTED INDEX (witness hashes):
 *      During sketch construction each register records the 64-bit hash of the
 *      k-mer that "won" it (the one that pushed it to its current level).
 *      These witness hashes are high-entropy identifiers — exactly like
 *      FastKMV keys — and support efficient inverted-index candidate generation.
 *      Candidates are verified with exact SIMD SetSketch Jaccard.
 *      Selected when the expected witness overlap is high enough (tight threshold).
 *
 *   2. ALL-PAIRS (SIMD batch + suffix-sum early abort):
 *      Direct O(N²) scan with aggressive per-pair early termination.
 *      Selected when the Jaccard threshold is too low for the inverted index.
 *
 * Usage:
 *   exe_test_SetSketch <file_list> <dist_threshold> <threads> [output_file]
 */

#include "Sketch.h"
#include "SetSketch.h"
#include "InvertedIndex.h"
#include "common.h"
#include "kseq.h"
#include "phmap.h"

#include <zlib.h>
#include <sys/time.h>
#include <err.h>
#include <omp.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <sys/stat.h>
#include <unistd.h>
#include <climits>

using namespace std;

KSEQ_INIT(gzFile, gzread)

static inline double mash_distance_from_jaccard(double jaccard, int kmerSize) {
    if (jaccard <= 0.0) return 1.0;
    if (jaccard >= 1.0) return 0.0;
    const double p = (2.0 * jaccard) / (1.0 + jaccard);
    return (p > 0.0) ? (-std::log(p) / static_cast<double>(kmerSize)) : 1.0;
}

static inline double min_jaccard_from_mash_distance(double mashDist, int kmerSize) {
    if (mashDist <= 0.0) return 1.0;
    const double p = std::exp(-mashDist * static_cast<double>(kmerSize));
    const double denom = 2.0 - p;
    if (denom <= 0.0) return 1.0;
    double j = p / denom;
    if (j < 0.0) return 0.0;
    if (j > 1.0) return 1.0;
    return j;
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cerr << "usage: " << argv[0]
             << " <file_list> <dist_threshold> <threads> [output_file]" << endl;
        return 1;
    }
    const string inputFile = argv[1];
    const double thres     = stod(argv[2]);
    int          nThreads  = stoi(argv[3]);
    if (nThreads < 1) nThreads = 1;
    string outPath = "res.dist.SetSketch";
    if (argc >= 5) outPath = argv[4];

    ifstream fs(inputFile);
    if (!fs) err(errno, "cannot open %s", inputFile.c_str());
    vector<string> fileList;
    { string line; while (getline(fs, line)) if (!line.empty()) fileList.push_back(line); }
    const int N = (int)fileList.size();
    cerr << "===== total files: " << N << "  (SetSketch)" << endl;

    // ── Sketch parameters ────────────────────────────────────────────────────
    static const int BITS = 13;
    static const int KMER_SIZE = 32;
    static const int WITNESS_STRIDE = 4;

    Sketch::SetSketch proto(BITS, 2.0, 5.0);
    const int    m      = proto.getM();
    const double factor = proto.getFactor();
    double bip_buf[64];
    memcpy(bip_buf, proto.getBaseInvPow(), 64 * sizeof(double));
    const double* bip = bip_buf;

    const int witnessesPerSketch = m / WITNESS_STRIDE;
    cerr << "registers=" << m << "  witnessStride=" << WITNESS_STRIDE
         << "  witnessKeys=" << witnessesPerSketch << endl;

    // ── Phase 1: Sketch construction + witness extraction ────────────────────
    vector<double>  sizes(N, 0.0);
    vector<uint8_t> flat_cores((size_t)N * m, 0);
    vector<vector<uint64_t>> skKeys(N);

    double t0 = get_sec();

    #pragma omp parallel for num_threads(nThreads) schedule(dynamic)
    for (int t = 0; t < N; t++) {
        Sketch::SetSketch sk(BITS, 2.0, 5.0);
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) continue;
        kseq_t* ks = kseq_init(fp);
        while (kseq_read(ks) >= 0)
            sk.update(ks->seq.s, static_cast<size_t>(ks->seq.l));
        kseq_destroy(ks);
        gzclose(fp);

        sizes[t] = sk.cardinality();
        memcpy(&flat_cores[(size_t)t * m], sk.getCore().data(), m);

        const uint64_t* wit = sk.getWitnesses().data();
        skKeys[t].reserve(witnessesPerSketch);
        for (int p = 0; p < m; p += WITNESS_STRIDE) {
            uint64_t key = wit[p];
            if (key != 0) skKeys[t].push_back(key);
        }
    }

    double t1 = get_sec();
    cerr << "sketch time: " << t1 - t0 << " s" << endl;

    // ── Sort by cardinality descending ───────────────────────────────────────
    vector<int> order(N);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return sizes[a] > sizes[b]; });

    {
        vector<string>  sF(N);
        vector<double>  sS(N, 0.0);
        vector<uint8_t> sC((size_t)N * m, 0);
        vector<vector<uint64_t>> sK(N);
        for (int ni = 0; ni < N; ++ni) {
            const int oi = order[ni];
            sF[ni] = fileList[oi];
            sS[ni] = sizes[oi];
            memcpy(&sC[(size_t)ni * m], &flat_cores[(size_t)oi * m], m);
            sK[ni] = std::move(skKeys[oi]);
        }
        fileList.swap(sF);
        sizes.swap(sS);
        flat_cores.swap(sC);
        skKeys.swap(sK);
    }

    // ── Precompute suffix sums ───────────────────────────────────────────────
    static const int TAIL_STEP = 64;
    const int nCP = m / TAIL_STEP + 1;
    vector<double> tailSums((size_t)N * nCP, 0.0);

    double tpre = get_sec();
    #pragma omp parallel for num_threads(nThreads) schedule(static)
    for (int t = 0; t < N; t++) {
        const uint8_t* core = &flat_cores[(size_t)t * m];
        double* ts = &tailSums[(size_t)t * nCP];
        ts[nCP - 1] = 0.0;
        for (int cp = nCP - 2; cp >= 0; --cp) {
            const int start = cp * TAIL_STEP;
            const int end   = min(start + TAIL_STEP, m);
            double blockSum = 0.0;
            for (int k = start; k < end; ++k)
                blockSum += bip[core[k]];
            ts[cp] = ts[cp + 1] + blockSum;
        }
    }
    cerr << "suffix sums: " << get_sec() - tpre << " s" << endl;

    // ── Decide mode ──────────────────────────────────────────────────────────
    const double minJac = min_jaccard_from_mash_distance(thres, KMER_SIZE);
    const double expectedOverlap = (double)witnessesPerSketch * minJac;

    // Switch at expectedOverlap >= 20: guarantees P(miss valid pair) < exp(-20) ≈ 2e-9
    const bool useInvIdx = (expectedOverlap >= 20.0);

    cerr << "minJac=" << minJac << "  expectedWitnessOverlap=" << expectedOverlap << endl;

    double t5 = get_sec();

    if (useInvIdx) {
        // ═════════════════════════════════════════════════════════════════════
        //  MODE A:  Inverted index (witness hashes) + exact verification
        //           Same clean flow as FastKMV: buildCSRIndex → computeDistancesExact
        // ═════════════════════════════════════════════════════════════════════
        cerr << "mode: INVERTED INDEX (witness hashes)\n";

        // Phase 2: build per-thread inverted index from (sorted) skKeys
        const int actualThreads = min(nThreads, N);
        vector<phmap::flat_hash_map<uint64_t, vector<uint32_t>>> threadIdx(actualThreads);

        double ti0 = get_sec();
        #pragma omp parallel num_threads(nThreads)
        {
            int tid = omp_get_thread_num();
            if (tid < actualThreads) {
                auto& localIdx = threadIdx[tid];
                #pragma omp for schedule(static)
                for (int t = 0; t < N; t++) {
                    for (uint64_t key : skKeys[t])
                        localIdx[key].push_back(static_cast<uint32_t>(t));
                }
            }
        }
        double ti1 = get_sec();
        cerr << "local index: " << ti1 - ti0 << " s" << endl;

        auto csrIdx = Sketch::buildCSRIndex<uint64_t>(threadIdx, nThreads);
        double ti2 = get_sec();
        cerr << "CSR build: " << ti2 - ti1 << " s" << endl;

        // ── Post-filter skKeys: drop singleton witnesses ──────────────────────
        // A singleton key (present in exactly one sketch) can never produce a
        // candidate pair — any lookup in postIdx for it returns an empty list.
        // All shared keys between any two sketches appear in ≥2 sketches and
        // therefore survive singleton removal in the CSR, so P(miss) is unchanged.
        // In practice ~75% of keys are singletons → saves ~245 MB and speeds up
        // candidate generation by the same fraction.
        {
            size_t kbefore = 0, kafter = 0;
            #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 64) \
                                     reduction(+:kbefore,kafter)
            for (int t = 0; t < N; t++) {
                kbefore += skKeys[t].size();
                vector<uint64_t> kept;
                kept.reserve(skKeys[t].size() / 5);
                for (uint64_t key : skKeys[t])
                    if (csrIdx.postIdx.count(key))
                        kept.push_back(key);
                kafter += kept.size();
                skKeys[t] = std::move(kept);
            }
            cerr << "skKeys singleton-filtered: " << kbefore << " -> " << kafter
                 << " (" << (kbefore > 0 ? 100.0 * kafter / kbefore : 0.0)
                 << "% retained, freed ~"
                 << (kbefore - kafter) * 8 / (1 << 20) << " MB)\n";
        }

        // Conservative minCommon: 8-sigma below expected
        const double sd = sqrt(expectedOverlap * (1.0 - minJac));
        const int minCommon = max(1, (int)floor(expectedOverlap - 8.0 * sd));
        cerr << "minCommon=" << minCommon << " (expected=" << expectedOverlap
             << ", 8σ=" << 8.0 * sd << ")\n";

        // Exact Jaccard lambda: captures core arrays, suffix sums, etc.
        auto exactJaccardFn = [&](int i, int j) -> double {
            const double si = sizes[i], sj = sizes[j];
            if (si > 0.0 && sj / si < minJac) return -1.0;
            return Sketch::SetSketch::jaccardFromCoresBatch(
                &flat_cores[(size_t)i * m], &flat_cores[(size_t)j * m],
                m, bip, factor, si, sj, minJac,
                &tailSums[(size_t)i * nCP], &tailSums[(size_t)j * nCP],
                TAIL_STEP);
        };

        auto minCommonFn = [minCommon](int /*nk*/) -> int { return minCommon; };

        Sketch::computeDistancesExact<uint64_t>(
            csrIdx, skKeys, fileList,
            N, KMER_SIZE, thres, exactJaccardFn, minCommonFn, outPath, nThreads);
    } else {
        // ═════════════════════════════════════════════════════════════════════
        //  MODE B:  All-pairs SIMD batch + suffix-sum early abort
        // ═════════════════════════════════════════════════════════════════════
        cerr << "mode: ALL-PAIRS (SIMD batch + suffix-sum early abort)\n";

        // Free witness keys — not needed for all-pairs mode
        { vector<vector<uint64_t>>().swap(skKeys); }

        // Resolve output path (same logic as computeDistancesExact)
        string finalPath = outPath;
        {
            struct stat st;
            if (stat(finalPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
                if (finalPath.back() != '/') finalPath += '/';
                finalPath += "res.dist";
            }
        }
        cerr << "output: " << finalPath << endl;

        const int TILE = 128;
        const int nTiles = (N + TILE - 1) / TILE;
        int progress = N / 20;
        if (progress < 1) progress = 1;

        vector<string> partPaths(nThreads);
        const int pid = static_cast<int>(getpid());
        for (int tid = 0; tid < nThreads; ++tid)
            partPaths[tid] = finalPath + ".part." + to_string(pid) + "." + to_string(tid);

        #pragma omp parallel num_threads(nThreads)
        {
            const int tid = omp_get_thread_num();
            FILE* tf = fopen(partPaths[tid].c_str(), "wb");
            if (!tf) {
                #pragma omp critical
                cerr << "ERROR: cannot open temp file: " << partPaths[tid] << "\n";
            } else {
                setvbuf(tf, nullptr, _IOFBF, 1 << 22);
                string buf;
                buf.reserve(1 << 22);

                #pragma omp for schedule(dynamic, 1)
                for (int tileI = 0; tileI < nTiles; ++tileI) {
                    const int iBeg = tileI * TILE;
                    const int iEnd = min(N, iBeg + TILE);
                    for (int tileJ = tileI; tileJ < nTiles; ++tileJ) {
                        const int jBeg = tileJ * TILE;
                        const int jEnd = min(N, jBeg + TILE);

                        if (tileJ > tileI && sizes[iEnd - 1] > 0.0) {
                            if (sizes[jBeg] / sizes[iEnd - 1] < minJac)
                                break;
                        }

                        for (int i = iBeg; i < iEnd; ++i) {
                            int jStart = (tileI == tileJ) ? max(i + 1, jBeg) : jBeg;
                            if (jStart >= jEnd) continue;

                            const uint8_t* core_i = &flat_cores[(size_t)i * m];
                            const double si = sizes[i];
                            const double* ts_i = &tailSums[(size_t)i * nCP];

                            for (int j = jStart; j < jEnd; ++j) {
                                const double sj = sizes[j];
                                const double maxJacBySize = (si > 0.0) ? (sj / si) : 0.0;
                                if (maxJacBySize < minJac) break;

                                if (j + 4 < jEnd) {
                                    __builtin_prefetch(&flat_cores[(size_t)(j + 4) * m], 0, 1);
                                    __builtin_prefetch(&tailSums[(size_t)(j + 4) * nCP], 0, 1);
                                }

                                const uint8_t* c2 = &flat_cores[(size_t)j * m];
                                const double* ts_j = &tailSums[(size_t)j * nCP];
                                double jaccard = Sketch::SetSketch::jaccardFromCoresBatch(
                                    core_i, c2, m, bip, factor, si, sj, minJac,
                                    ts_i, ts_j, TAIL_STEP);
                                if (jaccard < minJac) continue;
                                double dist = mash_distance_from_jaccard(jaccard, KMER_SIZE);

                                if (dist < thres) {
                                    char line[1024];
                                    int len = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                                                       fileList[i].c_str(), fileList[j].c_str(), dist);
                                    buf.append(line, static_cast<size_t>(len));
                                }

                                if (buf.size() > (1 << 22)) {
                                    fwrite(buf.data(), 1, buf.size(), tf);
                                    buf.clear();
                                }
                            }
                        }
                    }
                    const int iProgress = min(iEnd, N);
                    if (iProgress % progress == 0) {
                        #pragma omp critical
                        cerr << "  dist " << iProgress << " / " << N << "\n";
                    }
                }
                if (!buf.empty()) fwrite(buf.data(), 1, buf.size(), tf);
                fclose(tf);
            }
        }

        // Merge temp files
        FILE* fout = fopen(finalPath.c_str(), "wb");
        if (!fout) { cerr << "ERROR: cannot open " << finalPath << endl; return 1; }
        setvbuf(fout, nullptr, _IOFBF, 1 << 24);
        vector<char> copyBuf(1 << 20);
        for (int tid = 0; tid < nThreads; ++tid) {
            FILE* part = fopen(partPaths[tid].c_str(), "rb");
            if (!part) continue;
            while (true) {
                size_t got = fread(copyBuf.data(), 1, copyBuf.size(), part);
                if (got == 0) break;
                fwrite(copyBuf.data(), 1, got, fout);
            }
            fclose(part);
            remove(partPaths[tid].c_str());
        }
        fclose(fout);
    }

    double t6 = get_sec();
    cerr << "dist time: " << t6 - t5 << " s" << endl;
    cerr << "total time: " << t6 - t0 << " s" << endl;

    return 0;
}
