/**
 * test_HLL – HyperLogLog all-to-all with inverted-index witness candidates.
 *
 * REPLACES the previous LSH band-streaming implementation, which suffered
 * from MAX_BUCKET skips (false negatives) and slow O(B²) candidate explosion
 * inside large buckets. The new approach mirrors test_SetSketch.cpp:
 *
 *   1. Each HyperLogLog is constructed with witness tracking enabled. Every
 *      register update also records the 64-bit hash that "won" the register
 *      (the k-mer with the longest leading-zero run for that bucket).
 *   2. Subsample witnesses with WITNESS_STRIDE=4 (one per 4 registers).
 *      For bits=13 → 2048 witness keys per sketch (similar to FastKMV K).
 *      Two sketches sharing a register-i witness ⇒ shared k-mer ∈ A∩B,
 *      so witness-overlap is a high-quality candidate generator.
 *   3. Sort by cardinality descending so the size-ratio bound becomes a
 *      `break` instead of `continue` in any inner loop.
 *   4. Bucket scatter/sort into a CSR index over witness hashes, singleton-filter,
 *      then verify each candidate with HyperLogLog::distance() (Ertl joint MLE).
 *
 *   Mode A (expectedOverlap >= 20)  → inverted index + exact verify.
 *   Mode B (otherwise)              → in-memory O(N²) all-pairs over the
 *                                     same already-built HLL vector with
 *                                     size-ratio break-pruning.
 *
 * Usage:
 *   exe_test_HLL <file_list> <dist_threshold> <threads> [output_file]
 *
 * Output (matches old test_HLL output schema):
 *   res_dir/res.dist.HLL    or [output_file] if given
 *   format:    file_a \t file_b \t mash_distance
 */

#include "Sketch.h"
#include "InvertedIndex.h"
#include "common.h"
#include "kseq.h"

#include <zlib.h>
#include <sys/time.h>
#include <sys/stat.h>
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
#include <memory>

using namespace std;

KSEQ_INIT(gzFile, gzread)

// HLL gives Jaccard, mash distance is derived (k-mer size hard-coded as in
// the original test_HLL: HLL update() uses KMERLEN=32 internally).
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
    if (j < 0.0) j = 0.0;
    if (j > 1.0) j = 1.0;
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

    string outPath = "res_dir/res.dist.HLL";
    if (argc >= 5) outPath = argv[4];

    // The HLL update() encodes 32-mers; mash distance and pruning use this k.
    static const int KMER_SIZE = 32;
    static const int BITS      = 13;
    static const int M         = 1 << BITS;

    static const int WITNESS_STRIDE = 4;
    const int witnessesPerSketch = M / WITNESS_STRIDE;

    ifstream fs(inputFile);
    if (!fs) err(errno, "cannot open %s", inputFile.c_str());
    vector<string> fileList;
    { string line; while (getline(fs, line)) if (!line.empty()) fileList.push_back(line); }
    const int N = static_cast<int>(fileList.size());
    cerr << "===== total files: " << N << " (HyperLogLog, inverted-index)" << endl;
    cerr << "registers=" << M << "  witnessStride=" << WITNESS_STRIDE
         << "  witnessKeys=" << witnessesPerSketch
         << "  bits=" << BITS << endl;

    // ── Phase 1: parallel sketch + witness extraction ──────────────────────
    vector<unique_ptr<Sketch::HyperLogLog>> vhll(N);
    vector<double>           sizes(N, 0.0);
    vector<vector<uint64_t>> skKeys(N);

    double t0 = get_sec();
    #pragma omp parallel for num_threads(nThreads) schedule(dynamic)
    for (int t = 0; t < N; ++t) {
        auto sk = std::make_unique<Sketch::HyperLogLog>(BITS, /*track_witnesses=*/true);
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) { vhll[t] = std::move(sk); continue; }
        kseq_t* ks = kseq_init(fp);
        while (kseq_read(ks) >= 0) sk->update(ks->seq.s);
        kseq_destroy(ks);
        gzclose(fp);

        sizes[t] = sk->cardinality();
        const uint64_t* wit = sk->getWitnesses().data();
        skKeys[t].reserve(witnessesPerSketch);
        for (int p = 0; p < M; p += WITNESS_STRIDE)
            if (wit[p] != 0) skKeys[t].push_back(wit[p]);

        vhll[t] = std::move(sk);
    }
    double t1 = get_sec();
    cerr << "sketch + witness extract time: " << t1 - t0 << " s" << endl;

    // ── Sort by cardinality descending + reorder all per-genome arrays ─────
    {
        vector<int> order(N);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(),
             [&](int a, int b) { return sizes[a] > sizes[b]; });

        vector<unique_ptr<Sketch::HyperLogLog>> sH(N);
        vector<string>                          sF(N);
        vector<double>                          sS(N, 0.0);
        vector<vector<uint64_t>>                sK(N);
        for (int ni = 0; ni < N; ++ni) {
            const int oi = order[ni];
            sH[ni] = std::move(vhll[oi]);
            sF[ni] = fileList[oi];
            sS[ni] = sizes[oi];
            sK[ni] = std::move(skKeys[oi]);
        }
        vhll.swap(sH);
        fileList.swap(sF);
        sizes.swap(sS);
        skKeys.swap(sK);
    }

    // ── Mode decision ───────────────────────────────────────────────────────
    const double minJac = min_jaccard_from_mash_distance(thres, KMER_SIZE);
    const double expectedOverlap = static_cast<double>(witnessesPerSketch) * minJac;
    const bool   useInvIdx = (expectedOverlap >= 20.0);
    cerr << "minJac=" << minJac
         << "  expectedWitnessOverlap=" << expectedOverlap
         << "  mashD<" << thres << "  k=" << KMER_SIZE << endl;

    // Make the output directory exist if outPath is the default
    {
        struct stat st;
        if (stat("res_dir", &st) != 0) system("mkdir -p res_dir");
    }

    double t5 = get_sec();

    if (useInvIdx) {
        cerr << "mode: INVERTED INDEX (HLL witness hashes)" << endl;

        double ti0 = get_sec();
        auto csrIdx = Sketch::buildCSRIndexFromKeys<uint64_t>(skKeys, N, nThreads);
        double ti2 = get_sec();
        cerr << "CSR build: " << ti2 - ti0 << " s" << endl;

        // Singleton-filter: drop keys present in only one sketch.
        {
            size_t kbefore = 0, kafter = 0;
            #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 64) \
                                     reduction(+:kbefore,kafter)
            for (int t = 0; t < N; ++t) {
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
                 << (kbefore - kafter) * 8 / (1 << 20) << " MB)" << endl;
        }

        const double sd = std::sqrt(std::max(0.0, expectedOverlap * (1.0 - minJac)));
        const int minCommon = max(1, static_cast<int>(std::floor(expectedOverlap - 8.0 * sd)));
        cerr << "minCommon=" << minCommon << " (expected=" << expectedOverlap
             << ", 8sigma=" << 8.0 * sd << ")" << endl;

        auto exactJaccardFn = [&](int i, int j) -> double {
            const double si = sizes[i], sj = sizes[j];
            if (si > 0.0 && sj / si < minJac) return -1.0;
            const double d = vhll[i]->distance(*vhll[j]);
            return 1.0 - d;
        };
        auto minCommonFn = [minCommon](int) { return minCommon; };

        Sketch::computeDistancesExact<uint64_t>(
            csrIdx, skKeys, fileList,
            N, KMER_SIZE, thres, exactJaccardFn, minCommonFn,
            outPath, nThreads);

    } else {
        // ───────────── MODE B: O(N²) all-pairs over already-built HLLs ─────
        cerr << "mode: ALL-PAIRS (sorted, size-ratio break-prune)" << endl;
        { vector<vector<uint64_t>>().swap(skKeys); }

        vector<string> bufs(nThreads);
        int progress = N / 20;
        if (progress < 1) progress = 1;

        #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
        for (int i = 0; i < N; ++i) {
            const int tid = omp_get_thread_num();
            string& buf = bufs[tid];
            const double si = sizes[i];
            for (int j = i + 1; j < N; ++j) {
                const double sj = sizes[j];
                if (si > 0.0 && sj / si < minJac) break;

                const double dist = vhll[i]->distance(*vhll[j]);
                if (dist < thres) {
                    char line[1024];
                    int len = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                        fileList[i].c_str(), fileList[j].c_str(), dist);
                    buf.append(line, static_cast<size_t>(len));
                }
            }
            if (i % progress == 0) {
                #pragma omp critical
                cerr << "  dist " << i << " / " << N << endl;
            }
        }

        FILE* fp = fopen(outPath.c_str(), "w");
        if (!fp) err(errno, "cannot open output: %s", outPath.c_str());
        setvbuf(fp, nullptr, _IOFBF, 1 << 22);
        for (auto& b : bufs) fwrite(b.data(), 1, b.size(), fp);
        fclose(fp);
        cerr << "output: " << outPath << endl;
    }

    double t6 = get_sec();
    cerr << "dist time: "  << t6 - t5 << " s" << endl;
    cerr << "total time: " << get_sec() - t0 << " s" << endl;
    return 0;
}
