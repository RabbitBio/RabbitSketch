/**
 * test_MinHash – MinHash all-to-all via inverted index.
 *
 * Bottom-k hash values are directly the inverted-index keys; every hash is
 * a unique k-mer id so no stride/subsampling is needed (unlike HLL witnesses).
 *
 *   1. Build sketches in parallel and extract their native hash keys.
 *   2. Bucket count/scatter/sort + singleton removal → flat CSR inverted index.
 *   3. computeDistancesExact() walks posting lists with stamp/epoch counting,
 *      applying minCommon pruning before verifying with union-bottom-K Jaccard.
 *   4. Mash distance output: -log(2J/(1+J)) / k.
 *
 * Usage:
 *   exe_test_MinHash <file_list> <dist_threshold> <threads> [output_file]
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

using namespace std;

KSEQ_INIT(gzFile, gzread)

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
    string outPath = "res_dir/res.dist.MinHash";
    if (argc >= 5) outPath = argv[4];

    static const int  KMER_SIZE   = 21;
    static const int  SKETCH_SIZE = 1000;
    static const uint32_t SEED    = 42;

    ifstream fs(inputFile);
    if (!fs) err(errno, "cannot open %s", inputFile.c_str());
    vector<string> fileList;
    { string line; while (getline(fs, line)) if (!line.empty()) fileList.push_back(line); }
    const int N = static_cast<int>(fileList.size());
    cerr << "===== total files: " << N << " (MinHash inverted-index)" << endl;
    cerr << "kmerSize=" << KMER_SIZE << "  sketchSize=" << SKETCH_SIZE << endl;

    // ── Phase 1: parallel sketch construction and key extraction ──────────
    vector<vector<uint64_t>>      skKeys(N);

    const int actualThreads = max(1, min(nThreads, N));

    double t0 = get_sec();
    #pragma omp parallel for num_threads(actualThreads) schedule(dynamic)
    for (int t = 0; t < N; ++t) {
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) continue;
        kseq_t* ks = kseq_init(fp);

        Sketch::MinHash sk(KMER_SIZE, SKETCH_SIZE, SEED, /*rc=*/true);
        while (kseq_read(ks) >= 0) sk.update(ks->seq.s);
        kseq_destroy(ks);
        gzclose(fp);

        skKeys[t] = sk.getHashesSorted(); // calls finalize() internally
    }
    double t1 = get_sec();
    cerr << "sketch time: " << t1 - t0 << " s" << endl;

    // ── Phase 2: bucket scatter/sort + CSR build ──────────────────────────
    auto csrIdx = Sketch::buildCSRIndexFromKeys<uint64_t>(skKeys, N, actualThreads);
    double t2 = get_sec();
    cerr << "build CSR index: " << t2 - t1 << " s" << endl;

    // ── Pruning thresholds ─────────────────────────────────────────────────
    const double p_exp  = exp(-static_cast<double>(KMER_SIZE) * thres);
    const double minJac = p_exp / (2.0 - p_exp);
    const int minCommon = max(1,
        static_cast<int>(ceil(minJac * static_cast<double>(SKETCH_SIZE))));
    cerr << "pruning: minCommon=" << minCommon << "/" << SKETCH_SIZE
         << "  (minJac=" << minJac << ", mashD<" << thres
         << ", k=" << KMER_SIZE << ")" << endl;

    // ── Phase 3: posting-list traversal + exact verification ──────────────
    // Exact union-K Jaccard: replicate MinHash::jaccard() union-K merge.
    // Early-exit every 32 steps: if c + min(remaining_A, remaining_B) < minCommon,
    // can never pass threshold → abort merge.
    const int K = SKETCH_SIZE;
    const int mc = minCommon;
    auto exactJaccardFn = [&, K_cap = K, mc_cap = mc](int i, int j) -> double {
        const auto& hi = skKeys[i];
        const auto& hj = skKeys[j];
        const int si = (int)hi.size(), sj = (int)hj.size();
        int ii = 0, jj = 0, c = 0, denom = 0;
        while (denom < K_cap && ii < si && jj < sj) {
            if      (hi[ii] < hj[jj]) { ii++; }
            else if (hi[ii] > hj[jj]) { jj++; }
            else                      { c++; ii++; jj++; }
            denom++;
            if (__builtin_expect((denom & 31) == 0, 0) &&
                c + std::min(si - ii, sj - jj) < mc_cap) return 0.0;
        }
        if (denom < K_cap) {
            denom += (si - ii) + (sj - jj);
            if (denom > K_cap) denom = K_cap;
        }
        return (denom <= 0) ? 0.0 : (double)c / denom;
    };
    auto minCommonFn = [minCommon](int) { return minCommon; };

    {
        struct stat st;
        if (stat("res_dir", &st) != 0) system("mkdir -p res_dir");
    }

    double t3 = get_sec();
    Sketch::computeDistancesExact<uint64_t>(csrIdx, skKeys, fileList,
        N, KMER_SIZE, thres, exactJaccardFn, minCommonFn, outPath, nThreads);
    double t4 = get_sec();
    cerr << "dist time: "  << t4 - t3 << " s" << endl;
    cerr << "total time: " << t4 - t0 << " s" << endl;
    return 0;
}
