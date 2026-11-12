/**
 * test_FastKMV – FastKMV all-to-all via inverted index (zero false negatives).
 *
 * Builds an inverted index from each bottom-k hash value to the genomes that
 * contain it, then computes pairwise intersection counts via posting-list
 * traversal.
 *
 * Usage:
 *   exe_test_FastKMV <file_list> <dist_threshold> <threads> [output_file]
 */

#include "fastkmv.h"
#include "InvertedIndex.h"
#include "common.h"
#include "kseq.h"

#include <zlib.h>
#include <sys/time.h>
#include <err.h>
#include <omp.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cerr << "usage: " << argv[0]
             << " <file_list> <dist_threshold> <threads> [output_file]" << endl;
        return 1;
    }
    const string inputFile  = argv[1];
    const double maxDist    = stod(argv[2]);
    int          nThreads   = stoi(argv[3]);
    if (nThreads < 1) nThreads = 1;
    const string outPath = (argc >= 5) ? argv[4] : "res.dist.FastKMV";

    ifstream fs(inputFile);
    if (!fs) err(errno, "cannot open %s", inputFile.c_str());
    vector<string> fileList;
    { string line; while (getline(fs, line)) if (!line.empty()) fileList.push_back(line); }
    const int N = (int)fileList.size();
    cerr << "===== total files: " << N << "  (FastKMV)" << endl;

    static const uint32_t K     = 1024;
    static const int      KSIZE = 21;
    static const uint64_t SEED  = 42;

    // ── Phase 1: Sketch + native index keys ──────────────────────────────
    vector<int> sketchSizes(N);
    vector<vector<uint64_t>> skKeys(N);

    const int actualThreads = max(1, min(nThreads, N));

    double t0 = get_sec();

    #pragma omp parallel for num_threads(actualThreads) schedule(dynamic)
    for (int t = 0; t < N; t++) {
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) continue;
        kseq_t* ks = kseq_init(fp);

        Sketch::FastKMV sk(K, KSIZE, SEED);
        while (kseq_read(ks) >= 0)
            sk.update(ks->seq.s, ks->seq.l);

        kseq_destroy(ks);
        gzclose(fp);

        const uint64_t* regs = sk.getRegisters();
        uint32_t sz = sk.size();
        sketchSizes[t] = static_cast<int>(sz);
        skKeys[t].assign(regs, regs + sz);
    }

    double t1 = get_sec();
    cerr << "sketch time: " << t1 - t0 << " s" << endl;

    // ── Phase 2: Build CSR inverted index ────────────────────────────────
    auto csrIdx = Sketch::buildCSRIndexFromKeys<uint64_t>(skKeys, N, actualThreads);
    double t2 = get_sec();
    cerr << "build CSR index: " << t2 - t1 << " s" << endl;

    // ── Phase 3: Distance via inverted index ─────────────────────────────

    // Derive minimum Jaccard from Mash distance threshold
    const double p_exp  = exp(-static_cast<double>(KSIZE) * maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const int minCommon = max(1, static_cast<int>(ceil(minJac * K)));
    cerr << "pruning: minCommon=" << minCommon << "/" << K
         << "  (minJac=" << minJac << ", mashD<" << maxDist
         << ", k=" << KSIZE << ")" << endl;

    auto setJaccard = [](int common, int s0, int s1) -> double {
        int denom = s0 + s1 - common;
        if (denom <= 0) return 0.0;
        return static_cast<double>(common) / denom;
    };

    auto minCommonFn = [minCommon](int /*s0*/) -> int {
        return minCommon;
    };

    double t3 = get_sec();
    Sketch::computeDistances<uint64_t>(
        csrIdx, skKeys, sketchSizes, fileList,
        N, KSIZE, maxDist, setJaccard, minCommonFn, outPath, nThreads);
    double t4 = get_sec();

    cerr << "dist time: " << t4 - t3 << " s" << endl;
    cerr << "total time: " << t4 - t0 << " s" << endl;
    cerr << "output: " << outPath << endl;

    return 0;
}
