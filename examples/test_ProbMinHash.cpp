/**
 * test_ProbMinHash – ProbMinHash4 (weighted) all-to-all via inverted index.
 *
 * Each register stores a double value that acts as a unique fingerprint.
 * Key = raw-bits(double) XOR hash(register_index).  Two genomes share a key
 * iff the same k-mer produced the minimum for the same register in both —
 * probability of cross-register collision is ~2^-64.
 *
 * Jaccard = matching_registers / M.   Distance = Mash distance.
 *
 * Usage:
 *   exe_test_ProbMinHash <file_list> <dist_threshold> <threads> [output_file]
 */

#include "probmh.h"
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
#include <cstring>
#include <limits>

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
    const string outPath = (argc >= 5) ? argv[4] : "res.dist.ProbMinHash";

    ifstream fs(inputFile);
    if (!fs) err(errno, "cannot open %s", inputFile.c_str());
    vector<string> fileList;
    { string line; while (getline(fs, line)) if (!line.empty()) fileList.push_back(line); }
    const int N = (int)fileList.size();
    cerr << "===== total files: " << N << "  (ProbMinHash4)" << endl;

    static const uint32_t M     = 1024;
    static const int      KSIZE = 21;
    static const uint64_t SEED  = 42;
    const int mSize = static_cast<int>(M);

    // ── Phase 1: Sketch + native index keys ──────────────────────────────────
    vector<int> sketchSizes(N, mSize);
    vector<vector<uint64_t>> skKeys(N);

    const int actualThreads = max(1, min(nThreads, N));

    double t0 = get_sec();

    #pragma omp parallel for num_threads(actualThreads) schedule(dynamic)
    for (int t = 0; t < N; t++) {
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) continue;
        kseq_t* ks = kseq_init(fp);

        Sketch::ProbMinHash4 sk(M, KSIZE, SEED);
        while (kseq_read(ks) >= 0) {
            uint64_t L = static_cast<uint64_t>(ks->seq.l);
            if (L >= static_cast<uint64_t>(KSIZE))
                sk.updateEntropy(ks->seq.s, L);
        }
        kseq_destroy(ks);
        gzclose(fp);

        sk.getInvertedIndexKeys(skKeys[t]);
    }

    double t1 = get_sec();
    cerr << "sketch time: " << t1 - t0 << " s" << endl;

    // ── Phase 2: Build CSR inverted index ────────────────────────────────────
    auto csrIdx = Sketch::buildCSRIndexFromKeys<uint64_t>(skKeys, N, actualThreads);
    double t2 = get_sec();
    cerr << "build CSR index: " << t2 - t1 << " s" << endl;

    // ── Phase 3: Distance via inverted index ─────────────────────────────────
    // ProbMinHash: Jaccard = common / M, dist = Mash distance:
    //   D = -(1/k) * ln(2J/(1+J))
    // Convert maxDist -> minimum Jaccard threshold for pruning:
    //   J_min = exp(-kD) / (2 - exp(-kD))
    const double p_exp  = std::exp(-static_cast<double>(KSIZE) * maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const int minCommon = std::max(1, static_cast<int>(std::ceil(minJac * M)));
    cerr << "pruning: minCommon=" << minCommon << "/" << mSize
         << "  (minJac=" << minJac << ", maxMashDist=" << maxDist << ")" << endl;

    auto jaccardFn = [mSize](int common, int /*s0*/, int /*s1*/) -> double {
        return Sketch::ProbMinHash4::jaccardFromCommon(
            common, static_cast<uint32_t>(mSize));
    };

    auto minCommonFn = [minCommon](int /*s0*/) -> int {
        return minCommon;
    };

    double t3 = get_sec();
    Sketch::computeDistances<uint64_t>(
        csrIdx, skKeys, sketchSizes, fileList,
        N, KSIZE /*k>0 enables Mash distance mode*/, maxDist,
        jaccardFn, minCommonFn, outPath, nThreads);
    double t4 = get_sec();

    cerr << "dist time: " << t4 - t3 << " s" << endl;
    cerr << "total time: " << t4 - t0 << " s" << endl;

    return 0;
}
