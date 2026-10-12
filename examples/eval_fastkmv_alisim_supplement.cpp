/**
 * FastKMV-only AliSim sweep (18 parameter combos matching eval_alisim_sweep grid).
 *
 * Outputs ALISIM_SWEEP lines (algo=fastkmv) for:
 *   k  = 11,13,15,17,19,21,23,25,27,29,31  @ size=1024
 *   size = 64,128,256,512,1024,2048,4096,8192 @ kmer=21
 *
 * Usage:
 *   ./exe_eval_fastkmv_alisim_supplement --dir PATH --model NAME [--threads T] [--branch-len BL]
 */

#include "fastkmv.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <omp.h>

namespace {

struct Combo {
    int kmer;
    uint32_t size;
};

// Same 18 points as in results/alisim_sweep.csv (bits=13 placeholder).
static const Combo COMBOS[] = {
    {11, 1024}, {13, 1024}, {15, 1024}, {17, 1024}, {19, 1024}, {21, 1024},
    {23, 1024}, {25, 1024}, {27, 1024}, {29, 1024}, {31, 1024},
    {21, 64},   {21, 128},  {21, 256},  {21, 512},
    {21, 2048}, {21, 4096}, {21, 8192},
};
static constexpr int NCOMBOS = (int)(sizeof(COMBOS) / sizeof(COMBOS[0]));
static constexpr int PLACEHOLDER_BITS = 13;

std::string read_fasta(const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) return {};
    std::string seq;
    char buf[65536];
    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == '>') continue;
        int n = (int)strlen(buf);
        while (n > 0 && (buf[n - 1] == '\n' || buf[n - 1] == '\r')) --n;
        if (n > 0) seq.append(buf, n);
    }
    fclose(f);
    return seq;
}

struct PairPaths {
    int id;
    std::string fa;
    std::string fb;
};

std::vector<PairPaths> list_pairs(const std::string& dir) {
    std::vector<PairPaths> out;
    for (int i = 0; i < 100000; ++i) {
        char pa[4096], pb[4096];
        snprintf(pa, sizeof(pa), "%s/pair_%d_A.fa", dir.c_str(), i);
        snprintf(pb, sizeof(pb), "%s/pair_%d_B.fa", dir.c_str(), i);
        FILE* f = fopen(pa, "r");
        if (!f) break;
        fclose(f);
        out.push_back({i, pa, pb});
    }
    return out;
}

double exact_jaccard_k(const std::string& a, const std::string& b, int kmer) {
    if (kmer <= 0 || kmer > 31) return 0.0;
    const uint64_t mask = (kmer >= 32) ? ~0ULL : ((1ULL << (2 * kmer)) - 1);
    const int hi = 2 * (kmer - 1);
    static const int comp[4] = {3, 2, 1, 0};

    auto collect = [&](const std::string& seq) {
        std::vector<uint64_t> v;
        v.reserve(seq.size());
        uint64_t fwd = 0, rev = 0;
        int valid = 0;
        for (unsigned char ch : seq) {
            int b2 = -1;
            switch (ch) {
                case 'A':
                case 'a':
                    b2 = 0;
                    break;
                case 'C':
                case 'c':
                    b2 = 1;
                    break;
                case 'G':
                case 'g':
                    b2 = 2;
                    break;
                case 'T':
                case 't':
                    b2 = 3;
                    break;
                default:
                    valid = 0;
                    fwd = rev = 0;
                    continue;
            }
            fwd = ((fwd << 2) | (uint64_t)b2) & mask;
            rev = (rev >> 2) | ((uint64_t)comp[b2] << hi);
            if (++valid >= kmer) v.push_back(fwd <= rev ? fwd : rev);
        }
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        return v;
    };

    auto va = collect(a);
    auto vb = collect(b);
    size_t i = 0, j = 0, inter = 0;
    while (i < va.size() && j < vb.size()) {
        if (va[i] < vb[j])
            ++i;
        else if (va[i] > vb[j])
            ++j;
        else {
            ++inter;
            ++i;
            ++j;
        }
    }
    size_t uni = va.size() + vb.size() - inter;
    return uni ? (double)inter / (double)uni : 0.0;
}

struct Agg {
    double sum_exact = 0;
    double sum_est = 0;
    double sum_ae = 0;
    double sum_se = 0;
    int n = 0;
};

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr,
                "Usage: %s --dir PATH --model NAME [--threads T] [--branch-len BL]\n",
                argv[0]);
        return 1;
    }

    std::string dir, model;
    double branch_len = -1.0;
    int threads = omp_get_max_threads();
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--dir") && i + 1 < argc) dir = argv[++i];
        else if (!strcmp(argv[i], "--model") && i + 1 < argc) model = argv[++i];
        else if (!strcmp(argv[i], "--threads") && i + 1 < argc)
            threads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--branch-len") && i + 1 < argc)
            branch_len = atof(argv[++i]);
    }
    if (dir.empty() || model.empty()) {
        fprintf(stderr, "missing --dir or --model\n");
        return 1;
    }

    auto pairs = list_pairs(dir);
    if (pairs.empty()) {
        fprintf(stderr, "# no pair_*_A.fa in %s\n", dir.c_str());
        return 2;
    }

    if (branch_len < 0) {
        char meta[4096];
        snprintf(meta, sizeof(meta), "%s/pairs_meta.tsv", dir.c_str());
        FILE* mf = fopen(meta, "r");
        if (mf) {
            char line[256];
            fgets(line, sizeof(line), mf);
            if (fgets(line, sizeof(line), mf)) {
                double bl = 0;
                sscanf(line, "%*d%lf", &bl);
                if (bl > 0) branch_len = bl;
            }
            fclose(mf);
        }
    }
    if (branch_len < 0) branch_len = 0.0;

    omp_set_num_threads(threads);

    // Pre-load all sequences so each is read from disk only once.
    const int N_PAIRS = (int)pairs.size();
    std::vector<std::string> seqA(N_PAIRS), seqB(N_PAIRS);
    for (int i = 0; i < N_PAIRS; ++i) {
        seqA[i] = read_fasta(pairs[i].fa.c_str());
        seqB[i] = read_fasta(pairs[i].fb.c_str());
    }

    // Flatten pair × combo → single loop so 64 threads stay busy even when
    // N_PAIRS is small (e.g. L8M: 20 pairs × 18 combos = 360 iterations).
    const long N_TOTAL = (long)N_PAIRS * NCOMBOS;

    // Per-thread accumulators — no critical section needed.
    std::vector<std::vector<Agg>> thread_agg(threads, std::vector<Agg>(NCOMBOS));

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (long idx = 0; idx < N_TOTAL; ++idx) {
        const int pi = (int)(idx / NCOMBOS);
        const int ci = (int)(idx % NCOMBOS);

        const std::string& sa = seqA[pi];
        const std::string& sb = seqB[pi];
        if (sa.empty() || sb.empty()) continue;

        const uint64_t La = sa.size();
        const uint64_t Lb = sb.size();
        char* ca = const_cast<char*>(sa.data());
        char* cb = const_cast<char*>(sb.data());

        const int kmer = COMBOS[ci].kmer;
        const uint32_t k = COMBOS[ci].size;
        const double ej = exact_jaccard_k(sa, sb, kmer);

        Sketch::FastKMV v1(k, kmer, 42), v2(k, kmer, 42);
        v1.update(ca, La);
        v2.update(cb, Lb);
        v1.finalize();
        v2.finalize();
        const double est = v1.jaccard(v2);
        const double ae = std::fabs(est - ej);
        const double se = (est - ej) * (est - ej);

        const int tid = omp_get_thread_num();
        thread_agg[tid][ci].sum_exact += ej;
        thread_agg[tid][ci].sum_est   += est;
        thread_agg[tid][ci].sum_ae    += ae;
        thread_agg[tid][ci].sum_se    += se;
        thread_agg[tid][ci].n         += 1;
    }

    // Reduce across threads.
    std::vector<Agg> agg(NCOMBOS);
    for (int t = 0; t < threads; ++t)
        for (int ci = 0; ci < NCOMBOS; ++ci) {
            agg[ci].sum_exact += thread_agg[t][ci].sum_exact;
            agg[ci].sum_est   += thread_agg[t][ci].sum_est;
            agg[ci].sum_ae    += thread_agg[t][ci].sum_ae;
            agg[ci].sum_se    += thread_agg[t][ci].sum_se;
            agg[ci].n         += thread_agg[t][ci].n;
        }

    for (int ci = 0; ci < NCOMBOS; ++ci) {
        if (agg[ci].n <= 0) continue;
        const int kmer = COMBOS[ci].kmer;
        const uint32_t k = COMBOS[ci].size;
        const double n = (double)agg[ci].n;
        const double mej = agg[ci].sum_exact / n;
        const double mest = agg[ci].sum_est / n;
        const double mae = agg[ci].sum_ae / n;
        const double rmse = std::sqrt(agg[ci].sum_se / n);
        printf(
            "ALISIM_SWEEP dir=%s model=%s branch_len=%.6f algo=fastkmv "
            "kmer=%d size=%u bits=%d ss_a=20.0 ss_base=2.00 bd_bbits=16 bd_s64=16 drlevel=3 "
            "n=%d mean_exact_J=%.6f mean_est_J=%.6f MAE=%.6f RMSE=%.6f\n",
            dir.c_str(), model.c_str(), branch_len, kmer, k, PLACEHOLDER_BITS,
            agg[ci].n, mej, mest, mae, rmse);
    }

    return 0;
}
