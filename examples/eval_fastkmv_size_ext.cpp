/**
 * FastKMV size extension sweep: k=21 only, sketch sizes doubling from 8k to 8M.
 *
 * Sizes (registers): 8192, 16384, ..., 8388608  (11 points).
 *
 * Usage:
 *   ./exe_eval_fastkmv_size_ext --dir PATH --model NAME [--threads T] [--branch-len BL]
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

static const uint32_t SIZE_GRID[] = {
    8192u, 16384u, 32768u, 65536u, 131072u, 262144u,
    524288u, 1048576u, 2097152u, 4194304u, 8388608u,
};
static constexpr int NSIZES = (int)(sizeof(SIZE_GRID) / sizeof(SIZE_GRID[0]));
static constexpr int KMER = 21;
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
                case 'A': case 'a': b2 = 0; break;
                case 'C': case 'c': b2 = 1; break;
                case 'G': case 'g': b2 = 2; break;
                case 'T': case 't': b2 = 3; break;
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
        if (va[i] < vb[j]) ++i;
        else if (va[i] > vb[j]) ++j;
        else { ++inter; ++i; ++j; }
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

    const int N_PAIRS = (int)pairs.size();
    std::vector<std::string> seqA(N_PAIRS), seqB(N_PAIRS);
    for (int i = 0; i < N_PAIRS; ++i) {
        seqA[i] = read_fasta(pairs[i].fa.c_str());
        seqB[i] = read_fasta(pairs[i].fb.c_str());
    }

    const long N_TOTAL = (long)N_PAIRS * NSIZES;
    std::vector<std::vector<Agg>> thread_agg(threads, std::vector<Agg>(NSIZES));

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (long idx = 0; idx < N_TOTAL; ++idx) {
        const int pi = (int)(idx / NSIZES);
        const int si = (int)(idx % NSIZES);

        const std::string& sa = seqA[pi];
        const std::string& sb = seqB[pi];
        if (sa.empty() || sb.empty()) continue;

        const uint32_t k = SIZE_GRID[si];
        const double ej = exact_jaccard_k(sa, sb, KMER);

        Sketch::FastKMV v1(k, KMER, 42), v2(k, KMER, 42);
        v1.update(const_cast<char*>(sa.data()), sa.size());
        v2.update(const_cast<char*>(sb.data()), sb.size());
        v1.finalize();
        v2.finalize();
        const double est = v1.jaccard(v2);
        const double ae = std::fabs(est - ej);
        const double se = (est - ej) * (est - ej);

        const int tid = omp_get_thread_num();
        thread_agg[tid][si].sum_exact += ej;
        thread_agg[tid][si].sum_est   += est;
        thread_agg[tid][si].sum_ae    += ae;
        thread_agg[tid][si].sum_se    += se;
        thread_agg[tid][si].n         += 1;
    }

    std::vector<Agg> agg(NSIZES);
    for (int t = 0; t < threads; ++t)
        for (int si = 0; si < NSIZES; ++si) {
            agg[si].sum_exact += thread_agg[t][si].sum_exact;
            agg[si].sum_est   += thread_agg[t][si].sum_est;
            agg[si].sum_ae    += thread_agg[t][si].sum_ae;
            agg[si].sum_se    += thread_agg[t][si].sum_se;
            agg[si].n         += thread_agg[t][si].n;
        }

    for (int si = 0; si < NSIZES; ++si) {
        if (agg[si].n <= 0) continue;
        const uint32_t k = SIZE_GRID[si];
        const double n = (double)agg[si].n;
        printf(
            "ALISIM_SWEEP dir=%s model=%s branch_len=%.6f algo=fastkmv "
            "kmer=%d size=%u bits=%d ss_a=20.0 ss_base=2.00 bd_bbits=16 bd_s64=16 drlevel=3 "
            "n=%d mean_exact_J=%.6f mean_est_J=%.6f MAE=%.6f RMSE=%.6f\n",
            dir.c_str(), model.c_str(), branch_len, KMER, k, PLACEHOLDER_BITS,
            agg[si].n, agg[si].sum_exact / n, agg[si].sum_est / n,
            agg[si].sum_ae / n, std::sqrt(agg[si].sum_se / n));
    }

    return 0;
}
