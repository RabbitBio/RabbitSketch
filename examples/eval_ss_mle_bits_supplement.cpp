/**
 * SetSketch MLE bits sweep at default a=5, base=2, k=21 on existing AliSim pair FASTAs.
 *
 * Outputs ALISIM_SWEEP lines (algo=setsketch_mle) for bits=8..16 only.
 * Append to results/alisim_sweep.csv, then: python3 build_master_csv.py
 */

#include "SetSketch.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <omp.h>

namespace {

constexpr int KMER = 21;
constexpr double SS_A = 5.0;
constexpr double SS_BASE = 2.0;
const int BITS_LEVELS[] = {8, 9, 10, 11, 12, 13, 14, 15, 16};

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

double exact_jaccard_k21(const std::string& a, const std::string& b) {
    const uint64_t mask = (1ULL << (2 * KMER)) - 1;
    const int hi = 2 * (KMER - 1);
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
            if (++valid >= KMER) v.push_back(fwd <= rev ? fwd : rev);
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

    omp_set_num_threads(threads);
    const int nbits = (int)(sizeof(BITS_LEVELS) / sizeof(BITS_LEVELS[0]));
    std::vector<Agg> agg(nbits);

#pragma omp parallel for schedule(dynamic, 4)
    for (size_t pi = 0; pi < pairs.size(); ++pi) {
        std::string sa = read_fasta(pairs[pi].fa.c_str());
        std::string sb = read_fasta(pairs[pi].fb.c_str());
        if (sa.empty() || sb.empty()) continue;

        const double ej = exact_jaccard_k21(sa, sb);
        char* ca = sa.data();
        char* cb = sb.data();

        for (int bi = 0; bi < nbits; ++bi) {
            const int bits = BITS_LEVELS[bi];
            Sketch::SetSketch sk1(bits, SS_BASE, SS_A, KMER);
            Sketch::SetSketch sk2(bits, SS_BASE, SS_A, KMER);
            sk1.update(ca, sa.size());
            sk2.update(cb, sb.size());
            const double est = sk1.jaccard_index_mle(sk2);
            const double ae = std::fabs(est - ej);
            const double se = (est - ej) * (est - ej);

#pragma omp critical
            {
                agg[bi].sum_exact += ej;
                agg[bi].sum_est += est;
                agg[bi].sum_ae += ae;
                agg[bi].sum_se += se;
                agg[bi].n += 1;
            }
        }
    }

    if (branch_len < 0 && !pairs.empty()) {
        // Read BL from pairs_meta if present
        char meta[4096];
        snprintf(meta, sizeof(meta), "%s/pairs_meta.tsv", dir.c_str());
        FILE* mf = fopen(meta, "r");
        if (mf) {
            char line[256];
            fgets(line, sizeof(line), mf);  // header
            if (fgets(line, sizeof(line), mf)) {
                double bl = 0;
                sscanf(line, "%*d%lf", &bl);
                if (bl > 0) branch_len = bl;
            }
            fclose(mf);
        }
    }
    if (branch_len < 0) branch_len = 0.0;

    for (int bi = 0; bi < nbits; ++bi) {
        if (agg[bi].n <= 0) continue;
        const int bits = BITS_LEVELS[bi];
        const double n = (double)agg[bi].n;
        const double mej = agg[bi].sum_exact / n;
        const double mest = agg[bi].sum_est / n;
        const double mae = agg[bi].sum_ae / n;
        const double rmse = std::sqrt(agg[bi].sum_se / n);
        printf(
            "ALISIM_SWEEP dir=%s model=%s branch_len=%.6f algo=setsketch_mle "
            "kmer=%d size=0 bits=%d ss_a=%.1f ss_base=%.1f bd_bbits=0 bd_s64=0 drlevel=0 "
            "n=%d mean_exact_J=%.6f mean_est_J=%.6f MAE=%.6f RMSE=%.6f\n",
            dir.c_str(), model.c_str(), branch_len, KMER, bits, SS_A, SS_BASE,
            agg[bi].n, mej, mest, mae, rmse);
    }

    return 0;
}
