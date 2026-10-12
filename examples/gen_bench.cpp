/*
 * gen_bench.cpp  –  Generate variable-length simulated sequence pairs
 *
 * For each (length, mutation_rate) combination, generates N_PAIRS pairs:
 *   - Random reference sequence of given length
 *   - Independently mutated A and B sequences
 *   - Computes exact Binary Jaccard (k=21) via canonical k-mer hash sets
 *   - Writes FASTA files + pairs_meta.tsv + pairs_exact_j.tsv
 *
 * Usage:  ./gen_bench <out_dir> [seed]
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <omp.h>

// ── Config ──────────────────────────────────────────────────────────────────
static const int K           = 21;
static const int N_PAIRS     = 500;       // pairs per (len, rate) combo
static const int WRAP        = 80;        // FASTA line width

static const std::vector<int>    LENGTHS = {100000, 500000, 1000000,
                                            2000000, 4000000, 8000000};
static const std::vector<double> RATES   = {0.001, 0.010, 0.050, 0.100, 0.200};

// ── Utilities ────────────────────────────────────────────────────────────────
static const char BASES[4] = {'A','C','G','T'};
static const int  COMP[4]  = { 3, 2, 1, 0 };   // complement index

static void write_fasta(const char* path, const char* name, const std::string& seq) {
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); exit(1); }
    fprintf(f, ">%s\n", name);
    int n = (int)seq.size();
    for (int i = 0; i < n; i += WRAP) {
        int end = std::min(i + WRAP, n);
        fwrite(seq.data() + i, 1, end - i, f);
        fputc('\n', f);
    }
    fclose(f);
}


// Sort-based Jaccard: faster and more cache-friendly than hash sets for large k
static double exact_jaccard(const std::string& a, const std::string& b) {
    // Collect k-mers into sorted vectors, then dedup and merge-count
    auto collect_sorted = [&](const std::string& seq) {
        std::vector<uint64_t> v;
        v.reserve(seq.size());
        const uint64_t mask     = (1ULL << (2 * K)) - 1;
        const int      hi_shift = 2 * (K - 1);
        uint64_t fwd = 0, rev = 0; int valid = 0;
        for (unsigned char ch : seq) {
            int b;
            switch (ch) {
                case 'A': case 'a': b=0; break;
                case 'C': case 'c': b=1; break;
                case 'G': case 'g': b=2; break;
                case 'T': case 't': b=3; break;
                default: valid=0; fwd=rev=0; continue;
            }
            fwd = ((fwd<<2)|b)&mask;
            rev = (rev>>2)|((uint64_t)COMP[b]<<hi_shift);
            if (++valid>=K) v.push_back(fwd<=rev?fwd:rev);
        }
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        return v;
    };

    auto va = collect_sorted(a);
    auto vb = collect_sorted(b);

    size_t i=0, j=0, inter=0;
    while (i<va.size() && j<vb.size()) {
        if      (va[i]<vb[j]) ++i;
        else if (va[i]>vb[j]) ++j;
        else { ++inter; ++i; ++j; }
    }
    size_t uni = va.size() + vb.size() - inter;
    return uni ? (double)inter/uni : 0.0;
}

// Mutate: each base independently mutated with probability `rate`
static std::string mutate(const std::string& seq, double rate, std::mt19937_64& rng) {
    std::bernoulli_distribution  mut_d(rate);
    std::uniform_int_distribution<int> base_d(0, 2);   // 3 choices ≠ original
    std::string out = seq;
    for (char& c : out) {
        if (mut_d(rng)) {
            int orig;
            switch (c) {
                case 'A': orig=0; break; case 'C': orig=1; break;
                case 'G': orig=2; break; default: orig=3; break;
            }
            c = BASES[(orig + 1 + base_d(rng)) % 4];
        }
    }
    return out;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <out_dir> [seed]\n", argv[0]); return 1; }

    const std::string out_dir(argv[1]);
    uint64_t seed = (argc >= 3) ? (uint64_t)atoll(argv[2]) : 12345ULL;

    // Create output directory
    mkdir(out_dir.c_str(), 0755);

    FILE* fmeta  = fopen((out_dir + "/pairs_meta.tsv").c_str(), "w");
    FILE* fexact = fopen((out_dir + "/pairs_exact_j.tsv").c_str(), "w");
    fprintf(fmeta,  "pair_id\trate\tlen_bp\ttheo_j_k21\n");
    fprintf(fexact, "pair_id\trate\tlen_bp\texact_j_k21\n");

    // File lists
    FILE* fall  = fopen((out_dir + "/all_list.txt").c_str(), "w");
    FILE* falist = fopen((out_dir + "/A_list.txt").c_str(), "w");
    FILE* fblist = fopen((out_dir + "/B_list.txt").c_str(), "w");

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> base_rand(0, 3);

    int pair_id = 0;
    int total_combos = (int)(LENGTHS.size() * RATES.size());
    int combo = 0;

    for (int L : LENGTHS) {
        // Generate one reference sequence per length
        std::string ref(L, 'A');
        for (char& c : ref) c = BASES[base_rand(rng)];

        for (double rate : RATES) {
            ++combo;
            double theo_j = std::pow(1.0 - rate, K) / (2.0 - std::pow(1.0 - rate, K));
            fprintf(stderr, "[%d/%d] len=%d rate=%.3f pairs=%d...\n",
                    combo, total_combos, L, rate, N_PAIRS);

            // Pre-generate all seeds for reproducibility under parallelism
            std::vector<uint64_t> seeds_A(N_PAIRS), seeds_B(N_PAIRS);
            for (int p = 0; p < N_PAIRS; p++) {
                seeds_A[p] = rng(); seeds_B[p] = rng();
            }

            // Results buffer
            struct PairResult { std::string seqA, seqB; double ej; };
            std::vector<PairResult> results(N_PAIRS);

            #pragma omp parallel for schedule(dynamic, 1) num_threads(32)
            for (int p = 0; p < N_PAIRS; p++) {
                std::mt19937_64 rng_A(seeds_A[p]), rng_B(seeds_B[p]);
                results[p].seqA = mutate(ref, rate, rng_A);
                results[p].seqB = mutate(ref, rate, rng_B);
                results[p].ej   = exact_jaccard(results[p].seqA, results[p].seqB);
            }

            // Write sequentially in order
            for (int p = 0; p < N_PAIRS; p++, pair_id++) {
                char fa[512], fb[512], na[64], nb[64];
                snprintf(fa, sizeof(fa), "%s/pair%06d_A.fa", out_dir.c_str(), pair_id);
                snprintf(fb, sizeof(fb), "%s/pair%06d_B.fa", out_dir.c_str(), pair_id);
                snprintf(na, sizeof(na), "pair%06d_A_L%d_r%.3f", pair_id, L, rate);
                snprintf(nb, sizeof(nb), "pair%06d_B_L%d_r%.3f", pair_id, L, rate);

                write_fasta(fa, na, results[p].seqA);
                write_fasta(fb, nb, results[p].seqB);

                fprintf(fmeta,  "%d\t%.4f\t%d\t%.6f\n", pair_id, rate, L, theo_j);
                fprintf(fexact, "%d\t%.4f\t%d\t%.8f\n", pair_id, rate, L, results[p].ej);
                fprintf(fall,   "%s\n%s\n", fa, fb);
                fprintf(falist, "%s\n", fa);
                fprintf(fblist, "%s\n", fb);
            }
        }
    }

    fclose(fmeta); fclose(fexact);
    fclose(fall); fclose(falist); fclose(fblist);

    fprintf(stderr, "Done. %d pairs written to %s\n", pair_id, out_dir.c_str());
    return 0;
}
