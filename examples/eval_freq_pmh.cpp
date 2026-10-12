/**
 * eval_freq_pmh.cpp
 * Evaluate frequency-weighted ProbMinHash4 (PMH-freq) on bench_data .fa files.
 *
 * For each pair{id}_A.fa / pair{id}_B.fa:
 *   1. Count k-mer frequencies via sort (1/8 hash-sampling for memory)
 *   2. Build ProbMinHash4 sketch with weight = count/total (PMH-freq)
 *   3. Estimate Jaccard via jaccard_weighted()
 *   4. Compare against exact_j_k21 ground truth
 *
 * Usage:
 *   exe_eval_freq_pmh <bench_dir> <threads>
 *
 * Output (stderr):  PARSE_RATE len=<bp> rate=<r> FreqPMH=<mae>
 */

#include "probmh.h"
#include "kseq.h"
#include <zlib.h>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <atomic>
#include <sys/time.h>

static double get_sec() {
    struct timeval tv; gettimeofday(&tv, nullptr);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

KSEQ_INIT(gzFile, gzread)

// ── Exact copy of constants from sub_command.cpp ────────────────────────────
static constexpr int      KSIZE      = 21;
static constexpr uint32_t M          = 1024;
static constexpr uint64_t SEED       = 42;
static constexpr uint64_t PMH_SDIV   = 8;
static constexpr uint64_t PMH_STHRESH = UINT64_MAX / PMH_SDIV;
static constexpr uint64_t PMH_SMIX   = 0x9E3779B97F4A7C15ULL;

static const uint8_t PMH_BASE2[256] = {
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,   0,0xFF,   1,0xFF,0xFF,0xFF,   2,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,   3,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,   0,0xFF,   1,0xFF,0xFF,0xFF,   2,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,   3,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
};

// Exact pmh_kmers logic from sub_command.cpp
static inline void pmh_kmers(const char* seq, int len,
                              const uint64_t mask, const uint64_t hi_shift,
                              std::vector<uint64_t>& sampled, uint64_t& total) {
    uint64_t fwd = 0, rev = 0; int valid = 0;
    for (int i = 0; i < len; ++i) {
        const uint8_t b = PMH_BASE2[(uint8_t)seq[i]];
        if (b == 0xFF) { valid = 0; fwd = 0; rev = 0; continue; }
        fwd = ((fwd << 2) | b) & mask;
        rev = (rev >> 2) | ((uint64_t)(3 ^ b) << hi_shift);
        if (++valid >= KSIZE) {
            const uint64_t can = (fwd <= rev) ? fwd : rev;
            ++total;
            if ((can * PMH_SMIX) <= PMH_STHRESH)
                sampled.push_back(can);
        }
    }
}

// ── Build frequency-weighted PMH sketch from a FASTA(.gz) file ──────────────
// Mirrors pmh_build_weighted in sub_command.cpp but WITHOUT singleton filter
// (assembly sequences: every k-mer appears exactly once, no true singletons).
static Sketch::ProbMinHash4* build_freq_pmh(const std::string& path) {
    const uint64_t mask     = (1ULL << (2 * KSIZE)) - 1;
    const uint64_t hi_shift = 2 * (KSIZE - 1);

    std::vector<uint64_t> sampled;
    sampled.reserve(1 << 20);
    uint64_t total = 0;

    gzFile fp = gzopen(path.c_str(), "r");
    if (!fp) return nullptr;
    kseq_t* ks = kseq_init(fp);
    while (kseq_read(ks) >= 0) {
        if ((int)ks->seq.l >= KSIZE)
            pmh_kmers(ks->seq.s, (int)ks->seq.l, mask, hi_shift, sampled, total);
    }
    kseq_destroy(ks);
    gzclose(fp);

    std::sort(sampled.begin(), sampled.end());

    auto* sk = new Sketch::ProbMinHash4(M, KSIZE, SEED, 0);
    if (total > 0 && !sampled.empty()) {
        const double inv = 1.0 / (double)total;
        size_t i = 0;
        while (i < sampled.size()) {
            size_t j = i + 1;
            while (j < sampled.size() && sampled[j] == sampled[i]) ++j;
            const uint32_t cnt = (uint32_t)(j - i);
            // No singleton filter: assembly has no error k-mers, every unique
            // k-mer is genuine regardless of count
            sk->addHash(sampled[i], (double)cnt * inv);
            i = j;
        }
    }
    return sk;
}

// ── Main ─────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <bench_dir> [threads]\n", argv[0]);
        return 1;
    }
    std::string bench_dir = argv[1];
    int nthreads = (argc >= 3) ? atoi(argv[2]) : 8;

    // ── Load pairs_meta.tsv ─────────────────────────────────────────────────
    struct PairMeta { int id; double rate; int len_bp; double gt_j21; };
    std::vector<PairMeta> pairs;
    {
        std::ifstream fin(bench_dir + "/pairs_meta.tsv");
        if (!fin) { fprintf(stderr, "ERROR: cannot open pairs_meta.tsv\n"); return 1; }
        std::string line;
        std::getline(fin, line); // header
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::string tok;
            std::vector<std::string> cols;
            while (std::getline(ss, tok, '\t')) cols.push_back(tok);
            if (cols.size() < 4) continue;
            pairs.push_back({std::stoi(cols[0]), std::stod(cols[1]),
                             std::stoi(cols[2]), std::stod(cols[3])});
        }
    }

    // Override gt_j21 with exact values if available
    {
        std::ifstream fe(bench_dir + "/pairs_exact_j.tsv");
        if (fe) {
            std::string line; std::getline(fe, line);
            while (std::getline(fe, line)) {
                if (line.empty()) continue;
                std::istringstream ss(line);
                std::string tok;
                std::vector<std::string> cols;
                while (std::getline(ss, tok, '\t')) cols.push_back(tok);
                if (cols.size() < 4) continue;
                int pid = std::stoi(cols[0]);
                double ej = std::stod(cols[3]);
                for (auto& p : pairs) if (p.id == pid) { p.gt_j21 = ej; break; }
            }
        }
    }

    int N = (int)pairs.size();
    fprintf(stderr, "=== FreqPMH eval: %s  (%d pairs, k=%d, m=%d, threads=%d) ===\n",
            bench_dir.c_str(), N, KSIZE, M, nthreads);

    // ── Collect distinct rates for per-rate reporting ────────────────────────
    std::vector<double> rates;
    for (auto& p : pairs) {
        if (std::find(rates.begin(), rates.end(), p.rate) == rates.end())
            rates.push_back(p.rate);
    }
    std::sort(rates.begin(), rates.end());
    int n_rates = (int)rates.size();

    // ── Per-pair evaluation ──────────────────────────────────────────────────
    struct Res { double est_j; double gt_j; double rate; int len_bp; };
    std::vector<Res> results(N);
    std::atomic<int> errors(0);

    double t0 = get_sec();
    #pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
    for (int i = 0; i < N; i++) {
        const PairMeta& pm = pairs[i];
        char id_buf[16];
        std::snprintf(id_buf, sizeof(id_buf), "%06d", pm.id);
        std::string fa = bench_dir + "/pair" + id_buf + "_A.fa";
        std::string fb = bench_dir + "/pair" + id_buf + "_B.fa";

        Sketch::ProbMinHash4* sa = build_freq_pmh(fa);
        Sketch::ProbMinHash4* sb = build_freq_pmh(fb);

        double j = 0.0;
        if (sa && sb) j = sa->jaccard_weighted(*sb);
        else ++errors;

        delete sa; delete sb;
        results[i] = {j, pm.gt_j21, pm.rate, pm.len_bp};
    }
    double t1 = get_sec();

    if (errors > 0)
        fprintf(stderr, "WARN: %d I/O errors\n", (int)errors);

    // ── Infer representative length ──────────────────────────────────────────
    int rep_len = 0;
    if (!results.empty()) rep_len = results[0].len_bp;
    if (rep_len == 0) {
        // parse from dir name (L1M, L100k, etc.)
        std::string dname = bench_dir;
        auto sl = dname.rfind('/');
        if (sl != std::string::npos) dname = dname.substr(sl + 1);
        if (!dname.empty() && (dname[0]=='L'||dname[0]=='l')) {
            std::string num = dname.substr(1);
            int mult = 1;
            if (!num.empty()&&(num.back()=='k'||num.back()=='K')){mult=1000;num.pop_back();}
            else if (!num.empty()&&(num.back()=='m'||num.back()=='M')){mult=1000000;num.pop_back();}
            try { rep_len = std::stoi(num)*mult; } catch(...) {}
        }
    }

    // ── Per-rate accumulation ────────────────────────────────────────────────
    fprintf(stderr,
        "\n  rate      FreqPMH_MAE  FreqPMH_RMSE\n"
        "  ─────────────────────────────────────\n");

    double g_err = 0; double g_err2 = 0; int g_n = 0;
    for (int ri = 0; ri < n_rates; ri++) {
        double r = rates[ri];
        double sum_e = 0, sum_e2 = 0; int cnt = 0;
        for (auto& res : results) {
            if (res.rate != r) continue;
            double e = std::fabs(res.est_j - res.gt_j);
            sum_e += e; sum_e2 += e*e; ++cnt;
        }
        if (cnt == 0) continue;
        double mae  = sum_e  / cnt;
        double rmse = std::sqrt(sum_e2 / cnt);
        fprintf(stderr, "  %.4f    %11.6f  %12.6f\n", r, mae, rmse);
        fprintf(stderr,
            "PARSE_RATE len=%d rate=%.4f FreqPMH=%.6f FreqPMH_RMSE=%.6f\n",
            rep_len, r, mae, rmse);
        g_err += sum_e; g_err2 += sum_e2; g_n += cnt;
    }
    double g_mae = g_n ? g_err/g_n : 0;
    fprintf(stderr,
        "  ─────────────────────────────────────\n"
        "  AVG      %11.6f\n"
        "  time: %.1f s\n",
        g_mae, t1 - t0);

    return 0;
}
