/**
 * eval_param_sweep – per-algorithm parameter sensitivity sweep.
 *
 * Evaluates ONE sketch algorithm over a dense grid of its key parameters.
 * Three evaluation modes:
 *
 *   jaccard     (default) – standard Jaccard accuracy on equal-length random DNA.
 *   containment           – B is a shorter prefix of A; measure containment C(B⊆A).
 *   lowcomp               – A/B have a low-complexity prefix contaminating the sketch.
 *
 * Usage:
 *   ./exe_eval_param_sweep --algo NAME [options]
 *
 * Algorithms:  minhash | fastkmv | hll | setsketch | bindash | kssd | probmh
 *
 * Common options:
 *   --mode   MODE    jaccard|containment|lowcomp     (default: jaccard)
 *   --pairs  N       pairs per mutation rate         (default: 300)
 *   --len    L       length of sequence A in bp      (default: 1000000)
 *   --len-b  LB      length of sequence B (containment mode, default: same as --len)
 *   --threads T      OpenMP threads                  (default: all available)
 *   --kmer   K       k-mer size                      (default: 21; HLL/SS: 32)
 *   --lowcomp-frac F fraction of low-complexity bases in lowcomp mode (default: 0.5)
 *
 * Algorithm-specific options:
 *   MinHash/FastKMV/PMH:
 *     --size  S      sketch size (s for MinHash/FKMV, m for PMH)  (default: 1024)
 *   HLL/SetSketch:
 *     --bits  B      log2(registers) i.e. np              (default: 13)
 *   SetSketch only:
 *     --ss-a  A      a parameter                           (default: 5.0)
 *     --ss-base X    base parameter                        (default: 2.0)
 *   BinDash:
 *     --bd-bbits BB  b-bits per bin                        (default: 16)
 *     --bd-s64 S     sketch64 (bins = S×64)                (default: 16)
 *   KSSD:
 *     --drlevel D    dimensionality reduction level        (default: 3)
 *
 * Output (stderr, machine-readable):
 *   SWEEP mode=M algo=NAME kmer=K size=S bits=B ... len=L len_b=LB rate=R n=N MAE=M RMSE=R
 */

#include "Sketch.h"
#include "probmh.h"
#include "fastkmv.h"
#include "BinDash.h"
#include "SetSketch.h"
#include "common.h"

#include <omp.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <random>
#include <string>
#include <vector>

// ─── sequence utilities ────────────────────────────────────────────────────────

static const char BASES[] = "ACGT";

static void gen_random_seq(char* buf, int len, std::mt19937_64& rng)
{
    std::uniform_int_distribution<int> d4(0, 3);
    for (int i = 0; i < len; i++) buf[i] = BASES[d4(rng)];
    buf[len] = '\0';
}

// Low-complexity biased sequence: A_frac fraction 'A', rest split equally C/G/T.
static void gen_biased_seq(char* buf, int len, std::mt19937_64& rng, double A_frac = 0.75)
{
    static const char OTHER[] = "CGT";
    std::uniform_real_distribution<double> coin(0.0, 1.0);
    std::uniform_int_distribution<int>     d3(0, 2);
    for (int i = 0; i < len; i++)
        buf[i] = (coin(rng) < A_frac) ? 'A' : OTHER[d3(rng)];
    buf[len] = '\0';
}

static void gen_mutant(const char* base, char* out, int len,
                       double rate, std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> coin(0.0, 1.0);
    std::uniform_int_distribution<int>     d3(0, 2);
    for (int i = 0; i < len; i++) {
        if (coin(rng) < rate) {
            char alt[3]; int idx = 0;
            for (int b = 0; b < 4; b++)
                if (BASES[b] != base[i]) alt[idx++] = BASES[b];
            out[i] = alt[d3(rng)];
        } else {
            out[i] = base[i];
        }
    }
    out[len] = '\0';
}

static double theo_jaccard(double p, int k)
{
    double q = std::pow(1.0 - p, k);
    return q / (2.0 - q);
}

// True containment C(B⊆A): fraction of B's k-mers unchanged → (1-p)^k
static double theo_containment(double p, int k)
{
    return std::pow(1.0 - p, k);
}

// ─── KSSD shared-param cache (copy from sub_command.cpp pattern) ───────────────

static std::shared_ptr<const Sketch::kssd_parameter_t>
get_kssd_params(int half_k, int half_subk, int drlevel) {
    try { return Sketch::getKssdParameters(half_k, half_subk, drlevel); }
    catch (const std::invalid_argument&) { return nullptr; }
}

// ─── per-algo evaluation ──────────────────────────────────────────────────────

struct EvalResult { double mae, rmse; int n; };

enum class EvalMode { JACCARD, CONTAINMENT, LOWCOMP };

// Evaluate one algorithm.
// mode=CONTAINMENT : seq_b_all[i] has length len_b (shorter than A's len_a).
//   True value = theo_containment(rate, k).
// mode=LOWCOMP     : sequences have a low-complexity prefix of length lowcomp_len.
//   True value = theo_jaccard(rate, k) of the CLEAN region only.
// mode=JACCARD     : standard.
static std::vector<EvalResult>
eval_algo(EvalMode mode,
          const std::string& algo,
          const std::vector<std::vector<char>>& seq_a_all,
          const std::vector<std::vector<char>>& seq_b_all,
          const std::vector<double>& rates,
          int pairs_per_rate,
          int kmer, int size, int bits,
          double ss_a, double ss_base,
          int bd_bbits, int bd_s64,
          int drlevel,
          int threads)
{
    int n_rates     = (int)rates.size();
    int total_pairs = n_rates * pairs_per_rate;

    std::shared_ptr<const Sketch::kssd_parameter_t> kssd_params;
    int kssd_k = kmer;
    if (algo == "kssd") {
        kssd_params = get_kssd_params(10, 6, drlevel);
        if (kssd_params) kssd_k = kssd_params->half_k * 2;
    }

    std::vector<std::vector<double>> tl_ae(threads, std::vector<double>(n_rates, 0.0));
    std::vector<std::vector<double>> tl_se(threads, std::vector<double>(n_rates, 0.0));
    std::vector<std::vector<int>>    tl_cnt(threads, std::vector<int>(n_rates, 0));

    #pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int idx = 0; idx < total_pairs; idx++) {
        int    ri    = idx / pairs_per_rate;
        double rate  = rates[ri];
        int    tid   = omp_get_thread_num();

        const char* sa = seq_a_all[idx].data();
        const char* sb = seq_b_all[idx].data();
        int La = (int)seq_a_all[idx].size() - 1;
        int Lb = (int)seq_b_all[idx].size() - 1;

        double est = 0.0;
        int    gt_k = (algo == "kssd") ? kssd_k : kmer;

        // ── build sketches ─────────────────────────────────────────────────
        if (algo == "minhash") {
            Sketch::MinHash m1(kmer, size, 42, true), m2(kmer, size, 42, true);
            m1.update(const_cast<char*>(sa));
            m2.update(const_cast<char*>(sb));
            m1.finalize(); m2.finalize();
            est = (mode == EvalMode::CONTAINMENT)
                  ? m1.containJaccard(&m2)   // |A∩B|/|A| via sorted intersection
                  : m1.jaccard(&m2);
        }
        else if (algo == "fastkmv") {
            Sketch::FastKMV v1((uint32_t)size, kmer, 42), v2((uint32_t)size, kmer, 42);
            v1.update(const_cast<char*>(sa), La);
            v2.update(const_cast<char*>(sb), Lb);
            est = (mode == EvalMode::CONTAINMENT)
                  ? v2.containment(v1)   // C(B⊆A): B.containment(A)
                  : v1.jaccard(v2);
        }
        else if (algo == "hll") {
            Sketch::HyperLogLog h1(bits, kmer), h2(bits, kmer);
            h1.update(const_cast<char*>(sa));
            h2.update(const_cast<char*>(sb));
            est = (mode == EvalMode::CONTAINMENT)
                  ? h2.containment(h1)
                  : h1.jaccard_index(h2);
        }
        else if (algo == "setsketch") {
            Sketch::SetSketch s1(bits, ss_base, ss_a, kmer);
            Sketch::SetSketch s2(bits, ss_base, ss_a, kmer);
            s1.update(const_cast<char*>(sa));
            s2.update(const_cast<char*>(sb));
            est = (mode == EvalMode::CONTAINMENT)
                  ? s2.containment(s1)
                  : s1.jaccard_index(s2);
        }
        else if (algo == "bindash") {
            Sketch::BinDash b1((uint32_t)bd_s64, kmer, (uint32_t)bd_bbits, 42);
            Sketch::BinDash b2((uint32_t)bd_s64, kmer, (uint32_t)bd_bbits, 42);
            b1.update(const_cast<char*>(sa), (size_t)La);
            b2.update(const_cast<char*>(sb), (size_t)Lb);
            est = (mode == EvalMode::CONTAINMENT)
                  ? b2.containment(b1)
                  : b1.jaccard(b2);
        }
        else if (algo == "kssd") {
            // KSSD has no containment method – always use jaccard
            Sketch::Kssd* k1 = new Sketch::Kssd(kssd_params);
            Sketch::Kssd* k2 = new Sketch::Kssd(kssd_params);
            k1->update(sa, static_cast<size_t>(La));
            k2->update(sb, static_cast<size_t>(Lb));
            est = k1->jaccard(k2);
            delete k1; delete k2;
        }
        else if (algo == "probmh") {
            Sketch::ProbMinHash4 p1((uint32_t)size, kmer, 42), p2((uint32_t)size, kmer, 42);
            p1.update(const_cast<char*>(sa), La);
            p2.update(const_cast<char*>(sb), Lb);
            est = (mode == EvalMode::CONTAINMENT)
                  ? p2.containment(p1)
                  : p1.jaccard(p2);
        }

        // ── true value ─────────────────────────────────────────────────────
        double true_val = (mode == EvalMode::CONTAINMENT)
                          ? theo_containment(rate, gt_k)
                          : theo_jaccard(rate, gt_k);

        double err = est - true_val;
        tl_ae[tid][ri]  += std::fabs(err);
        tl_se[tid][ri]  += err * err;
        tl_cnt[tid][ri] += 1;
    }

    std::vector<EvalResult> res(n_rates);
    for (int ri = 0; ri < n_rates; ri++) {
        double ae = 0, se = 0; int n = 0;
        for (int t = 0; t < threads; t++) {
            ae += tl_ae[t][ri];
            se += tl_se[t][ri];
            n  += tl_cnt[t][ri];
        }
        res[ri].mae  = n > 0 ? ae / n : 0.0;
        res[ri].rmse = n > 0 ? std::sqrt(se / n) : 0.0;
        res[ri].n    = n;
    }
    return res;
}

// ─── CLI helpers ──────────────────────────────────────────────────────────────

static const char* get_arg(int argc, char** argv, const char* flag, const char* def)
{
    for (int i = 1; i < argc - 1; i++)
        if (!strcmp(argv[i], flag)) return argv[i + 1];
    return def;
}
static int    geti(int argc, char** argv, const char* flag, int def)    { auto s=get_arg(argc,argv,flag,nullptr); return s?atoi(s):def; }
static double getd(int argc, char** argv, const char* flag, double def) { auto s=get_arg(argc,argv,flag,nullptr); return s?atof(s):def; }
static std::string gets(int argc, char** argv, const char* flag, const char* def) {
    auto s = get_arg(argc, argv, flag, def); return s ? std::string(s) : std::string(def ? def : "");
}

// ─── main ─────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[])
{
    if (argc < 2 || !strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) {
        fprintf(stderr,
            "Usage: exe_eval_param_sweep --algo NAME [options]\n"
            "  --mode   MODE    jaccard|containment|lowcomp  (default: jaccard)\n"
            "  --algo   NAME    minhash|fastkmv|hll|setsketch|bindash|kssd|probmh\n"
            "  --pairs  N       pairs per rate              (default: 300)\n"
            "  --len    L       length of seq A in bp       (default: 1000000)\n"
            "  --len-b  LB      length of seq B (containment, default = len)\n"
            "  --lowcomp-frac F low-complexity fraction     (default: 0.5)\n"
            "  --threads T      OpenMP threads              (default: all)\n"
            "  --kmer   K       k-mer size                  (default: 21)\n"
            "  --size   S       sketch size (MH/FKMV/PMH)   (default: 1024)\n"
            "  --bits   B       np bits (HLL/SS)            (default: 13)\n"
            "  --ss-a   A       SetSketch a                 (default: 5.0)\n"
            "  --ss-base X      SetSketch base              (default: 2.0)\n"
            "  --bd-bbits BB    BinDash b-bits              (default: 16)\n"
            "  --bd-s64 S       BinDash sketch64            (default: 16)\n"
            "  --drlevel D      KSSD drlevel                (default: 3)\n");
        return 1;
    }

    const std::string mode_str = gets(argc, argv, "--mode", "jaccard");
    const std::string algo     = gets(argc, argv, "--algo", "minhash");
    const int pairs      = geti(argc, argv, "--pairs",   300);
    const int seq_len_a  = geti(argc, argv, "--len",     1000000);
    const int seq_len_b_arg = geti(argc, argv, "--len-b", -1);
    const int seq_len_b  = (seq_len_b_arg > 0) ? seq_len_b_arg : seq_len_a;
    const double lowcomp_frac = getd(argc, argv, "--lowcomp-frac", 0.5);
    const int kmer       = geti(argc, argv, "--kmer",
                                (algo=="hll"||algo=="setsketch") ? 32 : 21);
    const int size       = geti(argc, argv, "--size",    1024);
    const int bits       = geti(argc, argv, "--bits",    13);
    const double ss_a    = getd(argc, argv, "--ss-a",    5.0);
    const double ss_base = getd(argc, argv, "--ss-base", 2.0);
    const int bd_bbits   = geti(argc, argv, "--bd-bbits", 16);
    const int bd_s64     = geti(argc, argv, "--bd-s64",  16);
    const int drlevel    = geti(argc, argv, "--drlevel", 3);
    int threads          = geti(argc, argv, "--threads", omp_get_max_threads());

    EvalMode mode = EvalMode::JACCARD;
    if (mode_str == "containment") mode = EvalMode::CONTAINMENT;
    else if (mode_str == "lowcomp") mode = EvalMode::LOWCOMP;

    const std::vector<std::string> valid_algos =
        {"minhash","fastkmv","hll","setsketch","bindash","kssd","probmh"};
    if (std::find(valid_algos.begin(), valid_algos.end(), algo) == valid_algos.end()) {
        fprintf(stderr, "Unknown algo '%s'\n", algo.c_str());
        return 1;
    }
    if (threads < 1) threads = 1;
    if (pairs < 1)   { fprintf(stderr, "pairs must be >= 1\n"); return 1; }

    // Mutation rates: from near-identical to highly divergent
    const std::vector<double> rates = {
        0.001, 0.01, 0.05, 0.10, 0.20, 0.40, 0.80
    };
    const int n_rates     = (int)rates.size();
    const int total_pairs = n_rates * pairs;

    // Low-complexity region length
    const int lc_len = (mode == EvalMode::LOWCOMP)
                       ? (int)(seq_len_a * lowcomp_frac)
                       : 0;

    fprintf(stderr,
        "# eval_param_sweep mode=%s algo=%s kmer=%d size=%d bits=%d "
        "ss_a=%.1f ss_base=%.1f bd_bbits=%d bd_s64=%d drlevel=%d "
        "len=%d len_b=%d lowcomp_frac=%.2f pairs=%d threads=%d\n",
        mode_str.c_str(), algo.c_str(), kmer, size, bits,
        ss_a, ss_base, bd_bbits, bd_s64, drlevel,
        seq_len_a, seq_len_b, lowcomp_frac, pairs, threads);

    double t0 = get_sec();

    // ── Pre-generate all sequence pairs ────────────────────────────────────
    std::vector<std::vector<char>> seq_a_all(total_pairs), seq_b_all(total_pairs);

    #pragma omp parallel for schedule(static) num_threads(threads)
    for (int idx = 0; idx < total_pairs; idx++) {
        int ri = idx / pairs;
        std::mt19937_64 rng(12345ULL + (uint64_t)idx * 999983ULL);

        seq_a_all[idx].resize(seq_len_a + 1);
        seq_b_all[idx].resize(seq_len_b + 1);

        if (mode == EvalMode::LOWCOMP) {
            // Biased prefix (identical in both), then clean region mutated
            if (lc_len > 0) {
                gen_biased_seq(seq_a_all[idx].data(), lc_len, rng);
                memcpy(seq_b_all[idx].data(), seq_a_all[idx].data(), (size_t)lc_len);
            }
            int clean = seq_len_a - lc_len;
            if (clean > 0) {
                gen_random_seq(seq_a_all[idx].data() + lc_len, clean, rng);
                gen_mutant(seq_a_all[idx].data() + lc_len,
                           seq_b_all[idx].data() + lc_len,
                           clean, rates[ri], rng);
            }
            seq_a_all[idx][seq_len_a] = '\0';
            seq_b_all[idx][seq_len_b] = '\0';
        }
        else if (mode == EvalMode::CONTAINMENT) {
            // A: full random, B: first seq_len_b bases of A mutated
            gen_random_seq(seq_a_all[idx].data(), seq_len_a, rng);
            gen_mutant(seq_a_all[idx].data(), seq_b_all[idx].data(),
                       seq_len_b, rates[ri], rng);
            seq_b_all[idx][seq_len_b] = '\0';
        }
        else {
            // Standard Jaccard: equal-length random pair
            gen_random_seq(seq_a_all[idx].data(), seq_len_a, rng);
            gen_mutant(seq_a_all[idx].data(), seq_b_all[idx].data(),
                       seq_len_a, rates[ri], rng);
        }
    }

    double t1 = get_sec();
    fprintf(stderr, "# seqgen: %.1f s\n", t1 - t0);

    // ── Evaluate ───────────────────────────────────────────────────────────
    auto res = eval_algo(mode, algo, seq_a_all, seq_b_all,
                         rates, pairs,
                         kmer, size, bits, ss_a, ss_base,
                         bd_bbits, bd_s64, drlevel, threads);

    double t2 = get_sec();
    fprintf(stderr, "# eval:   %.1f s\n", t2 - t1);

    // ── Output SWEEP lines ─────────────────────────────────────────────────
    for (int ri = 0; ri < n_rates; ri++) {
        fprintf(stderr,
            "SWEEP mode=%s algo=%s kmer=%d size=%d bits=%d ss_a=%.1f ss_base=%.1f "
            "bd_bbits=%d bd_s64=%d drlevel=%d len=%d len_b=%d "
            "rate=%.4f n=%d MAE=%.6f RMSE=%.6f\n",
            mode_str.c_str(), algo.c_str(),
            kmer, size, bits, ss_a, ss_base,
            bd_bbits, bd_s64, drlevel,
            seq_len_a, seq_len_b,
            rates[ri], res[ri].n, res[ri].mae, res[ri].rmse);
    }

    double avg_mae = 0; int total_n = 0;
    for (auto& r : res) { avg_mae += r.mae * r.n; total_n += r.n; }
    if (total_n > 0) avg_mae /= total_n;
    fprintf(stderr,
        "SWEEP_SUMMARY mode=%s algo=%s kmer=%d size=%d bits=%d ss_a=%.1f ss_base=%.1f "
        "bd_bbits=%d bd_s64=%d drlevel=%d len=%d len_b=%d avg_MAE=%.6f\n",
        mode_str.c_str(), algo.c_str(),
        kmer, size, bits, ss_a, ss_base,
        bd_bbits, bd_s64, drlevel, seq_len_a, seq_len_b, avg_mae);

    fprintf(stderr, "# total:  %.1f s\n", get_sec() - t0);
    return 0;
}
