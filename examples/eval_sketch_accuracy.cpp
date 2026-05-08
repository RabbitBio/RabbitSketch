/**
 * eval_sketch_accuracy – synthetic benchmark for sketch Jaccard accuracy.
 *
 * FastKMV（单轮 fmix / 无 fmix）对 MinHash：定长+变长见 eval_fkmv_accuracy.cpp：
 *   make -C examples eval_fkmv_all
 *   ./eval_fkmv_{single,nofmix} <pairs_per_rate> <fixed_len_bp> <threads>
 *
 * Generates pairs of random DNA sequences with controlled per-base substitution
 * rates.  For each pair, builds all sketch types and compares the estimated
 * Jaccard against the theoretical Mash formula (unweighted) or the exact
 * sliding-window weighted Jaccard (weighted PMH evaluation).
 *
 * Sketch methods & their k-mer sizes:
 *   HyperLogLog  (k=32 hardcoded, 16384 registers, np=14)
 *   SetSketch    (k=32 hardcoded, 16384 registers, np=14)
 *   KSSD         (k=20, half_k=10 half_subk=6 drlevel=3)
 *   MinHash      (k=21, sketch size 1024)
 *   ProbMinHash4 (k=21, 1024 registers)
 *   WeightedPMH  (k=21, 1024 registers, independent random weights per sequence)
 *
 * §1  Fixed-length benchmark:  all pairs use seq_length (default 4 Mbp).
 * §2  Variable-length benchmark: lengths drawn log-uniformly from
 *     [50 Kbp, 8 Mbp].
 *
 * Weighted PMH ground truth:
 *   Weights = Shannon entropy of each k-mer, normalized to [0.1, 1.0].
 *   Same k-mer → same entropy → same weight in both sequences (ProbMinHash ✓).
 *   Exact weighted Jaccard is computed via O(L) sliding-window substitution
 *   tracking—no hash map needed. Low-complexity k-mers get low weight.
 *
 * Usage:
 *   exe_eval_sketch_acc [pairs_per_rate] [seq_length] [threads] [mode]
 *
 * mode (optional):
 *   both   — 定长 + 变长（默认）
 *   fixed  — 仅定长
 *   var    — 仅变长
 *
 * 机器可读汇总行（stderr）：
 *   PARSE_FIXED len=<bp> HLL=... SetSketch=... KSSD=... MinHash=... ProbMH=... OnePerm=... KMV=... TL1=... TL2=... TL4=... WgtPMH_vs_WJ=... UnwgtPMH_vs_WJ=...
 *   PARSE_VARLEN HLL=... ... WgtPMH_vs_WJ=... UnwgtPMH_vs_WJ=...
 */

#include "Sketch.h"
#include "probmh.h"
#include "fastkmv.h"
#include "BinDash.h"
#include "common.h"

#include <omp.h>
#include <sys/time.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <random>
#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <sstream>

using namespace std;

enum EvalMode { EVAL_BOTH = 0, EVAL_FIXED_ONLY, EVAL_VAR_ONLY, EVAL_LOWCOMP, EVAL_FILE };

static EvalMode parse_mode(int argc, char** argv) {
    if (argc <= 4)
        return EVAL_BOTH;
    if (strcmp(argv[4], "fixed") == 0)
        return EVAL_FIXED_ONLY;
    if (strcmp(argv[4], "var") == 0 || strcmp(argv[4], "varlen") == 0)
        return EVAL_VAR_ONLY;
    if (strcmp(argv[4], "lowcomp") == 0)
        return EVAL_LOWCOMP;
    if (strcmp(argv[4], "file") == 0)
        return EVAL_FILE;
    return EVAL_BOTH;
}

static const char BASES[] = "ACGT";

static double theo_jaccard(double p, int k) {
    double q = pow(1.0 - p, k);
    return q / (2.0 - q);
}

static void gen_random_seq(char* buf, int len, mt19937_64& rng) {
    uniform_int_distribution<int> d4(0, 3);
    for (int i = 0; i < len; i++)
        buf[i] = BASES[d4(rng)];
    buf[len] = '\0';
}

static void gen_mutant(const char* base, char* out, int len,
                       double rate, mt19937_64& rng) {
    uniform_real_distribution<double> coin(0.0, 1.0);
    uniform_int_distribution<int>     d3(0, 2);
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

// Like gen_mutant, but also records which base positions were substituted.
// sub_at[i] = true  ⟹  position i received a substitution.
// sub_at is resized to len before filling.
static void gen_mutant_tracked(const char* base, char* out, int len,
                                double rate, mt19937_64& rng,
                                vector<bool>& sub_at)
{
    sub_at.assign(len, false);
    uniform_real_distribution<double> coin(0.0, 1.0);
    uniform_int_distribution<int>     d3(0, 2);
    for (int i = 0; i < len; i++) {
        if (coin(rng) < rate) {
            char alt[3]; int idx = 0;
            for (int b = 0; b < 4; b++)
                if (BASES[b] != base[i]) alt[idx++] = BASES[b];
            out[i] = alt[d3(rng)];
            sub_at[i] = true;
        } else {
            out[i] = base[i];
        }
    }
    out[len] = '\0';
}

// Generates a sequence with biased base composition: A_frac fraction 'A',
// the remaining (1-A_frac) split equally among C, G, T.
// This creates low-complexity k-mers (weight ≈ 0.1–0.25) while keeping
// enough sequence diversity to produce many distinct k-mers.
static void gen_biased_seq(char* buf, int len, mt19937_64& rng, double A_frac) {
    uniform_real_distribution<double> coin(0.0, 1.0);
    uniform_int_distribution<int>     d3(0, 2);
    static const char OTHER[] = "CGT";
    for (int i = 0; i < len; i++)
        buf[i] = (coin(rng) < A_frac) ? 'A' : OTHER[d3(rng)];
    buf[len] = '\0';
}

// Exact unweighted Jaccard of only the first (L_clean - k + 1) k-mers,
// using the substitution tracking array for the clean region.
// Excludes any k-mers that touch the boundary with the biased region.
static double exact_j_clean(const vector<bool>& sub_at, int L_clean, int k) {
    const int N = max(0, L_clean - k + 1);
    int subs = 0, matched = 0;
    for (int j = 0; j < k && j < (int)sub_at.size(); j++)
        subs += sub_at[j] ? 1 : 0;
    for (int i = 0; i < N; i++) {
        if (i > 0) {
            subs -= sub_at[i - 1] ? 1 : 0;
            int tail = i + k - 1;
            if (tail < (int)sub_at.size())
                subs += sub_at[tail] ? 1 : 0;
        }
        if (subs == 0) matched++;
    }
    return (N > 0) ? (double)matched / (2 * N - matched) : 0.0;
}

// Per-k-mer Shannon entropy weights: see fill_kmer_entropy_weights() in common.h

// Compute exact weighted Jaccard in O(La + Lb) time using a sliding-window
// substitution counter — no hash map needed.
//
// Assumptions (hold for random DNA with k ≥ 21):
//   • k-mers are (almost surely) unique within each sequence.
//   • The substitution model has no indels, so position i in A aligns to
//     position i in B for positions 0 .. S-1 (S = min(La, Lb)).
//   • A k-mer starting at position i in A matches the k-mer at position i in B
//     iff zero substitutions occurred in bases [i, i+k-1].
//
// sub_at : substitutions in the shared prefix [0, S-1], size must equal S.
// wA     : per-k-mer weight array for A, length = La - k + 1.
// wB     : per-k-mer weight array for B, length = Lb - k + 1.
static double exact_wj_from_subs(const vector<bool>& sub_at,
                                  const double* wA, int La,
                                  const double* wB, int Lb,
                                  int k)
{
    const int S        = (int)sub_at.size();   // = min(La, Lb)
    const int N_shared = max(0, S - k + 1);    // k-mer positions fully inside shared prefix

    // Sliding window: count substitutions in [i, i+k-1].
    int subs = 0;
    for (int j = 0; j < k && j < S; j++)
        subs += sub_at[j] ? 1 : 0;

    double num = 0.0, den = 0.0;

    for (int i = 0; i < N_shared; i++) {
        if (i > 0) {
            subs -= sub_at[i - 1] ? 1 : 0;
            subs += sub_at[i + k - 1] ? 1 : 0;
        }
        double wa = wA[i], wb = wB[i];
        if (subs == 0) {
            // k-mers match exactly
            num += min(wa, wb);
            den += max(wa, wb);
        } else {
            // k-mers differ → each contributes only to its own side of the union
            den += wa + wb;
        }
    }

    // Remaining k-mers in A beyond the shared region (only if La > Lb)
    for (int i = N_shared; i <= La - k; i++)
        den += wA[i];

    // Remaining k-mers in B beyond the shared region (only if Lb > La, or
    // boundary-spanning k-mers whose tail falls in the random extension)
    for (int i = N_shared; i <= Lb - k; i++)
        den += wB[i];

    return den > 0.0 ? num / den : 0.0;
}

// ── Result structs ─────────────────────────────────────────────────────────────

struct PairResult {
    double rate;
    double theo_j32, theo_j21, theo_j20;
    double hll_j, ss_j, kssd_j, mh_j, pmh_j, oph_j, kmv_j, bd_j;
    double tl1_j, tl2_j, tl4_j;
    // Weighted PMH evaluation (independent random weights per sequence)
    double wpmh_j;   // weighted PMH estimate
    double gt_wj;    // exact weighted Jaccard (sliding-window, O(L))
};

int main(int argc, char* argv[])
{
    int pairs_per_rate = (argc > 1) ? atoi(argv[1]) : 500;
    int seq_length     = (argc > 2) ? atoi(argv[2]) : 4000000;
    int numThreads     = (argc > 3) ? atoi(argv[3]) : 8;
    EvalMode mode      = parse_mode(argc, argv);
    if (pairs_per_rate < 1) pairs_per_rate = 1;
    if (numThreads < 1) numThreads = 1;

    // ── KSSD parameters (shuffle dictionary generated in memory) ───────────
    Sketch::kssd_parameter_t kssdPara; // half_k=10, half_subk=6, drlevel=3

    vector<double> rates = {0.001, 0.005, 0.01, 0.02, 0.05,
                            0.1,   0.15,  0.2,  0.25, 0.3};
    int n_rates     = (int)rates.size();
    int total_pairs = n_rates * pairs_per_rate;

    double t0 = get_sec();

    // ═══════════════════════════════════════════════════════════════════════
    // §1  FIXED-LENGTH BENCHMARK
    // ═══════════════════════════════════════════════════════════════════════
    if (mode != EVAL_VAR_ONLY && mode != EVAL_LOWCOMP && mode != EVAL_FILE) {
    vector<PairResult> results(total_pairs);

    #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
    for (int idx = 0; idx < total_pairs; idx++) {
        int    ri = idx / pairs_per_rate;
        double p  = rates[ri];

        mt19937_64 rng(42ULL + (uint64_t)idx * 1000003ULL);

        vector<char> seq_a(seq_length + 1);
        vector<char> seq_b(seq_length + 1);
        gen_random_seq(seq_a.data(), seq_length, rng);

        // Generate mutant AND track which positions were substituted.
        vector<bool> sub_at;
        gen_mutant_tracked(seq_a.data(), seq_b.data(), seq_length, p, rng, sub_at);

        // Shannon entropy weights from seq_a (same k-mer → same entropy → same
        // weight in A and B when the k-mer matches, satisfying ProbMinHash's
        // requirement that the weight is intrinsic to the k-mer sequence).
        // For unmatched positions the k-mers differ, but the weight array index
        // [i] is the same; exact_wj_from_subs uses wA[i]+wB[i] = 2*w[i] there,
        // which is correct since both the A-side and B-side k-mers get weight w[i].
        const int k21    = 21;
        uint64_t  nwin21 = (uint64_t)(seq_length - k21 + 1);
        vector<double> w;
        fill_kmer_entropy_weights(w, seq_a.data(), (uint64_t)seq_length, k21);

        // Exact weighted Jaccard (O(L) sliding-window, no hash map).
        double gt_wj = exact_wj_from_subs(sub_at,
                                           w.data(), seq_length,
                                           w.data(), seq_length,
                                           k21);

        // ── HLL (k=32, np=14 → 16384 regs) ────────────────────────────
        Sketch::HyperLogLog h1(14), h2(14);
        h1.update(seq_a.data());
        h2.update(seq_b.data());
        double hll_j = h1.jaccard_index(h2);

        // ── SetSketch (k=32, np=14 → 16384 regs) ──────────────────────
        Sketch::SetSketch s1(14), s2(14);
        s1.update(seq_a.data());
        s2.update(seq_b.data());
        double ss_j = s1.jaccard_index(s2);

        // ── KSSD (k=20, half_k=10, half_subk=6, drlevel=3) ────────────
        Sketch::Kssd* k1 = new Sketch::Kssd(kssdPara);
        Sketch::Kssd* k2 = new Sketch::Kssd(kssdPara);
        k1->update(seq_a.data());
        k2->update(seq_b.data());
        double kssd_j = k1->jaccard(k2);
        delete k1;
        delete k2;

        // ── MinHash (k=21, size=1024) ──────────────────────────────────
        Sketch::MinHash m1(21, 1024, 42, true), m2(21, 1024, 42, true);
        m1.update(seq_a.data());
        m2.update(seq_b.data());
        m1.finalize();
        m2.finalize();
        double mh_j = m1.jaccard(&m2);

        // ── ProbMinHash4 unweighted (m=1024, k=21) ────────────────────
        Sketch::ProbMinHash4 pm1(1024, 21, 42), pm2(1024, 21, 42);
        pm1.update(seq_a.data(), seq_length);
        pm2.update(seq_b.data(), seq_length);
        double pmh_j = pm1.jaccard(pm2);

        // ── Route C: Top-L ProbMinHash (L=1, 2, 4) ─────────────────
        Sketch::ProbMinHash4 tl1a(1024,21,42,1), tl1b(1024,21,42,1);
        tl1a.update(seq_a.data(), seq_length);
        tl1b.update(seq_b.data(), seq_length);
        double tl1_j = tl1a.jaccard(tl1b);

        Sketch::ProbMinHash4 tl2a(1024,21,42,2), tl2b(1024,21,42,2);
        tl2a.update(seq_a.data(), seq_length);
        tl2b.update(seq_b.data(), seq_length);
        double tl2_j = tl2a.jaccard(tl2b);

        Sketch::ProbMinHash4 tl4a(1024,21,42,4), tl4b(1024,21,42,4);
        tl4a.update(seq_a.data(), seq_length);
        tl4b.update(seq_b.data(), seq_length);
        double tl4_j = tl4a.jaccard(tl4b);

        // ── ProbMinHash4OP – One-Permutation (m=1024, k=21) ─────────
        Sketch::ProbMinHash4OP op1(1024, 21, 42), op2(1024, 21, 42);
        op1.update(seq_a.data(), seq_length);
        op2.update(seq_b.data(), seq_length);
        double oph_j = op1.jaccard(op2);

        // ── FastKMV – bottom-k (k=1024, kmer=21) ───────────────────
        Sketch::FastKMV kv1(1024, 21, 42), kv2(1024, 21, 42);
        kv1.update(seq_a.data(), seq_length);
        kv2.update(seq_b.data(), seq_length);
        double kmv_j = kv1.jaccard(kv2);

        // ── BinDash (s=2^16, k=21, b=16) ────────────────────────────────
        Sketch::BinDash bd1(16, 21, 16, 42), bd2(16, 21, 16, 42);
        bd1.update(seq_a.data(), seq_length);
        bd2.update(seq_b.data(), seq_length);
        double bd_j = bd1.jaccard(bd2);

        // ── WeightedPMH (m=1024, k=21, Shannon entropy weights) ──────────
        Sketch::ProbMinHash4 wp1(1024, 21, 42), wp2(1024, 21, 42);
        wp1.updateEntropy(seq_a.data(), seq_length);
        wp2.updateEntropy(seq_b.data(), seq_length);
        double wpmh_j = wp1.jaccard(wp2);

        PairResult& r = results[idx];
        r.rate     = p;
        r.theo_j32 = theo_jaccard(p, 32);
        r.theo_j21 = theo_jaccard(p, 21);
        r.theo_j20 = theo_jaccard(p, 20);
        r.hll_j    = hll_j;
        r.ss_j     = ss_j;
        r.kssd_j   = kssd_j;
        r.mh_j     = mh_j;
        r.pmh_j    = pmh_j;
        r.oph_j    = oph_j;
        r.kmv_j    = kmv_j;
        r.bd_j     = bd_j;
        r.tl1_j    = tl1_j;
        r.tl2_j    = tl2_j;
        r.tl4_j    = tl4_j;
        r.wpmh_j   = wpmh_j;
        r.gt_wj    = gt_wj;
    }

    double t1 = get_sec();

    // ── Accumulate errors ──────────────────────────────────────────────────
    // Per-rate accumulators
    struct RateAcc {
        double hll=0, ss=0, kd=0, mh=0, pmh=0, oph=0, kmv=0, bd=0, wpmh=0;
        double hll2=0, ss2=0, kd2=0, mh2=0, pmh2=0, oph2=0, kmv2=0, bd2=0, wpmh2=0;
    };
    vector<RateAcc> acc(n_rates);
    for (int i = 0; i < total_pairs; i++) {
        int ri = i / pairs_per_rate;
        const PairResult& r = results[i];
        double j32 = theo_jaccard(rates[ri], 32);
        double j21 = theo_jaccard(rates[ri], 21);
        double j20 = theo_jaccard(rates[ri], 20);
        auto sq = [](double x){ return x*x; };
        double eh=r.hll_j-j32,  es=r.ss_j-j32,  ek=r.kssd_j-j20;
        double em=r.mh_j-j21,   ep=r.pmh_j-j21, eo=r.oph_j-j21, eb=r.bd_j-j21;
        double ev=r.kmv_j-j21,  ew=r.wpmh_j-r.gt_wj;
        acc[ri].hll  +=fabs(eh); acc[ri].hll2  +=sq(eh);
        acc[ri].ss   +=fabs(es); acc[ri].ss2   +=sq(es);
        acc[ri].kd   +=fabs(ek); acc[ri].kd2   +=sq(ek);
        acc[ri].mh   +=fabs(em); acc[ri].mh2   +=sq(em);
        acc[ri].pmh  +=fabs(ep); acc[ri].pmh2  +=sq(ep);
        acc[ri].oph  +=fabs(eo); acc[ri].oph2  +=sq(eo);
        acc[ri].bd   +=fabs(eb); acc[ri].bd2   +=sq(eb);
        acc[ri].kmv  +=fabs(ev); acc[ri].kmv2  +=sq(ev);
        acc[ri].wpmh +=fabs(ew); acc[ri].wpmh2 +=sq(ew);
    }

    // ── Fixed-length summary table ─────────────────────────────────────────
    // ground truth: theo_jaccard(p,k) for each sketch's k; wPMH vs exact WJ
    fprintf(stderr,
        "\n=== Fixed-Length (L=%d bp, m=1024, %d pairs×%d rates) ===\n"
        "        ──── k=32 ────  ─ k=20 ─  ───────────────────────── k=21 ─────────────────────────  ── wPMH ──\n"
        "rate      HLL    SS     KSSD     MinHash  ProbMH  OnePerm    KMV   BinDash   (vs exact WJ)\n"
        "        [MAE]  [MAE]   [MAE]     [MAE]    [MAE]    [MAE]   [MAE]   [MAE]       [MAE]\n"
        "───────  ──────────────────────────────────────────────────────────────────────────────────────\n",
        seq_length, pairs_per_rate, n_rates);

    RateAcc g;   // global totals
    for (int ri = 0; ri < n_rates; ri++) {
        double N = (double)pairs_per_rate;
        const RateAcc& a = acc[ri];
        fprintf(stderr,
            "%.3f   %6.4f %6.4f  %6.4f   %7.4f  %7.4f  %7.4f  %7.4f  %7.4f    %7.4f\n",
            rates[ri],
            a.hll/N, a.ss/N, a.kd/N,
            a.mh/N, a.pmh/N, a.oph/N, a.kmv/N, a.bd/N,
            a.wpmh/N);
        // Machine-readable per-rate line for parsing
        fprintf(stderr,
            "PARSE_RATE len=%d rate=%.4f HLL=%.6f SetSketch=%.6f KSSD=%.6f"
            " MinHash=%.6f ProbMH=%.6f OnePerm=%.6f KMV=%.6f BinDash=%.6f WgtPMH=%.6f\n",
            seq_length, rates[ri],
            a.hll/N, a.ss/N, a.kd/N,
            a.mh/N, a.pmh/N, a.oph/N, a.kmv/N, a.bd/N, a.wpmh/N);
        g.hll+=a.hll; g.ss+=a.ss; g.kd+=a.kd;
        g.mh+=a.mh;   g.pmh+=a.pmh; g.oph+=a.oph; g.kmv+=a.kmv; g.bd+=a.bd;
        g.wpmh+=a.wpmh;
        g.hll2+=a.hll2; g.ss2+=a.ss2; g.kd2+=a.kd2;
        g.mh2+=a.mh2; g.pmh2+=a.pmh2; g.oph2+=a.oph2; g.kmv2+=a.kmv2; g.bd2+=a.bd2;
        g.wpmh2+=a.wpmh2;
    }
    double T = (double)total_pairs;
    fprintf(stderr,
        "───────  ──────────────────────────────────────────────────────────────────────────────────────\n"
        "MAE    %6.4f %6.4f  %6.4f   %7.4f  %7.4f  %7.4f  %7.4f  %7.4f    %7.4f\n"
        "RMSE   %6.4f %6.4f  %6.4f   %7.4f  %7.4f  %7.4f  %7.4f  %7.4f    %7.4f\n"
        "time: %.1f s\n",
        g.hll/T, g.ss/T, g.kd/T,
        g.mh/T, g.pmh/T, g.oph/T, g.kmv/T, g.bd/T, g.wpmh/T,
        sqrt(g.hll2/T), sqrt(g.ss2/T), sqrt(g.kd2/T),
        sqrt(g.mh2/T), sqrt(g.pmh2/T), sqrt(g.oph2/T), sqrt(g.kmv2/T), sqrt(g.bd2/T), sqrt(g.wpmh2/T),
        t1 - t0);

    if (mode == EVAL_FIXED_ONLY)
        return 0;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // §2  VARIABLE-LENGTH BENCHMARK  (mixed-size pairs)
    // ═══════════════════════════════════════════════════════════════════════
    if (mode != EVAL_LOWCOMP && mode != EVAL_FILE) {
    struct VarPairResult {
        double rate;
        int    len_a, len_b;
        double theo_j32, theo_j21, theo_j20;
        double hll_j, ss_j, kssd_j, mh_j, pmh_j, oph_j, kmv_j, bd_j;
        double tl1_j, tl2_j, tl4_j;
        double wpmh_j;   // weighted PMH estimate
        double gt_wj;    // exact weighted Jaccard
    };

    auto theo_jaccard_vl = [](double p, int k, int La, int Lb) -> double {
        double S      = (double)min(La, Lb);
        double shared = S * pow(1.0 - p, k);
        return shared / (La + Lb - shared);
    };

    vector<VarPairResult> vr(total_pairs);

    {
        mt19937_64 lrng(0xDEADBEEFCAFEULL);
        uniform_real_distribution<double> log_ld(log(50000.0), log(8000000.0));
        for (int i = 0; i < total_pairs; i++) {
            vr[i].len_a = (int)round(exp(log_ld(lrng)));
            vr[i].len_b = (int)round(exp(log_ld(lrng)));
        }
    }

    // (lengths already set above)

    double t2 = get_sec();

    #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
    for (int idx = 0; idx < total_pairs; idx++) {
        int    ri2 = idx / pairs_per_rate;
        double p   = rates[ri2];
        int    La  = vr[idx].len_a;
        int    Lb  = vr[idx].len_b;
        int    S   = min(La, Lb);

        mt19937_64 rng(0xCAFEBABEULL + (uint64_t)idx * 1000033ULL);

        vector<char> seq_a(La + 1), seq_b(Lb + 1);

        gen_random_seq(seq_a.data(), La, rng);

        // Mutate the shared prefix S and track substitutions.
        vector<bool> sub_at;
        gen_mutant_tracked(seq_a.data(), seq_b.data(), S, p, rng, sub_at);

        // Tail extension (when Lb > La): fresh random bases.
        if (Lb > S)
            gen_random_seq(seq_b.data() + S, Lb - S, rng);
        seq_b[Lb] = '\0';

        // Shannon entropy weights from seq_a (A and B use the same weight for
        // each position; the weight of a matched k-mer is identical in both).
        // For positions where A is longer than B, we compute entropy from A
        // for the full La-k+1 positions; B uses the same array up to Lb-k.
        const int k21    = 21;
        uint64_t  nwinA  = (La >= k21) ? (uint64_t)(La - k21 + 1) : 0;
        uint64_t  nwinB  = (Lb >= k21) ? (uint64_t)(Lb - k21 + 1) : 0;
        vector<double> w;
        if (nwinA > 0) fill_kmer_entropy_weights(w, seq_a.data(), (uint64_t)La, k21);

        // Exact weighted Jaccard — needs the correct weight array for B.
        // When Lb > La, w only has (La-k+1) elements, so we must use a
        // separate wB array of size (Lb-k+1) to avoid out-of-bounds reads.
        // wB_ext is also reused below for the weighted PMH sketch.
        vector<double> wB_ext;
        if (nwinB > nwinA && nwinA > 0)
            fill_kmer_entropy_weights(wB_ext, seq_b.data(), (uint64_t)Lb, k21);
        const double* wB_ptr = (nwinB > nwinA && nwinA > 0)
                               ? wB_ext.data() : w.data();

        double gt_wj = 0.0;
        if (nwinA > 0 && nwinB > 0)
            gt_wj = exact_wj_from_subs(sub_at, w.data(), La, wB_ptr, Lb, k21);

        Sketch::HyperLogLog  h1(14), h2(14);
        h1.update(seq_a.data()); h2.update(seq_b.data());
        double hll_j = h1.jaccard_index(h2);

        Sketch::SetSketch    s1(14), s2(14);
        s1.update(seq_a.data()); s2.update(seq_b.data());
        double ss_j = s1.jaccard_index(s2);

        Sketch::Kssd* k1 = new Sketch::Kssd(kssdPara);
        Sketch::Kssd* k2 = new Sketch::Kssd(kssdPara);
        k1->update(seq_a.data()); k2->update(seq_b.data());
        double kssd_j = k1->jaccard(k2);
        delete k1; delete k2;

        Sketch::MinHash m1(21, 1024, 42, true), m2(21, 1024, 42, true);
        m1.update(seq_a.data()); m2.update(seq_b.data());
        m1.finalize(); m2.finalize();
        double mh_j = m1.jaccard(&m2);

        Sketch::ProbMinHash4 pm1(1024, 21, 42), pm2(1024, 21, 42);
        pm1.update(seq_a.data(), La); pm2.update(seq_b.data(), Lb);
        double pmh_j = pm1.jaccard(pm2);

        Sketch::ProbMinHash4 vtl1a(1024,21,42,1), vtl1b(1024,21,42,1);
        vtl1a.update(seq_a.data(), La); vtl1b.update(seq_b.data(), Lb);
        double vtl1_j = vtl1a.jaccard(vtl1b);

        Sketch::ProbMinHash4 vtl2a(1024,21,42,2), vtl2b(1024,21,42,2);
        vtl2a.update(seq_a.data(), La); vtl2b.update(seq_b.data(), Lb);
        double vtl2_j = vtl2a.jaccard(vtl2b);

        Sketch::ProbMinHash4 vtl4a(1024,21,42,4), vtl4b(1024,21,42,4);
        vtl4a.update(seq_a.data(), La); vtl4b.update(seq_b.data(), Lb);
        double vtl4_j = vtl4a.jaccard(vtl4b);

        Sketch::ProbMinHash4OP op1(1024, 21, 42), op2(1024, 21, 42);
        op1.update(seq_a.data(), La); op2.update(seq_b.data(), Lb);
        double oph_j = op1.jaccard(op2);

        Sketch::FastKMV kv1(1024, 21, 42), kv2(1024, 21, 42);
        kv1.update(seq_a.data(), La); kv2.update(seq_b.data(), Lb);
        double kmv_j = kv1.jaccard(kv2);

        Sketch::BinDash bd1(16, 21, 16, 42), bd2(16, 21, 16, 42);
        bd1.update(seq_a.data(), La); bd2.update(seq_b.data(), Lb);
        double bd_j = bd1.jaccard(bd2);

        // Weighted PMH — fused entropy-weighted update (no extra allocation).
        double wpmh_j = 0.0;
        if (nwinA > 0 && nwinB > 0) {
            Sketch::ProbMinHash4 wp1(1024, 21, 42), wp2(1024, 21, 42);
            wp1.updateEntropy(seq_a.data(), La);
            wp2.updateEntropy(seq_b.data(), Lb);
            wpmh_j = wp1.jaccard(wp2);
        }

        VarPairResult& r = vr[idx];
        r.rate     = p;
        r.len_a    = La;
        r.len_b    = Lb;
        r.theo_j32 = theo_jaccard_vl(p, 32, La, Lb);
        r.theo_j21 = theo_jaccard_vl(p, 21, La, Lb);
        r.theo_j20 = theo_jaccard_vl(p, 20, La, Lb);
        r.hll_j    = hll_j;
        r.ss_j     = ss_j;
        r.kssd_j   = kssd_j;
        r.mh_j     = mh_j;
        r.pmh_j    = pmh_j;
        r.oph_j    = oph_j;
        r.kmv_j    = kmv_j;
        r.bd_j     = bd_j;
        r.tl1_j    = vtl1_j;
        r.tl2_j    = vtl2_j;
        r.tl4_j    = vtl4_j;
        r.wpmh_j   = wpmh_j;
        r.gt_wj    = gt_wj;
    }

    double t3 = get_sec();

    // ── Variable-length summary table ─────────────────────────────────────
    double v_hll=0,v_ss=0,v_kd=0,v_mh=0,v_pmh=0,v_oph=0,v_kmv=0,v_bd=0,v_wpmh=0;
    double v_hll2=0,v_ss2=0,v_kd2=0,v_mh2=0,v_pmh2=0,v_oph2=0,v_kmv2=0,v_bd2=0,v_wpmh2=0;
    for (int i = 0; i < total_pairs; i++) {
        const VarPairResult& r = vr[i];
        auto sq = [](double x){ return x*x; };
        double eh=r.hll_j-r.theo_j32,  es=r.ss_j-r.theo_j32,  ek=r.kssd_j-r.theo_j20;
        double em=r.mh_j-r.theo_j21,   ep=r.pmh_j-r.theo_j21, eo=r.oph_j-r.theo_j21, eb=r.bd_j-r.theo_j21;
        double ev=r.kmv_j-r.theo_j21,  ew=r.wpmh_j-r.gt_wj;
        v_hll+=fabs(eh); v_ss+=fabs(es); v_kd+=fabs(ek);
        v_mh+=fabs(em);  v_pmh+=fabs(ep); v_oph+=fabs(eo);
        v_kmv+=fabs(ev); v_bd+=fabs(eb); v_wpmh+=fabs(ew);
        v_hll2+=sq(eh); v_ss2+=sq(es); v_kd2+=sq(ek);
        v_mh2+=sq(em);  v_pmh2+=sq(ep); v_oph2+=sq(eo);
        v_kmv2+=sq(ev); v_bd2+=sq(eb); v_wpmh2+=sq(ew);
    }
    double Nd = (double)total_pairs;
    fprintf(stderr,
        "\n=== Variable-Length ([50k–8M bp], m=1024, %d pairs×%d rates) ===\n"
        "        ──── k=32 ────  ─ k=20 ─  ───────────────────────── k=21 ─────────────────────────  ── wPMH ──\n"
        "          HLL    SS     KSSD     MinHash  ProbMH  OnePerm    KMV   BinDash   (vs exact WJ)\n"
        "───────  ──────────────────────────────────────────────────────────────────────────────────────\n"
        "MAE    %6.4f %6.4f  %6.4f   %7.4f  %7.4f  %7.4f  %7.4f  %7.4f    %7.4f\n"
        "RMSE   %6.4f %6.4f  %6.4f   %7.4f  %7.4f  %7.4f  %7.4f  %7.4f    %7.4f\n"
        "time: %.1f s\n",
        pairs_per_rate, n_rates,
        v_hll/Nd, v_ss/Nd, v_kd/Nd,
        v_mh/Nd, v_pmh/Nd, v_oph/Nd, v_kmv/Nd, v_bd/Nd, v_wpmh/Nd,
        sqrt(v_hll2/Nd), sqrt(v_ss2/Nd), sqrt(v_kd2/Nd),
        sqrt(v_mh2/Nd), sqrt(v_pmh2/Nd), sqrt(v_oph2/Nd), sqrt(v_kmv2/Nd), sqrt(v_bd2/Nd), sqrt(v_wpmh2/Nd),
        t3 - t2);

    }   // end §2

    // ═══════════════════════════════════════════════════════════════════════
    // §3  LOW-COMPLEXITY CONTAMINATION BENCHMARK
    //
    // Evaluates the advantage of Shannon-entropy-weighted Jaccard over
    // unweighted Jaccard when sequences contain a shared low-complexity region.
    //
    // Sequence structure:
    //   seq_A = [biased region (L_low bp, IDENTICAL in A and B)]
    //           + [clean ACGT region (L_clean bp, mutated at rate p)]
    //   seq_B = same structure
    //
    // Biased region: 90% A, 10% {C,G,T}  → k-mer Shannon entropy ≈ 0.25 bit
    //                                     → weight ≈ 0.1–0.22 (vs ~1.0 for clean)
    //
    // Ground truth: unweighted Jaccard of ONLY the clean-region k-mers.
    //
    // Demonstration:
    //   ProbMH (unweighted) overestimates J because it gives full weight (1.0)
    //   to the shared biased-region k-mers.
    //   wPMH (weighted)    is much closer to the clean-only truth because the
    //   biased-region k-mers have low entropy → low weight → small contribution.
    // ═══════════════════════════════════════════════════════════════════════

    if (mode != EVAL_LOWCOMP && mode != EVAL_FILE)
        return 0;

    // ── §3 runs only in EVAL_LOWCOMP mode ─────────────────────────────────
    if (mode == EVAL_LOWCOMP)
    {
    const double  A_frac    = 0.90;   // 90% A in the biased region
    const int     k21       = 21;
    // Test biased region fractions: 0% (baseline), 25%, 50%, 75%
    const vector<double> low_fracs = {0.0, 0.25, 0.50, 0.75};

    fprintf(stderr,
        "\n=== Low-complexity Contamination (L=%d bp, m=1024, %d pairs×%d rates) ===\n"
        "  Biased region: %.0f%% A + %.0f%% {{C,G,T}} (entropy≈0.25bit, weight≈0.12-0.22)\n"
        "  Clean region:  uniform ACGT (entropy≈2.0bit, weight≈1.0)\n"
        "  Ground truth:  Jaccard of clean-region k-mers only\n"
        "  threads: %d\n",
        seq_length, pairs_per_rate, n_rates,
        A_frac * 100, (1.0 - A_frac) * 100, numThreads);

    double t4 = get_sec();

    for (double low_frac : low_fracs) {
        int L_low   = (int)(seq_length * low_frac);
        int L_clean = seq_length - L_low;

        struct LCResult {
            double rate, j_clean, pmh_j, wpmh_j;
        };
        vector<LCResult> lc(total_pairs);

        #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
        for (int idx = 0; idx < total_pairs; idx++) {
            int    ri = idx / pairs_per_rate;
            double p  = rates[ri];

            mt19937_64 rng(0xBEEFCAFEULL + (uint64_t)idx * 999983ULL);

            vector<char> seq_a(seq_length + 1, 'A');
            vector<char> seq_b(seq_length + 1, 'A');

            // Biased region: IDENTICAL in A and B (shared low-complexity prefix)
            if (L_low > 0) {
                gen_biased_seq(seq_a.data(), L_low, rng, A_frac);
                memcpy(seq_b.data(), seq_a.data(), (size_t)L_low);
            }

            // Clean region: A is random, B is mutated from A at rate p
            vector<bool> sub_at_clean;
            if (L_clean >= k21) {
                gen_random_seq(seq_a.data() + L_low, L_clean, rng);
                gen_mutant_tracked(seq_a.data() + L_low,
                                   seq_b.data() + L_low,
                                   L_clean, p, rng, sub_at_clean);
            }
            seq_a[seq_length] = '\0';
            seq_b[seq_length] = '\0';

            // Ground truth: exact Jaccard of clean k-mers only
            double j_clean = (L_clean >= k21)
                             ? exact_j_clean(sub_at_clean, L_clean, k21)
                             : 0.0;

            // Unweighted ProbMH on full sequence
            Sketch::ProbMinHash4 pm1(1024, k21, 42), pm2(1024, k21, 42);
            pm1.update(seq_a.data(), seq_length);
            pm2.update(seq_b.data(), seq_length);
            double pmh_j = pm1.jaccard(pm2);

            // Weighted PMH on full sequence (low-complexity k-mers downweighted)
            Sketch::ProbMinHash4 wp1(1024, k21, 42), wp2(1024, k21, 42);
            wp1.updateEntropy(seq_a.data(), seq_length);
            wp2.updateEntropy(seq_b.data(), seq_length);
            double wpmh_j = wp1.jaccard(wp2);

            lc[idx] = {p, j_clean, pmh_j, wpmh_j};
        }

        // ── Per-rate summary ────────────────────────────────────────────────
        fprintf(stderr,
            "\n  biased=%.0f%% (L_low=%d bp), clean=%.0f%% (L_clean=%d bp)\n"
            "  rate    J_clean   PMH_MAE   wPMH_MAE  PMH_bias  wPMH_bias  Advantage\n"
            "  ─────────────────────────────────────────────────────────────────────\n",
            low_frac * 100, L_low, (1.0 - low_frac) * 100, L_clean);

        double g_pmh = 0, g_wpmh = 0;
        for (int ri = 0; ri < n_rates; ri++) {
            int N = pairs_per_rate;
            double j_mean = 0, pmh_ae = 0, wpmh_ae = 0, pmh_bi = 0, wpmh_bi = 0;
            for (int pi = 0; pi < N; pi++) {
                const LCResult& r = lc[ri * N + pi];
                j_mean  += r.j_clean;
                double ep = r.pmh_j  - r.j_clean;
                double ew = r.wpmh_j - r.j_clean;
                pmh_ae  += fabs(ep); pmh_bi  += ep;
                wpmh_ae += fabs(ew); wpmh_bi += ew;
            }
            g_pmh += pmh_ae; g_wpmh += wpmh_ae;
            double adv = (pmh_ae > 1e-9)
                         ? 100.0 * (1.0 - wpmh_ae / pmh_ae) : 0.0;
            fprintf(stderr,
                "  %.3f  %7.4f  %8.4f  %9.4f  %+8.4f  %+9.4f  %+.1f%%\n",
                rates[ri], j_mean / N,
                pmh_ae / N, wpmh_ae / N,
                pmh_bi / N, wpmh_bi / N,
                adv);
        }
        double adv_g = (g_pmh > 1e-9) ? 100.0 * (1.0 - g_wpmh / g_pmh) : 0.0;
        fprintf(stderr,
            "  ─────────────────────────────────────────────────────────────────────\n"
            "  GLOBAL    ---   %8.4f  %9.4f  --- wPMH saves %+.1f%% error\n",
            g_pmh / total_pairs, g_wpmh / total_pairs, adv_g);
    }

    fprintf(stderr, "\nLow-complexity total time: %.1f s\n", get_sec() - t4);
    }   // end §3

    // ═══════════════════════════════════════════════════════════════════════
    // §4  FILE-BASED BENCHMARK
    //
    // Reads actual FASTA files generated by gen_bench_data from disk.
    // Usage:  exe_eval_sketch_acc <pairs_per_rate> 0 <threads> file <bench_dir>
    //
    // bench_dir must contain:
    //   pairs_meta.tsv        — pair_id, rate, len_bp, theo_j_k21, theo_j_k32, theo_j_k20
    //   pairs_exact_j.tsv     — pair_id, rate, len_bp, exact_j_k21
    //   pairNNNNNN_A.fa / pairNNNNNN_B.fa  — individual FASTA files
    //
    // Ground truth:
    //   k=21 algorithms (MinHash, PMH, OnePerm, FastKMV): exact_j_k21
    //   KSSD (k=20): theo_j_k20
    //   HLL / SetSketch (k=32): theo_j_k32
    // ═══════════════════════════════════════════════════════════════════════
    if (mode != EVAL_FILE)
        return 0;

    // ── parse bench_dir ───────────────────────────────────────────────────
    if (argc < 6) {
        fprintf(stderr,
            "ERROR: file mode requires bench_dir as 6th argument.\n"
            "Usage: %s <pairs_per_rate> 0 <threads> file <bench_dir>\n", argv[0]);
        return 1;
    }
    std::string bench_dir = argv[5];
    if (bench_dir.back() == '/') bench_dir.pop_back();

    // ── helper: read one FASTA sequence from file ─────────────────────────
    auto read_fasta = [](const std::string& path) -> std::string {
        FILE* f = fopen(path.c_str(), "r");
        if (!f) return "";
        std::string seq;
        seq.reserve(8 * 1024 * 1024);
        char buf[65536];
        while (fgets(buf, sizeof(buf), f)) {
            if (buf[0] == '>') continue;
            int len = (int)strlen(buf);
            while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) len--;
            seq.append(buf, (size_t)len);
        }
        fclose(f);
        return seq;
    };

    // ── load pairs_meta.tsv ───────────────────────────────────────────────
    struct PairMeta {
        double rate, theo_j21, theo_j32, theo_j20;
    };
    std::unordered_map<int, PairMeta> meta_map;
    {
        std::ifstream fin(bench_dir + "/pairs_meta.tsv");
        if (!fin.is_open()) {
            fprintf(stderr, "ERROR: cannot open %s/pairs_meta.tsv\n", bench_dir.c_str());
            return 1;
        }
        std::string line;
        std::getline(fin, line);  // skip header
        while (std::getline(fin, line)) {
            std::istringstream ss(line);
            int pid; double rate; int len;
            double j21, j32, j20;
            char sep;
            ss >> pid >> sep >> rate >> sep >> len >> sep >> j21 >> sep >> j32 >> sep >> j20;
            if (meta_map.find(pid) == meta_map.end())
                meta_map[pid] = {rate, j21, j32, j20};
        }
    }

    // ── load pairs_exact_j.tsv ────────────────────────────────────────────
    std::unordered_map<int, double> exact_j21_map;
    {
        std::ifstream fin(bench_dir + "/pairs_exact_j.tsv");
        if (fin.is_open()) {
            std::string line;
            std::getline(fin, line);  // skip header
            while (std::getline(fin, line)) {
                std::istringstream ss(line);
                int pid; double rate; int len; double ej21;
                char sep;
                ss >> pid >> sep >> rate >> sep >> len >> sep >> ej21;
                exact_j21_map[pid] = ej21;
            }
        } else {
            fprintf(stderr, "WARN: pairs_exact_j.tsv not found, using theo_j_k21\n");
        }
    }

    // ── enumerate pairs to process ────────────────────────────────────────
    // pairs_per_rate pairs for each of the 10 rates
    // rate_i=0 → pair_id 0..pairs_per_rate-1 (rate=0.001)
    // rate_i=1 → pair_id 1000..1999 (rate=0.005), etc.
    int n_file_pairs = n_rates * pairs_per_rate;

    struct FileResult {
        int    pair_id;
        double rate;
        double gt_j21, gt_j20, gt_j32;
        double hll_j, ss_j, kssd_j, mh_j, pmh_j, oph_j, kmv_j, bd_j;
    };
    std::vector<FileResult> fres(n_file_pairs);

    fprintf(stderr,
        "\n=== File-Based Benchmark: %s (%d pairs×%d rates, threads=%d) ===\n",
        bench_dir.c_str(), pairs_per_rate, n_rates, numThreads);

    double tf0 = get_sec();
    std::atomic<int> io_errors(0);

    #pragma omp parallel for schedule(dynamic, 4) num_threads(numThreads)
    for (int idx = 0; idx < n_file_pairs; idx++) {
        int ri       = idx / pairs_per_rate;
        int pair_id  = ri * 1000 + (idx % pairs_per_rate);  // 1000 pairs per rate in bench data

        // ── read sequences ──────────────────────────────────────────────
        char pa[512], pb[512];
        snprintf(pa, sizeof(pa), "%s/pair%06d_A.fa", bench_dir.c_str(), pair_id);
        snprintf(pb, sizeof(pb), "%s/pair%06d_B.fa", bench_dir.c_str(), pair_id);

        std::string sa = read_fasta(pa);
        std::string sb = read_fasta(pb);
        if (sa.empty() || sb.empty()) {
            io_errors++;
            fres[idx] = {pair_id, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            continue;
        }

        // ── ground truth ────────────────────────────────────────────────
        double gt_j21 = 0, gt_j20 = 0, gt_j32 = 0;
        auto mit = meta_map.find(pair_id);
        if (mit != meta_map.end()) {
            gt_j21 = mit->second.theo_j21;
            gt_j20 = mit->second.theo_j20;
            gt_j32 = mit->second.theo_j32;
        }
        auto eit = exact_j21_map.find(pair_id);
        if (eit != exact_j21_map.end())
            gt_j21 = eit->second;  // prefer exact over theoretical

        double rate = (mit != meta_map.end()) ? mit->second.rate : rates[ri];

        // ── sketch computations ─────────────────────────────────────────
        Sketch::HyperLogLog h1(14), h2(14);
        h1.update(sa.data()); h2.update(sb.data());
        double hll_j = h1.jaccard_index(h2);

        Sketch::SetSketch ss1(14), ss2(14);
        ss1.update(sa.data()); ss2.update(sb.data());
        double ss_j = ss1.jaccard_index(ss2);

        {
            Sketch::Kssd* k1 = new Sketch::Kssd(kssdPara);
            Sketch::Kssd* k2 = new Sketch::Kssd(kssdPara);
            k1->update(sa.data()); k2->update(sb.data());
            fres[idx].kssd_j = k1->jaccard(k2);
            delete k1; delete k2;
        }

        Sketch::MinHash m1(21, 1024, 42, true), m2(21, 1024, 42, true);
        m1.update(sa.data()); m2.update(sb.data());
        m1.finalize(); m2.finalize();
        double mh_j = m1.jaccard(&m2);

        Sketch::ProbMinHash4 pm1(1024, 21, 42), pm2(1024, 21, 42);
        pm1.updateEntropy(sa.data(), sa.size());
        pm2.updateEntropy(sb.data(), sb.size());
        double pmh_j = pm1.jaccard(pm2);

        Sketch::ProbMinHash4OP op1(1024, 21, 42), op2(1024, 21, 42);
        op1.update(sa.data(), (int)sa.size());
        op2.update(sb.data(), (int)sb.size());
        double oph_j = op1.jaccard(op2);

        Sketch::FastKMV kv1(1024, 21, 42), kv2(1024, 21, 42);
        kv1.update(sa.data(), (int)sa.size());
        kv2.update(sb.data(), (int)sb.size());
        double kmv_j = kv1.jaccard(kv2);

        Sketch::BinDash bd1(16, 21, 16, 42), bd2(16, 21, 16, 42);
        bd1.update(sa.data(), sa.size());
        bd2.update(sb.data(), sb.size());
        double bd_j = bd1.jaccard(bd2);

        fres[idx] = {pair_id, rate, gt_j21, gt_j20, gt_j32,
                     hll_j, ss_j, fres[idx].kssd_j,
                     mh_j, pmh_j, oph_j, kmv_j, bd_j};
    }

    double tf1 = get_sec();

    if (io_errors > 0)
        fprintf(stderr, "WARN: %d I/O errors (missing FASTA files)\n", (int)io_errors);

    // ── accumulate errors per rate ────────────────────────────────────────
    struct FAcc {
        double hll=0,ss=0,kd=0,mh=0,pmh=0,oph=0,kmv=0,bd=0;
        double hll2=0,ss2=0,kd2=0,mh2=0,pmh2=0,oph2=0,kmv2=0,bd2=0;
        int n=0;
    };
    std::vector<FAcc> facc(n_rates);
    auto sq = [](double x){ return x*x; };

    for (int idx = 0; idx < n_file_pairs; idx++) {
        const FileResult& r = fres[idx];
        if (r.gt_j21 == 0 && r.hll_j == 0) continue;  // skip failed reads
        int ri = idx / pairs_per_rate;
        FAcc& a = facc[ri];
        double eh = r.hll_j  - r.gt_j32;
        double es = r.ss_j   - r.gt_j32;
        double ek = r.kssd_j - r.gt_j20;
        double em = r.mh_j   - r.gt_j21;
        double ep = r.pmh_j  - r.gt_j21;
        double eo = r.oph_j  - r.gt_j21;
        double ev = r.kmv_j  - r.gt_j21;
        double eb = r.bd_j   - r.gt_j21;
        a.hll+=fabs(eh);  a.hll2+=sq(eh);
        a.ss +=fabs(es);  a.ss2 +=sq(es);
        a.kd +=fabs(ek);  a.kd2 +=sq(ek);
        a.mh +=fabs(em);  a.mh2 +=sq(em);
        a.pmh+=fabs(ep);  a.pmh2+=sq(ep);
        a.oph+=fabs(eo);  a.oph2+=sq(eo);
        a.kmv+=fabs(ev);  a.kmv2+=sq(ev);
        a.bd +=fabs(eb);  a.bd2 +=sq(eb);
        a.n++;
    }

    // ── print table ───────────────────────────────────────────────────────
    fprintf(stderr,
        "        ──── k=32 ────  ─ k=20 ─  ───────────────────────── k=21 (exact GT) ────────────\n"
        "rate      HLL    SS     KSSD     MinHash  ProbMH  OnePerm    KMV   BinDash\n"
        "        [MAE]  [MAE]   [MAE]     [MAE]    [MAE]    [MAE]   [MAE]   [MAE]\n"
        "───────  ────────────────────────────────────────────────────────────────────────────────\n");

    // infer representative length from bench_dir name (e.g. ".../L1M" → 1000000)
    int rep_len = seq_length;  // seq_length is 0 in file mode; parse from dir name
    if (rep_len == 0) {
        std::string dname = bench_dir;
        auto slash = dname.rfind('/');
        if (slash != std::string::npos) dname = dname.substr(slash + 1);
        // patterns: L100k, L500k, L1M, L2M, L4M, L8M
        if (!dname.empty() && (dname[0] == 'L' || dname[0] == 'l')) {
            std::string num = dname.substr(1);
            int mult = 1;
            if (!num.empty() && (num.back() == 'k' || num.back() == 'K'))
                { mult = 1000;    num.pop_back(); }
            else if (!num.empty() && (num.back() == 'm' || num.back() == 'M'))
                { mult = 1000000; num.pop_back(); }
            try { rep_len = std::stoi(num) * mult; } catch (...) { rep_len = 0; }
        }
        // fallback: read from first pair in pairs_meta
        if (rep_len == 0) {
            std::ifstream fm(bench_dir + "/pairs_meta.tsv");
            std::string hdr, row;
            if (fm && std::getline(fm, hdr) && std::getline(fm, row)) {
                auto cols = [&](const std::string& s) {
                    std::vector<std::string> v; std::istringstream ss(s);
                    std::string t; while (std::getline(ss, t, '\t')) v.push_back(t);
                    return v;
                }(row);
                if (cols.size() >= 3) try { rep_len = std::stoi(cols[2]); } catch (...) {}
            }
        }
    }

    FAcc gf;
    for (int ri = 0; ri < n_rates; ri++) {
        const FAcc& a = facc[ri];
        if (a.n == 0) continue;
        double N = (double)a.n;
        fprintf(stderr,
            "%.3f   %6.4f %6.4f  %6.4f   %7.4f  %7.4f  %7.4f  %7.4f  %7.4f\n",
            rates[ri],
            a.hll/N, a.ss/N, a.kd/N,
            a.mh/N, a.pmh/N, a.oph/N, a.kmv/N, a.bd/N);
        // machine-readable per-rate line (WgtPMH not available in file mode → 0)
        fprintf(stderr,
            "PARSE_RATE len=%d rate=%.4f HLL=%.6f SetSketch=%.6f KSSD=%.6f"
            " MinHash=%.6f ProbMH=%.6f OnePerm=%.6f KMV=%.6f BinDash=%.6f WgtPMH=%.6f\n",
            rep_len, rates[ri],
            a.hll/N, a.ss/N, a.kd/N,
            a.mh/N, a.pmh/N, a.oph/N, a.kmv/N, a.bd/N, 0.0);
        gf.hll+=a.hll; gf.ss+=a.ss; gf.kd+=a.kd;
        gf.mh+=a.mh;   gf.pmh+=a.pmh; gf.oph+=a.oph; gf.kmv+=a.kmv; gf.bd+=a.bd;
        gf.hll2+=a.hll2; gf.ss2+=a.ss2; gf.kd2+=a.kd2;
        gf.mh2+=a.mh2; gf.pmh2+=a.pmh2; gf.oph2+=a.oph2; gf.kmv2+=a.kmv2; gf.bd2+=a.bd2;
        gf.n+=a.n;
    }
    double TF = (double)gf.n;
    fprintf(stderr,
        "───────  ────────────────────────────────────────────────────────────────────────────────\n"
        "MAE    %6.4f %6.4f  %6.4f   %7.4f  %7.4f  %7.4f  %7.4f  %7.4f\n"
        "RMSE   %6.4f %6.4f  %6.4f   %7.4f  %7.4f  %7.4f  %7.4f  %7.4f\n"
        "time: %.1f s\n",
        gf.hll/TF, gf.ss/TF, gf.kd/TF,
        gf.mh/TF, gf.pmh/TF, gf.oph/TF, gf.kmv/TF, gf.bd/TF,
        sqrt(gf.hll2/TF), sqrt(gf.ss2/TF), sqrt(gf.kd2/TF),
        sqrt(gf.mh2/TF), sqrt(gf.pmh2/TF), sqrt(gf.oph2/TF), sqrt(gf.kmv2/TF), sqrt(gf.bd2/TF),
        tf1 - tf0);

    // machine-readable line for script parsing
    fprintf(stdout,
        "PARSE_FILE dir=%s n=%d MAE_HLL=%.4f MAE_SS=%.4f MAE_KSSD=%.4f "
        "MAE_MH=%.4f MAE_PMH=%.4f MAE_OPH=%.4f MAE_KMV=%.4f MAE_BD=%.4f\n",
        bench_dir.c_str(), gf.n,
        gf.hll/TF, gf.ss/TF, gf.kd/TF,
        gf.mh/TF, gf.pmh/TF, gf.oph/TF, gf.kmv/TF, gf.bd/TF);

    return 0;
}   // end §4
