/**
 * eval_alisim_sweep – Full algorithm × parameter sweep on alisim FASTA pairs.
 *
 * For one alisim dataset directory (model/length/branch_len):
 *   1. Load all pair_N_A.fa / pair_N_B.fa into memory.
 *   2. Read branch_len from pairs_meta.tsv.
 *   3. Precompute exact canonical-k-mer Jaccard for every unique k used in sweeps.
 *   4. Sweep every algorithm × every 1-D parameter grid (matching run_param_sweep_all.sh).
 *   5. Print one ALISIM_SWEEP line per (algo, params) combo to stdout.
 *
 * Usage:
 *   ./exe_eval_alisim_sweep --dir PATH --threads T [--model MODEL_NAME]
 *
 * Output columns (space-separated key=value):
 *   ALISIM_SWEEP dir model branch_len algo kmer size bits ss_a ss_base
 *                bd_bbits bd_s64 drlevel n mean_exact_J mean_est_J MAE RMSE
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
#include <fstream>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

// ─── FASTA reader ──────────────────────────────────────────────────────────────

static std::string read_fasta_seq(const std::string& path)
{
    std::ifstream f(path);
    if (!f) return {};
    std::string line, seq;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '>') continue;
        seq += line;
    }
    return seq;
}

// ─── Exact canonical k-mer Jaccard ─────────────────────────────────────────────

static double exact_kmer_jaccard(const std::string& sa, const std::string& sb, int k)
{
    if ((int)sa.size() < k || (int)sb.size() < k) return -1.0;

    auto get_canonical_kmers = [k](const std::string& seq) -> std::vector<uint64_t> {
        std::vector<uint64_t> out;
        out.reserve(seq.size());
        uint64_t fwd = 0, rc = 0;
        const uint64_t mask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
        int valid = 0;
        for (char c : seq) {
            int b = -1;
            switch (c | 32) {
                case 'a': b = 0; break; case 'c': b = 1; break;
                case 'g': b = 2; break; case 't': b = 3; break;
            }
            if (b < 0) { valid = 0; fwd = rc = 0; continue; }
            fwd = ((fwd << 2) | (uint64_t)b) & mask;
            rc  = (rc >> 2) | ((uint64_t)(3 - b) << (2 * (k - 1)));
            if (++valid >= k)
                out.push_back(std::min(fwd, rc));
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    };

    auto va = get_canonical_kmers(sa);
    auto vb = get_canonical_kmers(sb);

    size_t inter = 0, i = 0, j = 0;
    while (i < va.size() && j < vb.size()) {
        if      (va[i] == vb[j]) { ++inter; ++i; ++j; }
        else if (va[i] <  vb[j]) ++i;
        else                     ++j;
    }
    size_t uni = va.size() + vb.size() - inter;
    return uni > 0 ? (double)inter / (double)uni : 0.0;
}

// ─── KSSD shared-param cache ────────────────────────────────────────────────────

static std::shared_ptr<const Sketch::kssd_parameter_t>
get_kssd_params(int half_k, int half_subk, int drlevel)
{
    static std::mutex mtx;
    static std::map<std::tuple<int,int,int>,
                    std::shared_ptr<Sketch::kssd_parameter_t>> cache;
    std::lock_guard<std::mutex> lk(mtx);
    auto key = std::make_tuple(half_k, half_subk, drlevel);
    auto it  = cache.find(key);
    if (it != cache.end()) return it->second;
    try {
        auto p = std::make_shared<Sketch::kssd_parameter_t>(half_k, half_subk, drlevel);
        cache[key] = p;
        return p;
    } catch (...) {
        cache[key] = nullptr; // invalid combo – remember and skip
        return nullptr;
    }
}

// ─── Parameter combo descriptor ────────────────────────────────────────────────

struct Combo {
    std::string algo;
    int kmer    = 21;
    int size    = 1024;
    int bits    = 13;
    double ss_a    = 20.0;
    double ss_base = 2.0;
    int bd_bbits  = 16;
    int bd_s64    = 16;
    int drlevel   = 3;
};

// Build the complete list of 1-D parameter sweep combos.
//
// All k-mer values are capped at 31 to avoid the 64-bit overflow at k>=32
// (a 2-bit-packed k-mer of length 32 needs exactly 64 bits and mask
// computation `(1ULL<<(2*k))-1` is undefined behavior). Each algorithm's
// "default" combo is registered exactly once so the analyzer can always
// find the reference row.
static std::vector<Combo> make_all_combos()
{
    std::vector<Combo> combos;
    auto push = [&](Combo c){ combos.push_back(c); };

    // Common k-mer grid (all algos) — odd values from 11 to 31
    const std::vector<int> KMERS = {11,13,15,17,19,21,23,25,27,29,31};

    // ── MinHash ────────────────────────────────────────────────────────────────
    // Default: kmer=21, size=1024
    {
        Combo c; c.algo="minhash"; c.kmer=21; c.size=1024; push(c);
    }
    for (int k : KMERS) {
        if (k == 21) continue;
        Combo c; c.algo="minhash"; c.kmer=k; c.size=1024; push(c);
    }
    for (int s : {64,128,256,512,2048,4096,8192}) {
        Combo c; c.algo="minhash"; c.kmer=21; c.size=s; push(c);
    }

    // ── FastKMV ───────────────────────────────────────────────────────────────
    // Default: kmer=21, size=1024
    {
        Combo c; c.algo="fastkmv"; c.kmer=21; c.size=1024; push(c);
    }
    for (int k : KMERS) {
        if (k == 21) continue;
        Combo c; c.algo="fastkmv"; c.kmer=k; c.size=1024; push(c);
    }
    for (int s : {64,128,256,512,2048,4096,8192}) {
        Combo c; c.algo="fastkmv"; c.kmer=21; c.size=s; push(c);
    }

    // ── HLL ───────────────────────────────────────────────────────────────────
    // Default: kmer=21, bits=13
    {
        Combo c; c.algo="hll"; c.kmer=21; c.bits=13; push(c);
    }
    for (int k : KMERS) {
        if (k == 21) continue;
        Combo c; c.algo="hll"; c.kmer=k; c.bits=13; push(c);
    }
    for (int b : {8,9,10,11,12,14,15,16}) {
        Combo c; c.algo="hll"; c.kmer=21; c.bits=b; push(c);
    }

    // ── SetSketch ─────────────────────────────────────────────────────────────
    // At each anchor a: default (k=21,bits=13) + k sweep + bits sweep (same grid as a=20).
    const std::vector<double> SS_A_GRID = {0.5, 1.0, 5.0, 10.0, 15.0, 20.0};
    for (double a : SS_A_GRID) {
        {
            Combo c; c.algo="setsketch"; c.kmer=21; c.bits=13;
            c.ss_a=a; c.ss_base=2.0; push(c);
        }
        for (int k : KMERS) {
            if (k == 21) continue;
            Combo c; c.algo="setsketch"; c.kmer=k; c.bits=13;
            c.ss_a=a; c.ss_base=2.0; push(c);
        }
        for (int b : {8,9,10,11,12,14,15,16}) {
            Combo c; c.algo="setsketch"; c.kmer=21; c.bits=b;
            c.ss_a=a; c.ss_base=2.0; push(c);
        }
    }
    // base (b) sweep at a=20 only
    for (double base : {1.2,1.5,1.8,2.5,3.0,4.0,5.0,8.0}) {
        Combo c; c.algo="setsketch"; c.kmer=21; c.bits=13;
        c.ss_a=20.0; c.ss_base=base; push(c);
    }

    // ── SetSketch (joint MLE estimator, Ertl 2021) ────────────────────────────
    // Mirror every SetSketch combo with its joint-MLE counterpart so a single
    // run produces both inclusion-exclusion and MLE rows for direct comparison.
    {
        size_t n = combos.size();
        for (size_t i = 0; i < n; ++i) {
            if (combos[i].algo == "setsketch") {
                Combo c = combos[i];
                c.algo = "setsketch_mle";
                push(c);
            }
        }
    }

    // ── BinDash ───────────────────────────────────────────────────────────────
    // Default: kmer=21, bbits=16, s64=16
    {
        Combo c; c.algo="bindash"; c.kmer=21; c.bd_bbits=16; c.bd_s64=16; push(c);
    }
    for (int k : KMERS) {
        if (k == 21) continue;
        Combo c; c.algo="bindash"; c.kmer=k; c.bd_bbits=16; c.bd_s64=16; push(c);
    }
    for (int bb : {2,4,6,8,10,12,14,20,24,32}) {
        Combo c; c.algo="bindash"; c.kmer=21; c.bd_bbits=bb; c.bd_s64=16; push(c);
    }
    for (int s : {2,4,6,8,10,12,20,24,32,64}) {
        Combo c; c.algo="bindash"; c.kmer=21; c.bd_bbits=16; c.bd_s64=s; push(c);
    }

    // ── KSSD ──────────────────────────────────────────────────────────────────
    // Default: drlevel=3 (constraint: half_subk - drlevel >= 3 → drlevel <= 3)
    {
        Combo c; c.algo="kssd"; c.drlevel=3; push(c);
    }
    for (int d : {1,2}) {
        Combo c; c.algo="kssd"; c.drlevel=d; push(c);
    }

    // ── ProbMinHash ───────────────────────────────────────────────────────────
    // Default: kmer=21, size=1024
    {
        Combo c; c.algo="probmh"; c.kmer=21; c.size=1024; push(c);
    }
    for (int k : KMERS) {
        if (k == 21) continue;
        Combo c; c.algo="probmh"; c.kmer=k; c.size=1024; push(c);
    }
    for (int s : {64,96,128,192,256,384,512,768,1536,2048,3072,4096}) {
        Combo c; c.algo="probmh"; c.kmer=21; c.size=s; push(c);
    }

    return combos;
}

// Supplemental combos for partial re-runs (FASTA regenerated, one new param point).
static std::vector<Combo> make_supplement_combos(const std::string& kind)
{
    std::vector<Combo> combos;
    if (kind == "bindash-8k" || kind == "bindash_8k") {
        Combo c; c.algo = "bindash"; c.kmer = 21; c.bd_bbits = 16; c.bd_s64 = 64;
        combos.push_back(c);
    }
    return combos;
}

// ─── Evaluate one combo on one pair → return estimated Jaccard ─────────────────

static double eval_one_pair(const Combo& co,
                             const std::string& sa, const std::string& sb,
                             std::shared_ptr<const Sketch::kssd_parameter_t> kp)
{
    int La = (int)sa.size();
    int Lb = (int)sb.size();
    const auto& algo = co.algo;

    if (algo == "minhash") {
        Sketch::MinHash m1(co.kmer, co.size, 42, true),
                        m2(co.kmer, co.size, 42, true);
        m1.update(const_cast<char*>(sa.c_str()));
        m2.update(const_cast<char*>(sb.c_str()));
        m1.finalize(); m2.finalize();
        return m1.jaccard(&m2);
    }
    if (algo == "fastkmv") {
        Sketch::FastKMV v1((uint32_t)co.size, co.kmer, 42),
                        v2((uint32_t)co.size, co.kmer, 42);
        v1.update(const_cast<char*>(sa.c_str()), La);
        v2.update(const_cast<char*>(sb.c_str()), Lb);
        return v1.jaccard(v2);
    }
    if (algo == "hll") {
        Sketch::HyperLogLog h1(co.bits, co.kmer), h2(co.bits, co.kmer);
        h1.update(const_cast<char*>(sa.c_str()));
        h2.update(const_cast<char*>(sb.c_str()));
        return h1.jaccard_index(h2);
    }
    if (algo == "setsketch") {
        // For head-to-head comparison only: keep the legacy inclusion-exclusion
        // estimator under the "setsketch" label. The default jaccard_index()
        // now resolves to MLE, so we go through the explicit alias.
        Sketch::SetSketch s1(co.bits, co.ss_base, co.ss_a, co.kmer);
        Sketch::SetSketch s2(co.bits, co.ss_base, co.ss_a, co.kmer);
        s1.update(const_cast<char*>(sa.c_str()));
        s2.update(const_cast<char*>(sb.c_str()));
        return s1.jaccard_index_inclexcl(s2);
    }
    if (algo == "setsketch_mle") {
        // Ertl 2021 joint-MLE estimator (= new default)
        Sketch::SetSketch s1(co.bits, co.ss_base, co.ss_a, co.kmer);
        Sketch::SetSketch s2(co.bits, co.ss_base, co.ss_a, co.kmer);
        s1.update(const_cast<char*>(sa.c_str()));
        s2.update(const_cast<char*>(sb.c_str()));
        return s1.jaccard_index_mle(s2);
    }
    if (algo == "bindash") {
        Sketch::BinDash b1((uint32_t)co.bd_s64, co.kmer, (uint32_t)co.bd_bbits, 42);
        Sketch::BinDash b2((uint32_t)co.bd_s64, co.kmer, (uint32_t)co.bd_bbits, 42);
        b1.update(const_cast<char*>(sa.c_str()), (size_t)La);
        b2.update(const_cast<char*>(sb.c_str()), (size_t)Lb);
        return b1.jaccard(b2);
    }
    if (algo == "kssd") {
        Sketch::Kssd* k1 = new Sketch::Kssd(kp);
        Sketch::Kssd* k2 = new Sketch::Kssd(kp);
        k1->update(const_cast<char*>(sa.c_str()));
        k2->update(const_cast<char*>(sb.c_str()));
        double r = k1->jaccard(k2);
        delete k1; delete k2;
        return r;
    }
    if (algo == "probmh") {
        Sketch::ProbMinHash4 p1((uint32_t)co.size, co.kmer, 42),
                             p2((uint32_t)co.size, co.kmer, 42);
        p1.update(const_cast<char*>(sa.c_str()), La);
        p2.update(const_cast<char*>(sb.c_str()), Lb);
        return p1.jaccard(p2);
    }
    return -1.0;
}

// ─── CLI helpers ───────────────────────────────────────────────────────────────

static const char* get_arg(int argc, char** argv, const char* flag, const char* def)
{
    for (int i = 1; i + 1 < argc; i++)
        if (!strcmp(argv[i], flag)) return argv[i + 1];
    return def;
}
static int         geti(int a,char**v,const char*f,int    d){ auto s=get_arg(a,v,f,nullptr); return s?atoi(s):d; }
static std::string gets(int a,char**v,const char*f,const char*d){
    auto s=get_arg(a,v,f,d); return s ? std::string(s) : std::string(d ? d : "");
}
static bool hasf(int a,char**v,const char*f){
    for(int i=1;i<a;i++) if(!strcmp(v[i],f)) return true; return false;
}

// ─── main ──────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[])
{
    if (argc < 3 || hasf(argc, argv, "--help") || hasf(argc, argv, "-h")) {
        fprintf(stderr,
            "Usage: exe_eval_alisim_sweep --dir PATH [--threads T] [--model MODEL]\n"
            "       [--algos LIST] [--supplement KIND]\n"
            "Sweeps all algorithm × parameter combos on alisim FASTA pairs.\n"
            "  --supplement bindash-8k   only BinDash k=21 bbits=16 s64=64 (~8 KB)\n"
            "Outputs ALISIM_SWEEP lines to stdout.\n");
        return 1;
    }

    const std::string dir     = gets(argc, argv, "--dir",    ".");
    const std::string model   = gets(argc, argv, "--model",  "unknown");
    const std::string supplement = gets(argc, argv, "--supplement", "");
    int threads = geti(argc, argv, "--threads", omp_get_max_threads());
    if (threads < 1) threads = 1;

    // --algos filter: comma-separated whitelist (e.g. "setsketch,setsketch_mle").
    // Empty / missing → all algorithms.
    std::set<std::string> algo_whitelist;
    {
        const char* s = get_arg(argc, argv, "--algos", nullptr);
        if (s) {
            std::string cur;
            for (const char* p = s; ; ++p) {
                if (*p == ',' || *p == '\0') {
                    if (!cur.empty()) algo_whitelist.insert(cur);
                    cur.clear();
                    if (*p == '\0') break;
                } else cur += *p;
            }
            fprintf(stderr, "# --algos filter active (%zu algos):", algo_whitelist.size());
            for (auto& a : algo_whitelist) fprintf(stderr, " %s", a.c_str());
            fprintf(stderr, "\n");
        }
    }

    // ── read pairs_meta.tsv ──────────────────────────────────────────────────
    // Format: pair_id  branch_len  len_A  len_B
    std::string meta_path = dir + "/pairs_meta.tsv";
    std::vector<double> branch_lens;
    {
        std::ifstream mf(meta_path);
        if (!mf) { fprintf(stderr, "Cannot open %s\n", meta_path.c_str()); return 1; }
        std::string line;
        std::getline(mf, line); // skip header
        while (std::getline(mf, line)) {
            if (line.empty()) continue;
            double bl = 0.0;
            int    id = 0;
            sscanf(line.c_str(), "%d %lf", &id, &bl);
            branch_lens.push_back(bl);
        }
    }
    int N = (int)branch_lens.size();
    if (N == 0) { fprintf(stderr, "No pairs in %s\n", meta_path.c_str()); return 1; }

    fprintf(stderr, "# dir=%s model=%s N=%d threads=%d\n",
            dir.c_str(), model.c_str(), N, threads);

    // ── load all sequences ────────────────────────────────────────────────────
    std::vector<std::string> seqs_a(N), seqs_b(N);

    double t0 = get_sec();
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int i = 0; i < N; i++) {
        char buf[4096];
        snprintf(buf, sizeof(buf), "%s/pair_%d_A.fa", dir.c_str(), i);
        seqs_a[i] = read_fasta_seq(buf);
        snprintf(buf, sizeof(buf), "%s/pair_%d_B.fa", dir.c_str(), i);
        seqs_b[i] = read_fasta_seq(buf);
    }
    fprintf(stderr, "# loaded seqs: %.1f s\n", get_sec() - t0);

    // Check at least one pair loaded
    {
        int loaded = 0;
        for (int i = 0; i < N; i++)
            if (!seqs_a[i].empty() && !seqs_b[i].empty()) loaded++;
        if (loaded == 0) {
            fprintf(stderr, "No FASTA files found in %s\n", dir.c_str());
            return 1;
        }
        if (loaded < N)
            fprintf(stderr, "# Warning: only %d/%d pairs have FASTA files\n", loaded, N);
        N = loaded;
    }

    // ── collect unique k-mer sizes needed by all combos ───────────────────────
    // IMPORTANT: include KSSD's effective k values (kp->half_k * 2) here so the
    // sorted unique_ks vector — and thus the row layout of exact_j — is
    // determined ONCE before any precomputation. Adding rows later via
    // vector::resize() corrupts row→k alignment (resize only appends; existing
    // rows stay in place even though the sorted unique_ks shifts them).
    auto combos = supplement.empty() ? make_all_combos() : make_supplement_combos(supplement);
    if (!supplement.empty() && combos.empty()) {
        fprintf(stderr, "Unknown --supplement kind: %s\n", supplement.c_str());
        return 1;
    }
    // Apply algo whitelist filter (if any).
    if (!algo_whitelist.empty()) {
        std::vector<Combo> filtered;
        filtered.reserve(combos.size());
        for (auto& co : combos)
            if (algo_whitelist.count(co.algo)) filtered.push_back(co);
        combos.swap(filtered);
        fprintf(stderr, "# combos after algo filter: %zu\n", combos.size());
    }
    std::vector<int> unique_ks;

    // 1) Hash-based / SetSketch / etc. k values: directly from each combo
    for (auto& co : combos) {
        if (co.algo == "kssd") continue;
        if (std::find(unique_ks.begin(), unique_ks.end(), co.kmer) == unique_ks.end())
            unique_ks.push_back(co.kmer);
    }
    // 2) KSSD effective k values (half_k * 2): pre-build params now and add k
    std::map<int, std::shared_ptr<const Sketch::kssd_parameter_t>> kssd_map;
    for (auto& co : combos) {
        if (co.algo != "kssd") continue;
        if (kssd_map.count(co.drlevel)) continue;
        auto kp = get_kssd_params(10, 6, co.drlevel);
        kssd_map[co.drlevel] = kp;
        if (!kp) continue;
        int eff_k = kp->half_k * 2;
        if (std::find(unique_ks.begin(), unique_ks.end(), eff_k) == unique_ks.end())
            unique_ks.push_back(eff_k);
    }
    std::sort(unique_ks.begin(), unique_ks.end());

    fprintf(stderr, "# %zu combos, %zu unique k values: ",
            combos.size(), unique_ks.size());
    for (int k : unique_ks) fprintf(stderr, "%d ", k);
    fprintf(stderr, "\n");

    // ── precompute exact k-mer Jaccard for each (pair, k) ─────────────────────
    // exact_j[ki][pair] = exact k-mer Jaccard for unique_ks[ki] at pair index
    const int Nk = (int)unique_ks.size();
    std::vector<std::vector<double>> exact_j(Nk, std::vector<double>(N, -1.0));

    double t_exj0 = get_sec();
    for (int ki = 0; ki < Nk; ki++) {
        int k = unique_ks[ki];
        #pragma omp parallel for schedule(dynamic,4) num_threads(threads)
        for (int i = 0; i < N; i++) {
            if (seqs_a[i].empty() || seqs_b[i].empty()) continue;
            exact_j[ki][i] = exact_kmer_jaccard(seqs_a[i], seqs_b[i], k);
        }
        fprintf(stderr, "# exact_J k=%d: %.1f s\n", k, get_sec() - t_exj0);
        t_exj0 = get_sec();
    }

    // ── helper: look up exact Jaccard ─────────────────────────────────────────
    // Returns exact_j for the given k; for KSSD, kp->half_k*2 is used.
    auto get_exact_j = [&](const Combo& co, int pair_idx,
                            const Sketch::kssd_parameter_t* kp) -> double {
        int eff_k = (co.algo == "kssd" && kp) ? kp->half_k * 2 : co.kmer;
        auto it = std::find(unique_ks.begin(), unique_ks.end(), eff_k);
        if (it == unique_ks.end()) return -1.0;
        int ki = (int)(it - unique_ks.begin());
        return exact_j[ki][pair_idx];
    };

    // (kssd_map was already populated above when collecting unique k values.)

    // ── sweep all combos ──────────────────────────────────────────────────────
    double branch_len = branch_lens.empty() ? 0.0 : branch_lens[0]; // all pairs same BL

    int combo_idx = 0;
    for (const auto& co : combos) {
        combo_idx++;

        // get kssd params if needed
        std::shared_ptr<const Sketch::kssd_parameter_t> kp;
        if (co.algo == "kssd") {
            auto it = kssd_map.find(co.drlevel);
            if (it != kssd_map.end()) kp = it->second;
            if (!kp) {
                fprintf(stderr, "# KSSD drlevel=%d: failed to build params, skipping\n",
                        co.drlevel);
                continue;
            }
        }
        const Sketch::kssd_parameter_t* kp_raw = kp.get();

        // per-thread accumulators
        std::vector<double> tl_ae(threads, 0.0);
        std::vector<double> tl_se(threads, 0.0);
        std::vector<double> tl_exact(threads, 0.0);
        std::vector<double> tl_est(threads, 0.0);
        std::vector<int>    tl_cnt(threads, 0);

        #pragma omp parallel for schedule(dynamic,2) num_threads(threads)
        for (int i = 0; i < N; i++) {
            if (seqs_a[i].empty() || seqs_b[i].empty()) continue;

            int tid = omp_get_thread_num();

            double ej = get_exact_j(co, i, kp_raw);
            if (ej < 0.0) continue; // sequence too short for this k

            double est = eval_one_pair(co, seqs_a[i], seqs_b[i], kp);
            if (est < 0.0) est = 0.0;

            double err = est - ej;
            tl_ae[tid]    += std::fabs(err);
            tl_se[tid]    += err * err;
            tl_exact[tid] += ej;
            tl_est[tid]   += est;
            tl_cnt[tid]   += 1;
        }

        double ae = 0, se = 0, sum_ej = 0, sum_est = 0;
        int    n  = 0;
        for (int t = 0; t < threads; t++) {
            ae      += tl_ae[t];
            se      += tl_se[t];
            sum_ej  += tl_exact[t];
            sum_est += tl_est[t];
            n       += tl_cnt[t];
        }
        if (n == 0) continue;

        double mae  = ae / n;
        double rmse = std::sqrt(se / n);
        double mej  = sum_ej  / n;
        double mest = sum_est / n;

        // Print ALISIM_SWEEP line
        printf("ALISIM_SWEEP dir=%s model=%s branch_len=%.6f "
               "algo=%s kmer=%d size=%d bits=%d "
               "ss_a=%.1f ss_base=%.2f bd_bbits=%d bd_s64=%d drlevel=%d "
               "n=%d mean_exact_J=%.6f mean_est_J=%.6f MAE=%.6f RMSE=%.6f\n",
               dir.c_str(), model.c_str(), branch_len,
               co.algo.c_str(), co.kmer, co.size, co.bits,
               co.ss_a, co.ss_base, co.bd_bbits, co.bd_s64, co.drlevel,
               n, mej, mest, mae, rmse);
        fflush(stdout);

        if (combo_idx % 20 == 0)
            fprintf(stderr, "# combo %d/%zu done\n", combo_idx, combos.size());
    }

    fprintf(stderr, "# all %zu combos done, total %.1f s\n",
            combos.size(), get_sec() - t0);
    return 0;
}
