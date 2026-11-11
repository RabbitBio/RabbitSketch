/**
 * eval_reads – Coverage-robustness benchmark across all 7 RabbitSketch algorithms.
 *
 * Ground-truth metrics:
 *   exact_BJ  = |A∩B| / |A∪B|               (binary Jaccard, k-mer presence/absence)
 *   exact_J_P = Σ min(pA,pB), Σp=1           (Probability Jaccard, Ertl 2020)
 *               coverage-robust: same genome at any depth → J_P ≈ constant
 *
 * Algorithms tested (all BJ estimators except PMH_norm which estimates J_P):
 *   MinHash    (k, skSize)           – estimates BJ via bottom-k hashes
 *   FastKMV    (k, skSize)           – estimates BJ via KMV sketch
 *   HyperLogLog(hllBits=13)          – estimates BJ via HLL Jaccard
 *   SetSketch  (ssBits=13,b=2,a=5)  – estimates BJ via SetSketch
 *   BinDash    (sk64=16,k,bb=16)     – estimates BJ via b-bit one-perm hashing
 *   KSSD       (half_k=k/2,hsub=6,dr=3) – estimates BJ via KSSD distance
 *   PMH_norm   (k, skSize)           – estimates J_P via weighted ProbMinHash ★
 *                                      weight = count/total_count (coverage-robust)
 *
 * ── Mode 1: reads-vs-reference (recommended for application evaluation) ──────
 *
 *   The biologically meaningful experiment: query reads at various coverage depths
 *   are searched against a reference database of assembled genomes.  For each
 *   query the rank of the correct reference is reported for every algorithm.
 *
 *   Binary-Jaccard estimators (MinHash, FastKMV, …) lose similarity as coverage
 *   drops because missing k-mers reduce intersection.  PMH_norm (J_P) normalises
 *   by frequency and retains more signal at low coverage.
 *
 *   Usage:
 *     eval_reads --mode ref \
 *                --db-list    db.list    \   # assembled reference genomes (one per line)
 *                --query-list query.list \   # reads files at various depths (one per line)
 *               [--truth      truth.tsv] \   # tab: query_basename <TAB> ref_basename
 *               [--k 21] [--sketch-size 1024] [--threads 8]
 *
 *   If --truth is omitted, the correct reference for a query is auto-detected:
 *   the db entry whose basename is the longest prefix of the query basename.
 *   Example: query "saureus_5x.fastq.gz" → ref "saureus.fna" (prefix "saureus").
 *
 *   Output (stdout):
 *     TSV with one row per (query × method):
 *       query  true_ref  depth  method  rank  sim_true  sim_best  best_ref  correct
 *     Followed by a Recall@1 summary table.
 *
 * ── Mode 2: all-pairs (original, backward-compatible) ────────────────────────
 *
 *   Usage:
 *     eval_reads [--mode pairs] <file_list> [k=21] [sketch_size=1024] [threads=8]
 */

#include "Sketch.h"
#include "probmh.h"
#include "fastkmv.h"
#include "SetSketch.h"
#include "BinDash.h"
#include "kseq.h"
#include "phmap.h"
#include "wyhash_v4_tmh.h"
#include "common.h"

#include <zlib.h>
#include <omp.h>
#include <cassert>
#include <mutex>
#include <memory>
#include <tuple>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <map>

KSEQ_INIT(gzFile, gzread)
using namespace std;

// ── 2-bit canonical k-mer helpers ─────────────────────────────────────────
// Only k ≤ 21 (fits in uint64_t with 2 bits per base, 42 bits total).
static const uint8_t BASE2[256] = {
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
static inline uint8_t comp2(uint8_t b) { return 3 ^ b; }

// Entropy weight matching ProbMinHash4::updateEntropy (escale=0.5, w_min=0.1).
// Input: canonical k-mer as 2-bit uint64_t with k bases packed in low bits.
static double kmer_entropy_weight(uint64_t can, int k, double w_min = 0.1) {
    int bcnt[4] = {0, 0, 0, 0};
    uint64_t tmp = can;
    for (int i = 0; i < k; ++i) {
        bcnt[tmp & 3]++;
        tmp >>= 2;
    }
    double log2k = log2((double)k);
    double inv_k = 1.0 / (double)k;
    double nlogn = 0.0;
    for (int b = 0; b < 4; ++b)
        if (bcnt[b] > 0)
            nlogn += (double)bcnt[b] * log2((double)bcnt[b]);
    double H  = log2k - inv_k * nlogn;  // Shannon entropy in bits ∈ [0, log2(4)=2]
    double wt = H * 0.5;                // scale to [0, 1]
    return (wt < w_min) ? w_min : wt;
}

// ── KSSD parameter cache (thread-safe) ────────────────────────────────────
static std::shared_ptr<const Sketch::kssd_parameter_t>
get_kssd_params(int half_k, int half_subk, int drlevel) {
    return Sketch::getKssdParameters(half_k, half_subk, drlevel);
}

// ── Per-file data ──────────────────────────────────────────────────────────
struct FileData {
    string                                   path;
    phmap::flat_hash_map<uint64_t,uint32_t>  freq;   // 2-bit canonical → count
    uint64_t total_count = 0;

    // BJ estimators
    Sketch::MinHash*        minhash  = nullptr;
    Sketch::FastKMV*        fastkmv  = nullptr;
    Sketch::HyperLogLog*    hll      = nullptr;
    Sketch::SetSketch*      ssk      = nullptr;
    Sketch::BinDash*        bindash  = nullptr;
    Sketch::Kssd*           kssd     = nullptr;

    // J_P estimator (coverage-robust) ★
    Sketch::ProbMinHash4*   pmh_norm = nullptr;

    ~FileData() {
        delete minhash; delete fastkmv; delete hll;
        delete ssk; delete bindash; delete kssd; delete pmh_norm;
    }
};

// ── Pass 1: canonical k-mer counting ──────────────────────────────────────
// Stores the 2-bit canonical k-mer itself as the map key (not its hash),
// so that we can decode its entropy later for the ground truth computation.
static void count_kmers(FileData& fd, int k) {
    const uint64_t mask     = (1ULL << (2*k)) - 1;
    const uint64_t hi_shift = 2ULL * (k - 1);

    gzFile fp = gzopen(fd.path.c_str(), "r");
    if (!fp) { cerr << "ERROR: cannot open " << fd.path << "\n"; return; }
    kseq_t* ks = kseq_init(fp);

    while (kseq_read(ks) >= 0) {
        const char* seq = ks->seq.s;
        const int   len = (int)ks->seq.l;
        if (len < k) continue;

        uint64_t fwd = 0, rev = 0;
        int valid = 0;

        for (int i = 0; i < len; ++i) {
            uint8_t b = BASE2[(uint8_t)seq[i]];
            if (b == 0xFF) { valid = 0; fwd = 0; rev = 0; continue; }
            fwd = ((fwd << 2) | b) & mask;
            rev = (rev >> 2) | ((uint64_t)comp2(b) << hi_shift);
            if (++valid >= k) {
                uint64_t can = (fwd <= rev) ? fwd : rev;
                ++fd.freq[can];
                ++fd.total_count;
            }
        }
    }
    kseq_destroy(ks); gzclose(fp);
}

// ── Pass 2: build all sketches ─────────────────────────────────────────────
static void build_sketches(FileData& fd, int k, int skSize) {
    // BJ estimators
    fd.minhash  = new Sketch::MinHash(k, skSize, 42, /*rc=*/true);
    fd.fastkmv  = new Sketch::FastKMV(static_cast<uint32_t>(skSize), k, 42ULL);
    fd.hll      = new Sketch::HyperLogLog(13);        // 2^13 = 8192 registers
    fd.ssk      = new Sketch::SetSketch(13, 2.0, 5.0);
    // BinDash: sk64*64 = total bins; sk64=skSize/64 to match skSize bins
    const uint32_t sk64 = static_cast<uint32_t>(skSize) / 64;
    fd.bindash  = new Sketch::BinDash(sk64 > 0 ? sk64 : 1, k, 16, 42ULL);
    fd.kssd     = new Sketch::Kssd(get_kssd_params(k/2, 6, 3));

    // PMH-norm: weight = count/total_count (coverage-robust J_P estimator)
    fd.pmh_norm = new Sketch::ProbMinHash4(static_cast<uint32_t>(skSize), k);
    const double inv_total = (fd.total_count > 0) ? 1.0 / (double)fd.total_count : 1.0;
    for (const auto& [can, cnt] : fd.freq)
        fd.pmh_norm->addHash(wyhash(&can, 8, 0), (double)cnt * inv_total);

    gzFile fp = gzopen(fd.path.c_str(), "r");
    if (!fp) return;
    kseq_t* ks = kseq_init(fp);

    while (kseq_read(ks) >= 0) {
        char*     seq = ks->seq.s;
        const int len = (int)ks->seq.l;
        if (len <= 0) continue;

        fd.minhash->update(seq);
        fd.fastkmv->update(seq, (uint64_t)len);
        fd.hll->update(seq);
        fd.ssk->update(seq);
        fd.bindash->update(seq, (uint64_t)len);
        fd.kssd->update(seq, static_cast<size_t>(len));
    }
    kseq_destroy(ks); gzclose(fp);

    fd.kssd->finalize();
    fd.minhash->finalize();
    fd.bindash->finalize();
}

// ── Exact BJ, WJ, EJ from freq tables ────────────────────────────────────
struct GroundTruth { double bj, wj, wj_norm, jp, ej, containment_ms; };

// Ground truth definitions:
//   BJ      = |A∩B| / |A∪B|                    (binary Jaccard)
//   WJ      = Σmin(fA,fB) / Σmax(fA,fB)        (freq-weighted Jaccard)
//   WJ_norm = Σmin(pA,pB) / Σmax(pA,pB)        (relative-freq WJ)
//   J_P     = Σmin(pA,pB)  with Σp=1           (Probability Jaccard, Ertl 2020)
//                                                 = ProbMinHash collision probability
//   EJ      = Σ_{x∈A∩B} H(x) / Σ_{x∈A∪B} H(x) (entropy-weighted Jaccard)
//
// J_P and WJ_norm are related by: J_P = 2·WJ_norm / (1 + WJ_norm)
// J_P is what ProbMinHash with normalized weights actually estimates.
static GroundTruth exact_jaccard(
    const phmap::flat_hash_map<uint64_t,uint32_t>& A, uint64_t totA,
    const phmap::flat_hash_map<uint64_t,uint32_t>& B, uint64_t totB,
    int k)
{
    uint64_t bj_inter = 0, bj_union = 0;
    double   wj_num = 0.0, wj_den = 0.0;
    double   wn_num = 0.0, wn_den = 0.0;
    double   jp_sum = 0.0;
    double   ej_num = 0.0, ej_den = 0.0;

    const double invA = (totA > 0) ? 1.0/totA : 1.0;
    const double invB = (totB > 0) ? 1.0/totB : 1.0;

    for (const auto& [can, fa] : A) {
        auto it = B.find(can);
        uint32_t fb = (it != B.end()) ? it->second : 0;
        bool in_B = (fb > 0);

        bj_inter += in_B ? 1 : 0;
        bj_union += 1;

        wj_num += (double)min(fa, fb);
        wj_den += (double)max(fa, fb);

        double pa = fa * invA, pb = fb * invB;
        wn_num += min(pa, pb);
        wn_den += max(pa, pb);
        jp_sum += min(pa, pb);  // J_P contribution (only nonzero if in both)

        double ew = kmer_entropy_weight(can, k);
        ej_num += in_B ? ew : 0.0;
        ej_den += ew;
    }
    for (const auto& [can, fb] : B) {
        if (A.find(can) == A.end()) {
            bj_union += 1;
            wj_den   += (double)fb;
            wn_den   += fb * invB;
            ej_den   += kmer_entropy_weight(can, k);
        }
    }

    const double c_ms = (totA > 0) ? wj_num / (double)totA : 0.0;
    return {
        (bj_union > 0) ? (double)bj_inter / bj_union : 0.0,
        (wj_den   > 0) ? wj_num / wj_den             : 0.0,
        (wn_den   > 0) ? wn_num / wn_den             : 0.0,
        jp_sum,
        (ej_den   > 0) ? ej_num / ej_den             : 0.0,
        c_ms
    };
}

// ── Helpers ────────────────────────────────────────────────────────────────

static auto shortname = [](const string& p) -> string {
    size_t s = p.rfind('/'); return (s == string::npos) ? p : p.substr(s+1);
};

static vector<string> read_list(const string& path) {
    vector<string> v;
    ifstream fs(path); string ln;
    while (getline(fs, ln)) { if (!ln.empty()) v.push_back(ln); }
    return v;
}

// Strip common extensions to get a "species prefix" for truth auto-detection.
// e.g. "saureus_5x.fastq.gz" → "saureus_5x", "saureus.fna" → "saureus"
static string strip_ext(const string& name) {
    string s = name;
    for (const char* ext : {".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz",
                             ".fastq", ".fq", ".fasta", ".fa", ".fna", ".gz"}) {
        size_t n = strlen(ext);
        if (s.size() > n && s.compare(s.size()-n, n, ext) == 0)
            s = s.substr(0, s.size()-n);
    }
    return s;
}

// Auto-detect truth: for a query basename, find the db ref whose stripped
// basename is the longest prefix of the query's stripped basename.
// e.g. query="saureus_5x" matches ref="saureus" (longest common prefix match).
static int auto_detect_truth(const string& query_short,
                              const vector<string>& db_shorts) {
    string qbase = strip_ext(query_short);
    int best_idx = -1;
    size_t best_len = 0;
    for (int i = 0; i < (int)db_shorts.size(); ++i) {
        string rbase = strip_ext(db_shorts[i]);
        // Check if qbase starts with rbase
        if (qbase.size() >= rbase.size() &&
            qbase.compare(0, rbase.size(), rbase) == 0 &&
            rbase.size() > best_len) {
            best_len = rbase.size();
            best_idx = i;
        }
    }
    return best_idx;
}

// Extract coverage depth string from filename, e.g. "saureus_5x" → "5x".
// Falls back to the full stripped basename if no pattern found.
static string extract_depth(const string& name) {
    string s = strip_ext(shortname(name));
    // Find last underscore-delimited token that looks like a depth: [0-9.]+x
    size_t pos = s.rfind('_');
    while (pos != string::npos) {
        string tok = s.substr(pos + 1);
        bool looks_like_depth = !tok.empty() &&
            tok.back() == 'x' &&
            tok.find_first_not_of("0123456789.x") == string::npos;
        if (looks_like_depth) return tok;
        pos = (pos == 0) ? string::npos : s.rfind('_', pos - 1);
    }
    // For ecoli-style names like "ecoli_30pct", return last segment
    pos = s.rfind('_');
    return (pos != string::npos) ? s.substr(pos + 1) : s;
}

// Build a vector of FileData objects from a file list (parallel).
static vector<FileData*> build_all(const vector<string>& files,
                                   int k, int skSize, int nThreads,
                                   const string& tag) {
    int N = (int)files.size();
    vector<FileData*> data(N);
    for (int i = 0; i < N; ++i) { data[i] = new FileData(); data[i]->path = files[i]; }

    cerr << "  [" << tag << "] Pass 1: k-mer counting (" << N << " files)\n";
    #pragma omp parallel for num_threads(nThreads) schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        count_kmers(*data[i], k);
        #pragma omp critical
        cerr << "    [pass1] " << shortname(files[i])
             << ": " << data[i]->freq.size() << " unique " << k << "-mers\n";
    }

    cerr << "  [" << tag << "] Pass 2: sketch building\n";
    #pragma omp parallel for num_threads(nThreads) schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        build_sketches(*data[i], k, skSize);
        #pragma omp critical
        cerr << "    [pass2] " << shortname(files[i]) << " done\n";
    }
    return data;
}

// ── Method names (7 algorithms) ────────────────────────────────────────────
static const vector<string> METHOD_NAMES = {
    "MinHash", "FastKMV", "HLL", "SetSketch", "BinDash", "KSSD", "PMH_norm"
};

static const vector<string> CONTAINMENT_METHOD_NAMES = {
    "MinHash", "FastKMV", "HLL", "SetSketch", "BinDash", "KSSD"
};

// C(query ⊆ ref) from J and sketch cardinalities / list sizes (KSSD has no containment()).
static double containment_from_jaccard(double j, double size_query, double size_ref) {
    if (size_query <= 0.0 || j <= 0.0) return 0.0;
    return j * (size_query + size_ref) / (size_query * (1.0 + j));
}

static double kssd_list_size(const Sketch::Kssd* sk) {
    return sk->getHashes64().empty()
        ? (double)sk->getHashes().size()
        : (double)sk->getHashes64().size();
}

static double kssd_containment_query_in_ref(Sketch::Kssd* q, Sketch::Kssd* r) {
    const double j = q->jaccard(r);
    return containment_from_jaccard(j, kssd_list_size(q), kssd_list_size(r));
}

// Six binary-sketch estimators: C(query reads ⊆ reference genome).
static vector<double> compute_containment_sims(const FileData* q, const FileData* r) {
  return {
      r->minhash->containJaccard(q->minhash),
      q->fastkmv->containment(*r->fastkmv),
      q->hll->containment(*r->hll),
      q->ssk->containment(*r->ssk),
      q->bindash->containment(*r->bindash),
      kssd_containment_query_in_ref(q->kssd, r->kssd)
  };
}

// Compute all 7 similarity scores between query q and reference r.
static vector<double> compute_sims(const FileData* q, const FileData* r) {
    return {
        q->minhash ->jaccard(r->minhash),
        q->fastkmv ->jaccard(*r->fastkmv),
        q->hll     ->jaccard_index(*r->hll),
        q->ssk     ->jaccard_index(*r->ssk),
        q->bindash ->jaccard(*r->bindash),
        q->kssd    ->jaccard(r->kssd),
        q->pmh_norm->jaccard_weighted(*r->pmh_norm)
    };
}

// ── Mode 1: reads-vs-reference ─────────────────────────────────────────────
static void run_ref_mode(
    const vector<string>& db_files,
    const vector<string>& query_files,
    const string& truth_path,
    int k, int skSize, int nThreads,
    bool use_containment)
{
    double t0 = get_sec();

    const vector<string>& active_methods =
        use_containment ? CONTAINMENT_METHOD_NAMES : METHOD_NAMES;
    const int M = (int)active_methods.size();

    cerr << "=== eval_reads (ref mode"
         << (use_containment ? ", containment" : ", BJ/J_P")
         << "): " << db_files.size() << " refs, "
         << query_files.size() << " queries"
         << "  k=" << k << "  skSize=" << skSize
         << "  threads=" << nThreads << "\n";

    // Build sketches
    auto db_data    = build_all(db_files,    k, skSize, nThreads, "db");
    auto query_data = build_all(query_files, k, skSize, nThreads, "query");

    const int Ndb = (int)db_files.size();
    const int Nq  = (int)query_files.size();

    // Short names for display
    vector<string> db_short(Ndb), q_short(Nq);
    for (int i = 0; i < Ndb; ++i) db_short[i] = shortname(db_files[i]);
    for (int i = 0; i < Nq;  ++i) q_short[i]  = shortname(query_files[i]);

    // Load or auto-detect truth: truth[qi] = index into db_data of correct ref
    vector<int> truth(Nq, -1);
    if (!truth_path.empty()) {
        // Load explicit mapping: query_basename <TAB> ref_basename
        ifstream fs(truth_path); string ln;
        while (getline(fs, ln)) {
            size_t tab = ln.find('\t');
            if (tab == string::npos) continue;
            string qname = ln.substr(0, tab);
            string rname = ln.substr(tab + 1);
            for (int qi = 0; qi < Nq;  ++qi)
                if (q_short[qi] == qname)
                    for (int ri = 0; ri < Ndb; ++ri)
                        if (db_short[ri] == rname) { truth[qi] = ri; break; }
        }
    } else {
        for (int qi = 0; qi < Nq; ++qi)
            truth[qi] = auto_detect_truth(q_short[qi], db_short);
    }

    // Warn about unresolved truth
    for (int qi = 0; qi < Nq; ++qi)
        if (truth[qi] < 0)
            cerr << "  WARNING: cannot find truth ref for query: " << q_short[qi] << "\n";

    // ── Evaluate: for each query, rank all db refs per method ───────────────
    // results[qi][m] = {rank_of_true (1-based), sim_to_true, best_sim, best_ref_idx}
    struct Result { int rank; double sim_true; double sim_best; int best_ref; };
    vector<vector<Result>> results(Nq, vector<Result>(M));

    double t1 = get_sec();
    cerr << "  Computing similarities...\n";

    #pragma omp parallel for num_threads(nThreads) schedule(dynamic)
    for (int qi = 0; qi < Nq; ++qi) {
        // sim_matrix[ri][m] = similarity of query qi to ref ri for method m
        vector<vector<double>> sims(Ndb);
        for (int ri = 0; ri < Ndb; ++ri) {
            sims[ri] = use_containment
                ? compute_containment_sims(query_data[qi], db_data[ri])
                : compute_sims(query_data[qi], db_data[ri]);
        }

        for (int m = 0; m < M; ++m) {
            // Sort refs by descending similarity for method m
            vector<int> order(Ndb);
            iota(order.begin(), order.end(), 0);
            sort(order.begin(), order.end(),
                 [&](int a, int b){ return sims[a][m] > sims[b][m]; });

            int true_ri = truth[qi];
            double sim_true = (true_ri >= 0) ? sims[true_ri][m] : 0.0;
            int best_ref    = order[0];
            double sim_best = sims[best_ref][m];

            // Rank = position of correct ref in sorted list (1-based)
            int rank = Ndb; // default: last
            for (int pos = 0; pos < Ndb; ++pos)
                if (order[pos] == true_ri) { rank = pos + 1; break; }

            results[qi][m] = {rank, sim_true, sim_best, best_ref};
        }
    }
    double t2 = get_sec();
    cerr << "  Similarity computation: " << t2-t1 << " s\n";

    // Also compute exact ground-truth for true pairs (query vs. correct ref)
    // Output exact BJ and J_P for reference
    vector<GroundTruth> exact(Nq);
    for (int qi = 0; qi < Nq; ++qi) {
        int ri = truth[qi];
        if (ri < 0) continue;
        exact[qi] = exact_jaccard(query_data[qi]->freq, query_data[qi]->total_count,
                                  db_data[ri]->freq,    db_data[ri]->total_count, k);
    }

    // ── Output TSV ──────────────────────────────────────────────────────────
    // Header
    cout << "# eval_reads ref mode: reads query vs reference database\n";
    cout << "# k=" << k << "  skSize=" << skSize << "  Ndb=" << Ndb << "  Nq=" << Nq << "\n";
    if (use_containment)
        cout << "# Metric: containment C(query⊆ref).  exact_C = Σmin(fQ,fR)/ΣfQ.\n";
    else
        cout << "# Rank: 1 = best match (correct).  sim_true = similarity to correct ref.\n";
    cout << "#\n";
    if (use_containment) {
        cout << "query\ttrue_ref\tdepth\texact_C";
    } else {
        cout << "query\ttrue_ref\tdepth\texact_BJ\texact_JP";
    }
    for (const auto& m : active_methods)
        cout << "\t" << m << "_rank\t" << m << "_sim_true\t" << m << "_sim_best\t" << m << "_best_ref";
    cout << "\n";

    for (int qi = 0; qi < Nq; ++qi) {
        int ri = truth[qi];
        string true_name = (ri >= 0) ? db_short[ri] : "?";
        string depth = extract_depth(query_files[qi]);

        cout << q_short[qi] << "\t" << true_name << "\t" << depth << "\t";
        if (use_containment)
            cout << (ri >= 0 ? exact[qi].containment_ms : -1.0);
        else
            cout << (ri >= 0 ? exact[qi].bj : -1.0) << "\t"
                 << (ri >= 0 ? exact[qi].jp : -1.0);

        for (int m = 0; m < M; ++m) {
            const auto& r = results[qi][m];
            cout << "\t" << r.rank
                 << "\t" << r.sim_true
                 << "\t" << r.sim_best
                 << "\t" << (r.best_ref >= 0 ? db_short[r.best_ref] : "?");
        }
        cout << "\n";
    }

    // ── Recall@1 summary ────────────────────────────────────────────────────
    cout << "\n# ── Recall@1 summary (fraction of queries where correct ref ranked 1st) ──\n";
    cout << "depth";
    for (const auto& m : active_methods) cout << "\t" << m;
    cout << "\n";

    // Group queries by depth
    map<string, vector<int>> by_depth;
    for (int qi = 0; qi < Nq; ++qi)
        by_depth[extract_depth(query_files[qi])].push_back(qi);

    // Sort depth keys naturally (0.5x < 1x < 2x < 5x < 10x < 20x < 40x)
    auto depth_val = [](const string& d) -> double {
        double v = 0.0;
        try { v = stod(d); } catch (...) {}
        return v;
    };
    vector<string> depth_keys;
    for (auto& [d, _] : by_depth) depth_keys.push_back(d);
    sort(depth_keys.begin(), depth_keys.end(),
         [&](const string& a, const string& b){ return depth_val(a) < depth_val(b); });

    for (const auto& dep : depth_keys) {
        const auto& idxs = by_depth[dep];
        int n = (int)idxs.size();
        cout << dep;
        for (int m = 0; m < M; ++m) {
            int correct = 0;
            for (int qi : idxs)
                if (results[qi][m].rank == 1) ++correct;
            cout << "\t" << fixed << setprecision(3)
                 << (n > 0 ? (double)correct / n : 0.0);
        }
        cout << "\n";
    }
    // Overall row
    cout << "ALL";
    for (int m = 0; m < M; ++m) {
        int correct = 0, total = 0;
        for (int qi = 0; qi < Nq; ++qi) {
            if (truth[qi] < 0) continue;
            ++total;
            if (results[qi][m].rank == 1) ++correct;
        }
        cout << "\t" << fixed << setprecision(3)
             << (total > 0 ? (double)correct / total : 0.0);
    }
    cout << "\n";

    cerr << "Total: " << get_sec()-t0 << " s\n";

    for (auto* d : db_data)    delete d;
    for (auto* d : query_data) delete d;
}

// ── Mode 2: all-pairs (original) ───────────────────────────────────────────
static void run_pairs_mode(const vector<string>& files,
                            int k, int skSize, int nThreads)
{
    const int N = (int)files.size();
    cerr << "=== eval_reads (pairs mode): N=" << N << "  k=" << k
         << "  skSize=" << skSize << "  threads=" << nThreads << "\n";

    auto data = build_all(files, k, skSize, nThreads, "pairs");

    cout << "# All 7 RabbitSketch algorithms vs exact_BJ and exact_J_P ground truth\n"
         << "# Algorithms 1-6 estimate Binary Jaccard (BJ); PMH_norm estimates J_P\n"
         << "#\n";
    cout << "fileA\tfileB\t"
            "exact_BJ\texact_J_P\t"
            "MinHash\tFastKMV\tHLL\tSetSketch\tBinDash\tKSSD\t"
            "PMH_norm\n";

    double t2 = get_sec();
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {
            auto gt = exact_jaccard(data[i]->freq, data[i]->total_count,
                                    data[j]->freq, data[j]->total_count, k);
            auto sims = compute_sims(data[i], data[j]);

            printf("%-35s\t%-35s\t%.4f\t%.4f\t"
                   "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                shortname(files[i]).c_str(), shortname(files[j]).c_str(),
                gt.bj, gt.jp,
                sims[0], sims[1], sims[2], sims[3], sims[4], sims[5], sims[6]);
        }
    }
    cerr << "Pairwise eval: " << get_sec()-t2 << " s\n";
    for (auto* d : data) delete d;
}

// ── main ───────────────────────────────────────────────────────────────────
int main(int argc, char* argv[])
{
    // Parse arguments
    string mode       = "pairs";  // default: backward-compatible
    string metric     = "bj";     // ref mode: bj | containment
    string db_list_path, query_list_path, truth_path;
    string pairs_list_path;
    int k        = 21;
    int skSize   = 1024;
    int nThreads = 8;

    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        if      (a == "--mode"         && i+1 < argc) { mode          = argv[++i]; }
        else if (a == "--metric"       && i+1 < argc) { metric        = argv[++i]; }
        else if (a == "--db-list"      && i+1 < argc) { db_list_path  = argv[++i]; mode = "ref"; }
        else if (a == "--query-list"   && i+1 < argc) { query_list_path = argv[++i]; mode = "ref"; }
        else if (a == "--truth"        && i+1 < argc) { truth_path    = argv[++i]; }
        else if ((a == "--k" || a == "-k") && i+1 < argc) { k        = atoi(argv[++i]); }
        else if (a == "--sketch-size"  && i+1 < argc) { skSize        = atoi(argv[++i]); }
        else if (a == "--threads"      && i+1 < argc) { nThreads      = atoi(argv[++i]); }
        else if (a[0] != '-' && pairs_list_path.empty()) { pairs_list_path = a; }
        else if (a[0] != '-' && k == 21  ) { k        = atoi(a.c_str()); }
        else if (a[0] != '-' && skSize == 1024) { skSize   = atoi(a.c_str()); }
        else if (a[0] != '-'             ) { nThreads  = atoi(a.c_str()); }
    }

    if (k < 1 || k > 21) { cerr << "k must be 1..21\n"; return 1; }

    if (mode == "ref") {
        if (db_list_path.empty() || query_list_path.empty()) {
            cerr << "Usage (ref mode): " << argv[0]
                 << " --db-list db.list --query-list queries.list"
                 << " [--truth truth.tsv] [--k 21] [--sketch-size 1024] [--threads 8]\n";
            return 1;
        }
        auto db_files    = read_list(db_list_path);
        auto query_files = read_list(query_list_path);
        if (db_files.empty() || query_files.empty()) {
            cerr << "ERROR: empty db or query list\n"; return 1;
        }
        const bool use_containment = (metric == "containment");
        if (metric != "bj" && metric != "jaccard" && !use_containment) {
            cerr << "ERROR: --metric must be bj or containment\n";
            return 1;
        }
        run_ref_mode(db_files, query_files, truth_path, k, skSize, nThreads,
                     use_containment);
    } else {
        // pairs mode (original)
        if (pairs_list_path.empty()) {
            cerr << "Usage: " << argv[0]
                 << " <file_list> [k=21] [sketch_size=1024] [threads=8]\n"
                 << "       " << argv[0]
                 << " --mode ref --db-list db.list --query-list q.list [options]\n";
            return 1;
        }
        auto files = read_list(pairs_list_path);
        if ((int)files.size() < 2) { cerr << "Need >= 2 files\n"; return 1; }
        run_pairs_mode(files, k, skSize, nThreads);
    }
    return 0;
}
