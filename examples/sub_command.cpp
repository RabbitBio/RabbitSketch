/**
 * sub_command.cpp – every rabbitsketch algorithm pipeline lives here.
 *
 * Layout
 * ──────
 *   §1  Helpers (file list IO, mash distance, seq streaming, KSSD param cache)
 *   §2  Per-algorithm sketch builders (build_*)        — pair / all-pairs mode
 *   §3  run_pair                                       — pairwise distance
 *   §4  list_allpairs<T>                               — generic O(N²)
 *   §5  run_index_fastkmv  (mirrors test_FastKMV.cpp)
 *   §6  run_index_probmh   (mirrors test_ProbMinHash.cpp)
 *   §7  run_index_setsketch (mirrors test_SetSketch.cpp inverted-index mode)
 *   §8  run_index_kssd     (mirrors test_Kssd.cpp – SIMD encode + 64-shard merge
 *                          + flat CSR index + stamp-based dist + 5-col output)
 *   §9  run_allpairs_bindash (mirrors test_BinDash.cpp – flat_usigs + sort by
 *                          cardinality + tile-based dist + per-thread temp files)
 *   §10 run_cli            — dispatch table
 */

#include "sub_command.h"

#include "Sketch.h"
#include "BinDash.h"
#include "fastkmv.h"
#include "probmh.h"
#include "SetSketch.h"
#include "InvertedIndex.h"
#include "shuffle.h"
#include "common.h"
#include "kseq.h"
#include "phmap.h"

#include <zlib.h>
#include <omp.h>
#include <err.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <signal.h>
#ifdef __x86_64__
#  include <immintrin.h>
#endif

#include <algorithm>
#include <atomic>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

KSEQ_INIT(gzFile, gzread)

// ═══════════════════════════════════════════════════════════════════════════
//  §1  Helpers
// ═══════════════════════════════════════════════════════════════════════════

static std::vector<std::string> load_file_list(const std::string& path) {
    std::ifstream fs(path);
    if (!fs) err(errno, "cannot open file list: %s", path.c_str());
    std::vector<std::string> out;
    std::string line;
    while (std::getline(fs, line))
        if (!line.empty()) out.push_back(line);
    return out;
}

// ── Frequency-weighted k-mer counting for PMH-norm ────────────────────────
// Estimates J_P = Σ min(p_A,p_B) / Σ max(p_A,p_B)  (coverage-robust metric).
//
// Key design choices vs naïve phmap approach:
//
//  1. Vector+sort instead of hash map for counting:
//     phmap::flat_hash_map grows by 2x with mmap(), causing hundreds of millions
//     of minor page faults per run (observed: 693M faults, 5179 s system time
//     for 208K files).  A std::vector<uint64_t> + std::sort uses one contiguous
//     allocation, making page faults rare and avoiding kernel overhead entirely.
//     Result is identical (addHash order does not affect the PMH sketch).
//
//  2. libdeflate for .gz decompression (2-3× faster than zlib):
//     Reads the full compressed file into a raw buffer, decompresses in one
//     shot, then parses FASTQ in-memory without zlib overhead.  Falls back
//     to kseq (zlib) for non-gz files or if decompression fails.
//
//  3. 1/PMH_SDIV hash-sampling keeps per-thread memory constant (~10 MB):
//     E[J_P_sampled] = J_P  (numerator and denominator both scale by 1/SDIV).
// ──────────────────────────────────────────────────────────────────────────

#include <libdeflate.h>
#include <sys/stat.h>
#include <fcntl.h>

static constexpr uint64_t PMH_SDIV      = 8;
static constexpr uint64_t PMH_STHRESH   = UINT64_MAX / PMH_SDIV;
static constexpr uint64_t PMH_SMIX      = 0x9E3779B97F4A7C15ULL;

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

// Process k-mers from a sequence buffer, accumulating sampled canonical k-mers
// into `sampled` and counting all k-mer instances in `total`.
static inline void pmh_kmers(const char* seq, int len, int k,
                              const uint64_t mask, const uint64_t hi_shift,
                              std::vector<uint64_t>& sampled, uint64_t& total) {
    uint64_t fwd = 0, rev = 0; int valid = 0;
    for (int i = 0; i < len; ++i) {
        const uint8_t b = PMH_BASE2[(uint8_t)seq[i]];
        if (b == 0xFF) { valid = 0; fwd = 0; rev = 0; continue; }
        fwd = ((fwd << 2) | b) & mask;
        rev = (rev >> 2) | ((uint64_t)(3 ^ b) << hi_shift);
        if (++valid >= k) {
            const uint64_t can = (fwd <= rev) ? fwd : rev;
            ++total;
            if ((can * PMH_SMIX) <= PMH_STHRESH)
                sampled.push_back(can);
        }
    }
}

// Fast in-memory FASTQ parser: only extracts sequence lines.
// Assumes single-line sequences per record (standard short-read FASTQ).
static void parse_fastq_mem(const char* buf, size_t len, int k,
                             const uint64_t mask, const uint64_t hi_shift,
                             std::vector<uint64_t>& sampled, uint64_t& total) {
    const char* p   = buf;
    const char* end = buf + len;
    while (p < end) {
        if (*p != '@') { ++p; continue; }
        // skip header line
        while (p < end && *p != '\n') ++p;
        if (p >= end) break; ++p;
        // sequence line
        const char* seq = p;
        while (p < end && *p != '\n') ++p;
        int slen = (int)(p - seq);
        if (slen >= k) pmh_kmers(seq, slen, k, mask, hi_shift, sampled, total);
        if (p >= end) break; ++p;
        // '+' line
        while (p < end && *p != '\n') ++p;
        if (p >= end) break; ++p;
        // quality line (skip)
        while (p < end && *p != '\n') ++p;
        if (p < end) ++p;
    }
}

// Build a PMH-norm sketch using libdeflate for .gz, kseq fallback otherwise.
// Uses vector+sort for counting (no hash map = no mmap page-fault explosion).
static Sketch::ProbMinHash4*
pmh_build_weighted(const std::string& path, uint32_t m, int k,
                   uint64_t seed, uint32_t maxL) {
    const uint64_t mask     = (1ULL << (2 * k)) - 1;
    const uint64_t hi_shift = 2ULL * (k - 1);

    // Reserve for ~1M sampled k-mers; a single 8 MB contiguous allocation.
    std::vector<uint64_t> sampled;
    sampled.reserve(1 << 20);
    uint64_t total = 0;

    bool used_libdeflate = false;

    // ── libdeflate path for .gz files ────────────────────────────────────
    // Only used for files < 50 MB compressed: larger files would need a
    // 200+ MB decompression buffer per thread, which wastes memory.
    // For large files, kseq (streaming zlib) is used instead.
    static constexpr size_t DEFLATE_MAX_CMP = 50ULL * 1024 * 1024;
    const bool is_gz = path.size() > 3 &&
                       path.compare(path.size() - 3, 3, ".gz") == 0;
    if (is_gz) {
        int fd = open(path.c_str(), O_RDONLY);
        if (fd >= 0) {
            struct stat st;
            if (fstat(fd, &st) == 0 && st.st_size > 0 &&
                static_cast<size_t>(st.st_size) <= DEFLATE_MAX_CMP) {
                std::vector<uint8_t> cmp(static_cast<size_t>(st.st_size));
                size_t nread = 0;
                while (nread < cmp.size()) {
                    ssize_t r = read(fd, cmp.data() + nread, cmp.size() - nread);
                    if (r <= 0) break;
                    nread += r;
                }
                close(fd);
                if (nread == cmp.size()) {
                    // Decompress: try 4×, grow to 8× if needed.
                    size_t out_cap = cmp.size() * 4;
                    std::vector<char> ubuf(out_cap);
                    auto* dc = libdeflate_alloc_decompressor();
                    size_t actual = 0, in_used = 0;
                    for (int attempt = 0; attempt < 4; ++attempt) {
                        auto rc = libdeflate_gzip_decompress_ex(
                            dc, cmp.data(), cmp.size(),
                            ubuf.data(), out_cap, &in_used, &actual);
                        if (rc == LIBDEFLATE_SUCCESS) { used_libdeflate = true; break; }
                        if (rc == LIBDEFLATE_INSUFFICIENT_SPACE) {
                            out_cap *= 2; ubuf.resize(out_cap);
                        } else { break; }
                    }
                    libdeflate_free_decompressor(dc);
                    if (used_libdeflate) {
                        // Handle multi-member gzip (rare but possible)
                        const char* p = ubuf.data();
                        size_t rem    = actual;
                        parse_fastq_mem(p, rem, k, mask, hi_shift, sampled, total);
                        // If there are more gzip members in `cmp`, process them
                        while (in_used < cmp.size()) {
                            size_t out2 = out_cap;
                            std::vector<char> ubuf2(out2);
                            size_t act2 = 0, in2 = 0;
                            dc = libdeflate_alloc_decompressor();
                            for (int attempt = 0; attempt < 4; ++attempt) {
                                auto rc = libdeflate_gzip_decompress_ex(
                                    dc, cmp.data() + in_used, cmp.size() - in_used,
                                    ubuf2.data(), out2, &in2, &act2);
                                if (rc == LIBDEFLATE_SUCCESS) break;
                                if (rc == LIBDEFLATE_INSUFFICIENT_SPACE) {
                                    out2 *= 2; ubuf2.resize(out2);
                                } else { act2 = 0; break; }
                            }
                            libdeflate_free_decompressor(dc);
                            if (act2 == 0) break;
                            parse_fastq_mem(ubuf2.data(), act2, k, mask, hi_shift, sampled, total);
                            in_used += in2;
                        }
                    }
                } else { close(fd); }
            } else { close(fd); }
        }
    }

    // ── kseq / zlib fallback ─────────────────────────────────────────────
    if (!used_libdeflate) {
        gzFile fp = gzopen(path.c_str(), "r");
        if (fp) {
            kseq_t* ks = kseq_init(fp);
            while (kseq_read(ks) >= 0) {
                const char* seq = ks->seq.s;
                const int   len = (int)ks->seq.l;
                if (len >= k)
                    pmh_kmers(seq, len, k, mask, hi_shift, sampled, total);
            }
            kseq_destroy(ks);
            gzclose(fp);
        } else {
            std::cerr << "WARNING: cannot open " << path << "\n";
        }
    }

    // ── sort + run-length count → PMH sketch ─────────────────────────────
    // Sorting gives identical results to hash-map counting (addHash is
    // order-independent) while using contiguous memory and zero mmap resizes.
    std::sort(sampled.begin(), sampled.end());

    auto* sk = new Sketch::ProbMinHash4(m, k, seed, maxL);
    if (total > 0 && !sampled.empty()) {
        const double inv = 1.0 / static_cast<double>(total);
        size_t i = 0;
        while (i < sampled.size()) {
            size_t j = i + 1;
            while (j < sampled.size() && sampled[j] == sampled[i]) ++j;
            const uint32_t cnt = static_cast<uint32_t>(j - i);
            // No singleton filter: assembly sequences have all counts=1 and
            // filtering would empty the sketch. For reads, low-count error
            // k-mers carry negligible weight (cnt/total ≈ 1/N) anyway.
            sk->addHash(sampled[i], static_cast<double>(cnt) * inv);
            i = j;
        }
    }
    return sk;
}

static inline double mash_distance(double jaccard, int kmer) {
    if (jaccard >= 1.0) return 0.0;
    if (jaccard <= 0.0) return 1.0;
    return -std::log(2.0 * jaccard / (1.0 + jaccard))
           / static_cast<double>(kmer);
}

static void write_pair(const std::string& outPath,
                       const std::string& a, const std::string& b, double d)
{
    FILE* fp = std::fopen(outPath.c_str(), "w");
    if (!fp) err(errno, "cannot open output: %s", outPath.c_str());
    std::fprintf(fp, "%s\t%s\t%.6f\n", a.c_str(), b.c_str(), d);
    std::fclose(fp);
    std::cerr << "output: " << outPath << "\n";
}

// Stream a FASTA/FASTQ file and invoke `cb(seq, len)` for every record.
template<class CB>
static void stream_seq(const std::string& path, CB&& cb) {
    gzFile fp = gzopen(path.c_str(), "r");
    if (!fp) { std::cerr << "WARN: cannot open " << path << "\n"; return; }
    kseq_t* ks = kseq_init(fp);
    while (kseq_read(ks) >= 0)
        cb(ks->seq.s, static_cast<uint64_t>(ks->seq.l));
    kseq_destroy(ks);
    gzclose(fp);
}

// kssd_parameter_t cache.  generate_shuffle_dim() inside the constructor uses
// libc srand/rand which is *not* thread-safe; concurrent constructions corrupt
// the dim shuffle map and produce zero-Jaccard sketches.  We cache one
// parameter object per (half_k, half_subk, drlevel) tuple under a mutex so
// every thread receives the same fully-initialised instance.
static const Sketch::kssd_parameter_t&
get_kssd_params(int half_k, int half_subk, int drlevel) {
    static std::mutex mu;
    static std::vector<std::pair<std::tuple<int,int,int>,
                                 std::shared_ptr<Sketch::kssd_parameter_t>>> cache;
    std::lock_guard<std::mutex> lk(mu);
    const auto key = std::make_tuple(half_k, half_subk, drlevel);
    for (auto& [k, v] : cache) if (k == key) return *v;
    cache.emplace_back(key,
        std::make_shared<Sketch::kssd_parameter_t>(half_k, half_subk, drlevel));
    return *cache.back().second;
}

// ═══════════════════════════════════════════════════════════════════════════
//  §2  Per-algorithm sketch builders – all return heap-allocated T*
// ═══════════════════════════════════════════════════════════════════════════

static Sketch::MinHash* build_minhash(const std::string& path, const Args& a) {
    auto* sk = new Sketch::MinHash(a.kmerSize, a.minhashSize,
                                   static_cast<uint32_t>(a.seed), true);
    sk->fileName = path;
    stream_seq(path, [&](char* seq, uint64_t /*len*/) { sk->update(seq); });
    sk->finalize();
    return sk;
}

static Sketch::HyperLogLog* build_hll(const std::string& path, const Args& a) {
    auto* sk = new Sketch::HyperLogLog(a.hllBits, a.kmerSize);
    stream_seq(path, [&](char* seq, uint64_t /*len*/) { sk->update(seq); });
    return sk;
}

static Sketch::Kssd* build_kssd(const std::string& path, const Args& a) {
    auto* sk = new Sketch::Kssd(get_kssd_params(/*half_k*/ a.kmerSize / 2,
                                                /*half_subk*/ 6,
                                                /*drlevel*/ a.kssdDrlevel));
    sk->fileName = path;
    stream_seq(path, [&](char* seq, uint64_t /*len*/) { sk->update(seq); });
    return sk;
}

static Sketch::ProbMinHash4* build_probminhash(const std::string& path, const Args& a) {
    // Frequency-weighted PMH-norm: estimates J_P = Σ min(p_A,p_B)/Σ max(p_A,p_B).
    // J_P is robust to coverage variation (relative k-mer frequencies are
    // preserved under uniform subsampling), unlike Binary Jaccard.
    return pmh_build_weighted(path, a.pmhM, a.kmerSize, a.seed, a.pmhMaxL);
}

static Sketch::BinDash* build_bindash(const std::string& path, const Args& a) {
    auto* sk = new Sketch::BinDash(a.bdSketch64, a.kmerSize, a.bdBbits, a.seed);
    stream_seq(path, [&](char* seq, uint64_t len) { sk->update(seq, len); });
    sk->finalize();
    return sk;
}

static Sketch::SetSketch* build_setsketch(const std::string& path, const Args& a) {
    auto* sk = new Sketch::SetSketch(a.ssBits, a.ssBase, a.ssA, a.kmerSize,
                                     /*track_witnesses=*/false);
    stream_seq(path, [&](char* seq, uint64_t len) {
        sk->update(seq, static_cast<size_t>(len));
    });
    return sk;
}

static Sketch::FastKMV* build_fastkmv(const std::string& path, const Args& a) {
    auto* sk = new Sketch::FastKMV(a.fkmvK, a.kmerSize, a.seed);
    stream_seq(path, [&](char* seq, uint64_t len) { sk->update(seq, len); });
    return sk;
}

// ═══════════════════════════════════════════════════════════════════════════
//  §3  Pair-mode (-i fileA fileB)
// ═══════════════════════════════════════════════════════════════════════════
static int run_pair(const Args& a) {
    const std::string& fA = a.inputs[0];
    const std::string& fB = a.inputs[1];
    double d = 1.0;

    switch (a.algo) {
    case Algo::MINHASH: {
        std::unique_ptr<Sketch::MinHash> sa(build_minhash(fA, a));
        std::unique_ptr<Sketch::MinHash> sb(build_minhash(fB, a));
        d = sa->distance(sb.get());
        break;
    }
    case Algo::HLL: {
        std::unique_ptr<Sketch::HyperLogLog> sa(build_hll(fA, a));
        std::unique_ptr<Sketch::HyperLogLog> sb(build_hll(fB, a));
        d = mash_distance(sa->jaccard_index(*sb), a.kmerSize);
        break;
    }
    case Algo::KSSD: {
        std::unique_ptr<Sketch::Kssd> sa(build_kssd(fA, a));
        std::unique_ptr<Sketch::Kssd> sb(build_kssd(fB, a));
        d = sa->distance(sb.get());
        break;
    }
    case Algo::PROBMINHASH: {
        std::unique_ptr<Sketch::ProbMinHash4> sa(build_probminhash(fA, a));
        std::unique_ptr<Sketch::ProbMinHash4> sb(build_probminhash(fB, a));
        d = mash_distance(sa->jaccard_weighted(*sb), a.kmerSize);
        break;
    }
    case Algo::BINDASH: {
        std::unique_ptr<Sketch::BinDash> sa(build_bindash(fA, a));
        std::unique_ptr<Sketch::BinDash> sb(build_bindash(fB, a));
        d = mash_distance(sa->jaccard(*sb), a.kmerSize);
        break;
    }
    case Algo::SETSKETCH: {
        std::unique_ptr<Sketch::SetSketch> sa(build_setsketch(fA, a));
        std::unique_ptr<Sketch::SetSketch> sb(build_setsketch(fB, a));
        d = mash_distance(sa->jaccard_index(*sb), a.kmerSize);
        break;
    }
    case Algo::FASTKMV: {
        std::unique_ptr<Sketch::FastKMV> sa(build_fastkmv(fA, a));
        std::unique_ptr<Sketch::FastKMV> sb(build_fastkmv(fB, a));
        d = mash_distance(sa->jaccard(*sb), a.kmerSize);
        break;
    }
    default: return 1;
    }

    write_pair(a.output, fA, fB, d);
    return 0;
}

// ═══════════════════════════════════════════════════════════════════════════
//  §3b MinHash --index  (mirrors run_index_fastkmv)
//
//      Bottom-k hash values are directly the inverted-index keys — no stride
//      or witness subsampling needed since every hash is a unique k-mer id.
//
//      Candidate generation uses minCommon = ceil(minJac × K) calibrated to
//      the library's union-k Jaccard estimator.  Exact verification runs the
//      same union-bottom-K merge as MinHash::jaccard() to ensure distances
//      are byte-identical to the brute-force distance() output.
// ═══════════════════════════════════════════════════════════════════════════
static void run_index_minhash(const Args& a, const std::vector<std::string>& files) {
    const int N = static_cast<int>(files.size());
    std::vector<int>                   sketchSizes(N);
    std::vector<std::vector<uint64_t>> skKeys(N);

    const int actualThreads = std::min(a.threads, N);
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>> threadIdx(actualThreads);

    double t0 = get_sec();
    #pragma omp parallel num_threads(a.threads)
    {
        int tid = omp_get_thread_num();
        auto& localIdx = (tid < actualThreads) ? threadIdx[tid] : threadIdx[0];

        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N; ++t) {
            gzFile fp = gzopen(files[t].c_str(), "r");
            if (!fp) continue;
            kseq_t* ks = kseq_init(fp);

            Sketch::MinHash sk(a.kmerSize, a.minhashSize,
                               static_cast<uint32_t>(a.seed), /*rc=*/true);
            while (kseq_read(ks) >= 0) sk.update(ks->seq.s);
            kseq_destroy(ks);
            gzclose(fp);

            const auto& hashes = sk.getHashesSorted();   // calls finalize() internally
            const int sz = static_cast<int>(hashes.size());
            sketchSizes[t] = sz;
            skKeys[t].resize(sz);
            for (int i = 0; i < sz; ++i) {
                skKeys[t][i] = hashes[i];
                localIdx[hashes[i]].push_back(static_cast<uint32_t>(t));
            }
        }
    }
    double t1 = get_sec();
    std::cerr << "sketch + local index: " << t1 - t0 << " s\n";

    auto csrIdx = Sketch::buildCSRIndex<uint64_t>(threadIdx, a.threads);
    double t2 = get_sec();
    std::cerr << "build CSR index: " << t2 - t1 << " s\n";

    // minCommon threshold: calibrated so that c ≥ minCommon implies c' ≥ minJac·K
    // where c' is the union-K intersection (what jaccard() computes).
    // Since c ≥ c' always (full intersection ≥ union-K intersection), using the
    // same threshold ensures no false negatives: if c' ≥ minJac·K then c ≥ c'.
    const double p_exp  = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const int    minCommon = std::max(1,
        static_cast<int>(std::ceil(minJac * static_cast<double>(a.minhashSize))));
    std::cerr << "pruning: minCommon=" << minCommon << "/" << a.minhashSize
              << "  (minJac=" << minJac << ", mashD<" << a.maxDist
              << ", k=" << a.kmerSize << ")\n";

    // Exact verification: replicate MinHash::jaccard() union-K merge intersection.
    // This guarantees distances byte-identical to the brute-force distance() path.
    //
    // Early-exit (every 32 steps): if c + min(remaining_A, remaining_B) < minCommon,
    // the union-K intersection can never reach minCommon → abort merge.  Checking
    // every 32 steps keeps per-step overhead ~1/32 while still catching borderline
    // pairs that passed the inverted-index filter but fail exact verification.
    const int K = a.minhashSize;
    const int mc = minCommon;
    auto exactJaccardFn = [&, K_cap = K, mc_cap = mc](int i, int j) -> double {
        const auto& hi = skKeys[i];
        const auto& hj = skKeys[j];
        const int si = static_cast<int>(hi.size());
        const int sj = static_cast<int>(hj.size());
        int ii = 0, jj = 0, c = 0, denom = 0;
        while (denom < K_cap && ii < si && jj < sj) {
            if      (hi[ii] < hj[jj]) { ii++; }
            else if (hi[ii] > hj[jj]) { jj++; }
            else                      { c++; ii++; jj++; }
            denom++;
            // Check every 32 steps: abort if impossible to reach minCommon
            if (__builtin_expect((denom & 31) == 0, 0) &&
                c + std::min(si - ii, sj - jj) < mc_cap) return 0.0;
        }
        if (denom < K_cap) {
            denom += (si - ii) + (sj - jj);
            if (denom > K_cap) denom = K_cap;
        }
        return (denom <= 0) ? 0.0 : static_cast<double>(c) / denom;
    };
    auto minCommonFn = [minCommon](int) { return minCommon; };

    double t3 = get_sec();
    Sketch::computeDistancesExact<uint64_t>(csrIdx, skKeys, files,
        N, a.kmerSize, a.maxDist, exactJaccardFn, minCommonFn, a.output, a.threads);
    double t4 = get_sec();
    std::cerr << "dist time: "  << t4 - t3 << " s\n";
    std::cerr << "total time: " << t4 - t0 << " s\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  §4  list_allpairs<T> – generic O(N²) baseline used by algorithms that have
//      no specialised list-mode path (MinHash, HLL, ProbMinHash, FastKMV in
//      no-index mode).  BinDash uses its own fast path (see §9).
// ═══════════════════════════════════════════════════════════════════════════
template<class T, class BuildFn, class DistFn>
static void list_allpairs(const Args& a,
                          const std::vector<std::string>& files,
                          BuildFn build, DistFn dist)
{
    const int N = static_cast<int>(files.size());
    std::vector<std::unique_ptr<T>> sks(N);

    double t0 = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(dynamic)
    for (int i = 0; i < N; ++i) sks[i].reset(build(files[i], a));
    std::cerr << "sketch time: " << get_sec() - t0 << " s\n";

    double t1 = get_sec();
    std::vector<std::string> bufs(a.threads);

    #pragma omp parallel for num_threads(a.threads) schedule(dynamic, 1)
    for (int i = 0; i < N; ++i) {
        const int tid = omp_get_thread_num();
        std::string& buf = bufs[tid];
        for (int j = i + 1; j < N; ++j) {
            const double d = dist(sks[i].get(), sks[j].get());
            if (d < a.maxDist) {
                char line[2048];
                int len = std::snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                                        files[i].c_str(), files[j].c_str(), d);
                buf.append(line, static_cast<size_t>(len));
            }
        }
    }

    FILE* fp = std::fopen(a.output.c_str(), "w");
    if (!fp) err(errno, "cannot open output: %s", a.output.c_str());
    std::setvbuf(fp, nullptr, _IOFBF, 1 << 22);
    for (auto& b : bufs) std::fwrite(b.data(), 1, b.size(), fp);
    std::fclose(fp);
    std::cerr << "dist time:   " << get_sec() - t1 << " s\n";
    std::cerr << "output: " << a.output << "\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  §5  FastKMV --index  (mirrors test_FastKMV.cpp exactly)
//
//      1. Each genome's bottom-k hash list is published into a per-thread
//         local inverted index (hash → list of genome IDs).
//      2. Local indices are merged into a single CSR-format inverted index.
//      3. computeDistancesExact() walks each genome's keys, follows the posting
//         list for each key, and uses stamp/epoch counting to accumulate the
//         intersection size c against every other genome in
//         O(N · K · avg_pl).
//
//  Exact verification replicates FastKMV::jaccard() union-K merge estimator:
//    Jaccard = c / min(|A_k ∪ B_k|, K)
//  This ensures distances from --index and pairwise paths are numerically
//  identical for the same sketch pair.
// ═══════════════════════════════════════════════════════════════════════════
static void run_index_fastkmv(const Args& a, const std::vector<std::string>& files) {
    const int N = static_cast<int>(files.size());
    std::vector<int>                   sketchSizes(N);
    std::vector<std::vector<uint64_t>> skKeys(N);

    const int actualThreads = std::min(a.threads, N);
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>> threadIdx(actualThreads);

    double t0 = get_sec();
    #pragma omp parallel num_threads(a.threads)
    {
        int tid = omp_get_thread_num();
        auto& localIdx = (tid < actualThreads) ? threadIdx[tid] : threadIdx[0];

        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N; ++t) {
            gzFile fp = gzopen(files[t].c_str(), "r");
            if (!fp) continue;
            kseq_t* ks = kseq_init(fp);

            Sketch::FastKMV sk(a.fkmvK, a.kmerSize, a.seed);
            while (kseq_read(ks) >= 0)
                sk.update(ks->seq.s, ks->seq.l);

            kseq_destroy(ks);
            gzclose(fp);

            const uint64_t* regs = sk.getRegisters();
            uint32_t sz = sk.size();
            sketchSizes[t] = static_cast<int>(sz);
            skKeys[t].resize(sz);
            for (uint32_t i = 0; i < sz; ++i) {
                skKeys[t][i] = regs[i];
                localIdx[regs[i]].push_back(static_cast<uint32_t>(t));
            }
        }
    }
    double t1 = get_sec();
    std::cerr << "sketch + local index: " << t1 - t0 << " s\n";

    auto csrIdx = Sketch::buildCSRIndex<uint64_t>(threadIdx, a.threads);
    double t2 = get_sec();
    std::cerr << "build CSR index: " << t2 - t1 << " s\n";

    // Mash-distance ↔ Jaccard pruning identical to test_FastKMV.cpp:
    //   minJac = exp(-k·D) / (2 - exp(-k·D)),  minCommon = ceil(minJac · K)
    const double p_exp  = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const int    minCommon = std::max(1,
        static_cast<int>(std::ceil(minJac * static_cast<double>(a.fkmvK))));
    std::cerr << "pruning: minCommon=" << minCommon << "/" << a.fkmvK
              << "  (minJac=" << minJac << ", mashD<" << a.maxDist
              << ", k=" << a.kmerSize << ")\n";

    // Exact verification: replicate FastKMV::jaccard() union-K merge estimator.
    // Jaccard = c / min(|A_k ∪ B_k|, K), identical to the pairwise path.
    const int K_fkmv = a.fkmvK;
    const int mc_fkmv = minCommon;
    auto exactJaccardFn = [&, K_cap = K_fkmv, mc_cap = mc_fkmv](int i, int j) -> double {
        const auto& hi = skKeys[i];
        const auto& hj = skKeys[j];
        const int si = static_cast<int>(hi.size());
        const int sj = static_cast<int>(hj.size());
        int ii = 0, jj = 0, c = 0, denom = 0;
        while (denom < K_cap && ii < si && jj < sj) {
            if      (hi[ii] < hj[jj]) { ii++; }
            else if (hi[ii] > hj[jj]) { jj++; }
            else                      { c++; ii++; jj++; }
            denom++;
            if (__builtin_expect((denom & 31) == 0, 0) &&
                c + std::min(si - ii, sj - jj) < mc_cap) return 0.0;
        }
        if (denom < K_cap) {
            denom += (si - ii) + (sj - jj);
            if (denom > K_cap) denom = K_cap;
        }
        return (denom <= 0) ? 0.0 : static_cast<double>(c) / denom;
    };
    auto minCommonFn = [minCommon](int) { return minCommon; };

    double t3 = get_sec();
    Sketch::computeDistancesExact<uint64_t>(csrIdx, skKeys, files,
        N, a.kmerSize, a.maxDist, exactJaccardFn, minCommonFn, a.output, a.threads);
    double t4 = get_sec();
    std::cerr << "dist time: "  << t4 - t3 << " s\n";
    std::cerr << "total time: " << t4 - t0 << " s\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  §6  ProbMinHash --index  (frequency-weighted J_P mode)
//
//  Two-pass frequency-weighted pipeline:
//    Pass-1: count k-mer occurrences (1/PMH_SDIV hash-sampling for memory).
//    Pass-2: insert each unique k-mer once with weight = count/total_kmers.
//  This estimates J_P = Σ min(p_A,p_B) / Σ max(p_A,p_B), which is
//  invariant to sequencing depth — the core PMH advantage over BJ estimators.
//  Index key = winners_[k] XOR (k × golden-ratio) via getWinnerIndexKeys().
// ═══════════════════════════════════════════════════════════════════════════

static void run_index_probmh(const Args& a, const std::vector<std::string>& files) {
    const int N = static_cast<int>(files.size());
    const int mSize = static_cast<int>(a.pmhM);
    std::vector<int> sketchSizes(N, mSize);
    std::vector<std::vector<uint64_t>> skKeys(N);

    const int actualThreads = std::min(a.threads, N);
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>> threadIdx(actualThreads);

    double t0 = get_sec();

    #pragma omp parallel num_threads(a.threads)
    {
        int tid = omp_get_thread_num();
        auto& localIdx = (tid < actualThreads) ? threadIdx[tid] : threadIdx[0];

        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N; ++t) {
            std::unique_ptr<Sketch::ProbMinHash4> sk(
                pmh_build_weighted(files[t], a.pmhM, a.kmerSize, a.seed, a.pmhMaxL));

            sk->getWinnerIndexKeys(skKeys[t]);
            for (uint64_t key : skKeys[t])
                localIdx[key].push_back(static_cast<uint32_t>(t));
        }
    }

    double t1 = get_sec();
    std::cerr << "sketch + local index: " << t1 - t0 << " s\n";

    auto csrIdx = Sketch::buildCSRIndex<uint64_t>(threadIdx, a.threads);
    double t2 = get_sec();
    std::cerr << "build CSR index: " << t2 - t1 << " s\n";

    const double p_exp  = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const int minCommon = std::max(1, static_cast<int>(std::ceil(minJac * a.pmhM)));
    std::cerr << "pruning: minCommon=" << minCommon << "/" << a.pmhM
              << "  (minJac=" << minJac << ", mashD<" << a.maxDist
              << ", k=" << a.kmerSize << ")\n";

    auto jaccardFn = [mSize](int c, int, int) -> double {
        return Sketch::ProbMinHash4::jaccardFromCommon(c, static_cast<uint32_t>(mSize));
    };
    auto minCommonFn = [minCommon](int) { return minCommon; };

    double t3 = get_sec();
    Sketch::computeDistances<uint64_t>(csrIdx, skKeys, sketchSizes, files,
        N, a.kmerSize, a.maxDist, jaccardFn, minCommonFn, a.output, a.threads);
    double t4 = get_sec();
    std::cerr << "dist time: "  << t4 - t3 << " s\n";
    std::cerr << "total time: " << t4 - t0 << " s\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  §7  SetSketch --index  (mirrors test_SetSketch.cpp – witnesses extracted
//      during sketching, sort by cardinality + reorder, build CSR inverted
//      index, singleton-filter skKeys, then either:
//        (A) inverted-index + exact verify   (expectedOverlap >= 20)
//        (B) all-pairs SIMD batch + size-ratio prune + suffix-sum early abort)
// ═══════════════════════════════════════════════════════════════════════════
static void run_index_setsketch(const Args& a, const std::vector<std::string>& files_in) {
    using namespace std;

    const int N = static_cast<int>(files_in.size());
    if (N == 0) { cerr << "ERROR: empty file list\n"; return; }

    Sketch::SetSketch proto(a.ssBits, a.ssBase, a.ssA, a.kmerSize,
                            /*track_witnesses=*/false);
    const int    m      = proto.getM();
    const double factor = proto.getFactor();
    double bip_buf[64];
    memcpy(bip_buf, proto.getBaseInvPow(), 64 * sizeof(double));
    const double* bip = bip_buf;

    // One inverted-index key per occupied register: key = p * 64 + core_[p]
    // P(key_A[p] == key_B[p]) ≈ J_SetSketch  ← same property as MinHash keys.
    // This eliminates the witness hash arrays entirely, saving 64 KB per sketch.
    cerr << "registers=" << m << "  keyEncoding=(reg_idx*64+core_val)\n";

    // ── Phase 1: Sketch construction + witness extraction ───────────────────
    vector<double>  sizes(N, 0.0);
    vector<uint8_t> flat_cores(static_cast<size_t>(N) * m, 0);
    vector<vector<uint64_t>> skKeys(N);
    vector<string>  fileList = files_in;     // mutable copy (sorted in place below)

    double t0 = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(dynamic)
    for (int t = 0; t < N; ++t) {
        Sketch::SetSketch sk(a.ssBits, a.ssBase, a.ssA, a.kmerSize,
                             /*track_witnesses=*/false);
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) continue;
        kseq_t* ks = kseq_init(fp);
        while (kseq_read(ks) >= 0)
            sk.update(ks->seq.s, static_cast<size_t>(ks->seq.l));
        kseq_destroy(ks);
        gzclose(fp);

        sizes[t] = sk.cardinality();
        const uint8_t* core = sk.getCore().data();
        memcpy(&flat_cores[static_cast<size_t>(t) * m], core, (size_t)m);
        skKeys[t].reserve(m);
        for (int p = 0; p < m; ++p)
            if (core[p] != 0)
                skKeys[t].push_back(static_cast<uint64_t>(p) * 64ULL + core[p]);
    }
    double t1 = get_sec();
    cerr << "sketch time: " << t1 - t0 << " s\n";

    // ── Sort by cardinality descending + reorder all per-genome arrays ──────
    // Sorting enables O(N log N) tile-level / inner-loop pruning later: when
    // sizes[j] / sizes[i] < minJac, the containment bound forbids any pair (i, j')
    // with j' >= j from clearing the Jaccard threshold, so we can `break` instead
    // of `continue`.
    {
        vector<int> order(N);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(),
             [&](int aa, int bb) { return sizes[aa] > sizes[bb]; });

        vector<string>  sF(N);
        vector<double>  sS(N, 0.0);
        vector<uint8_t> sC(static_cast<size_t>(N) * m, 0);
        vector<vector<uint64_t>> sK(N);
        for (int ni = 0; ni < N; ++ni) {
            const int oi = order[ni];
            sF[ni] = fileList[oi];
            sS[ni] = sizes[oi];
            memcpy(&sC[static_cast<size_t>(ni) * m],
                   &flat_cores[static_cast<size_t>(oi) * m], (size_t)m);
            sK[ni] = std::move(skKeys[oi]);
        }
        fileList.swap(sF);
        sizes.swap(sS);
        flat_cores.swap(sC);
        skKeys.swap(sK);
    }

    // ── Precompute per-register suffix sums (for early abort in jaccardBatch) ──
    static const int TAIL_STEP = 64;
    const int nCP = m / TAIL_STEP + 1;
    vector<double> tailSums(static_cast<size_t>(N) * nCP, 0.0);

    double tpre = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(static)
    for (int t = 0; t < N; ++t) {
        const uint8_t* core = &flat_cores[static_cast<size_t>(t) * m];
        double* ts = &tailSums[static_cast<size_t>(t) * nCP];
        ts[nCP - 1] = 0.0;
        for (int cp = nCP - 2; cp >= 0; --cp) {
            const int start = cp * TAIL_STEP;
            const int end   = min(start + TAIL_STEP, m);
            double blockSum = 0.0;
            for (int kk = start; kk < end; ++kk) blockSum += bip[core[kk]];
            ts[cp] = ts[cp + 1] + blockSum;
        }
    }
    cerr << "suffix sums: " << get_sec() - tpre << " s\n";

    // ── Decide mode: inverted index needs ≥20 expected witness overlap ──────
    // expectedOverlap >= 20 ⇒ P(missed valid pair) < exp(-20) ≈ 2e-9.
    const double p_exp  = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    // Expected number of matching register-value keys between two sequences at
    // the Jaccard threshold: E = m * minJac  (P(match per register) ≈ J).
    const double expectedOverlap = static_cast<double>(m) * minJac;
    const bool   useInvIdx = (expectedOverlap >= 20.0);
    cerr << "minJac=" << minJac << "  expectedRegValOverlap=" << expectedOverlap << "\n";

    double t5 = get_sec();

    if (useInvIdx) {
        // ───────────── MODE A: inverted index (witness hashes) + exact ─────
        cerr << "mode: INVERTED INDEX (register-value keys)\n";

        const int actualThreads = min(a.threads, N);
        vector<phmap::flat_hash_map<uint64_t, vector<uint32_t>>> threadIdx(actualThreads);

        double ti0 = get_sec();
        #pragma omp parallel num_threads(a.threads)
        {
            int tid = omp_get_thread_num();
            if (tid < actualThreads) {
                auto& localIdx = threadIdx[tid];
                #pragma omp for schedule(static)
                for (int t = 0; t < N; ++t)
                    for (uint64_t key : skKeys[t])
                        localIdx[key].push_back(static_cast<uint32_t>(t));
            }
        }
        double ti1 = get_sec();
        cerr << "local index: " << ti1 - ti0 << " s\n";

        auto csrIdx = Sketch::buildCSRIndex<uint64_t>(threadIdx, a.threads);
        double ti2 = get_sec();
        cerr << "CSR build: " << ti2 - ti1 << " s\n";

        // Singleton-filter skKeys: a key present in exactly one sketch can never
        // produce a candidate pair. Dropping them cuts skKeys footprint by ~75%
        // and speeds up candidate generation by the same factor; P(miss) is
        // unchanged because every shared pair shares only non-singleton keys.
        {
            size_t kbefore = 0, kafter = 0;
            #pragma omp parallel for num_threads(a.threads) schedule(dynamic, 64) \
                                     reduction(+:kbefore,kafter)
            for (int t = 0; t < N; ++t) {
                kbefore += skKeys[t].size();
                vector<uint64_t> kept;
                kept.reserve(skKeys[t].size() / 5);
                for (uint64_t key : skKeys[t])
                    if (csrIdx.postIdx.count(key))
                        kept.push_back(key);
                kafter += kept.size();
                skKeys[t] = std::move(kept);
            }
            cerr << "skKeys singleton-filtered: " << kbefore << " -> " << kafter
                 << " (" << (kbefore > 0 ? 100.0 * kafter / kbefore : 0.0)
                 << "% retained, freed ~"
                 << (kbefore - kafter) * 8 / (1 << 20) << " MB)\n";
        }

        const double sd = std::sqrt(std::max(0.0, expectedOverlap * (1.0 - minJac)));
        const int minCommon = max(1, static_cast<int>(std::floor(expectedOverlap - 8.0 * sd)));
        cerr << "minCommon=" << minCommon << " (expected=" << expectedOverlap
             << ", 8sigma=" << 8.0 * sd << ")\n";

        auto exactJaccardFn = [&](int i, int j) -> double {
            const double si = sizes[i], sj = sizes[j];
            if (si > 0.0 && sj / si < minJac) return -1.0;
            return Sketch::SetSketch::jaccardFromCoresBatch(
                &flat_cores[static_cast<size_t>(i) * m],
                &flat_cores[static_cast<size_t>(j) * m],
                m, bip, factor, si, sj, minJac,
                &tailSums[static_cast<size_t>(i) * nCP],
                &tailSums[static_cast<size_t>(j) * nCP],
                TAIL_STEP);
        };
        auto minCommonFn = [minCommon](int) { return minCommon; };

        Sketch::computeDistancesExact<uint64_t>(
            csrIdx, skKeys, fileList,
            N, a.kmerSize, a.maxDist, exactJaccardFn, minCommonFn,
            a.output, a.threads);
    } else {
        // ───────────── MODE B: all-pairs SIMD batch + suffix-sum early abort ─
        cerr << "mode: ALL-PAIRS (SIMD batch + suffix-sum early abort)\n";
        { vector<vector<uint64_t>>().swap(skKeys); }   // free witnesses

        // Resolve output path (handle directory)
        string finalPath = a.output;
        {
            struct stat st;
            if (stat(finalPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
                if (finalPath.back() != '/') finalPath += '/';
                finalPath += "rabbitsketch.dist";
            }
        }
        cerr << "output: " << finalPath << "\n";

        const int TILE = 128;
        const int nTiles = (N + TILE - 1) / TILE;
        int progress = N / 20;
        if (progress < 1) progress = 1;

        const int actualThreads = min(a.threads, N);
        const int procId = static_cast<int>(getpid());
        vector<string> partPaths(actualThreads);
        for (int tid = 0; tid < actualThreads; ++tid)
            partPaths[tid] = finalPath + ".part." + to_string(procId)
                                       + "." + to_string(tid);

        #pragma omp parallel num_threads(actualThreads)
        {
            const int tid = omp_get_thread_num();
            FILE* tf = fopen(partPaths[tid].c_str(), "wb");
            if (!tf) {
                #pragma omp critical
                cerr << "ERROR: cannot open temp file: " << partPaths[tid] << "\n";
            } else {
                setvbuf(tf, nullptr, _IOFBF, 1 << 22);
                string buf;
                buf.reserve(1 << 22);

                #pragma omp for schedule(dynamic, 1)
                for (int tileI = 0; tileI < nTiles; ++tileI) {
                    const int iBeg = tileI * TILE;
                    const int iEnd = min(N, iBeg + TILE);
                    for (int tileJ = tileI; tileJ < nTiles; ++tileJ) {
                        const int jBeg = tileJ * TILE;
                        const int jEnd = min(N, jBeg + TILE);

                        if (tileJ > tileI && sizes[iEnd - 1] > 0.0 &&
                            sizes[jBeg] / sizes[iEnd - 1] < minJac)
                            break;

                        for (int i = iBeg; i < iEnd; ++i) {
                            int jStart = (tileI == tileJ) ? max(i + 1, jBeg) : jBeg;
                            if (jStart >= jEnd) continue;

                            const uint8_t* core_i = &flat_cores[static_cast<size_t>(i) * m];
                            const double si  = sizes[i];
                            const double* ts_i = &tailSums[static_cast<size_t>(i) * nCP];

                            for (int j = jStart; j < jEnd; ++j) {
                                const double sj = sizes[j];
                                const double maxJacBySize = (si > 0.0) ? (sj / si) : 0.0;
                                if (maxJacBySize < minJac) break;

                                if (j + 4 < jEnd) {
                                    __builtin_prefetch(&flat_cores[static_cast<size_t>(j + 4) * m], 0, 1);
                                    __builtin_prefetch(&tailSums[static_cast<size_t>(j + 4) * nCP], 0, 1);
                                }

                                const uint8_t* c2 = &flat_cores[static_cast<size_t>(j) * m];
                                const double* ts_j = &tailSums[static_cast<size_t>(j) * nCP];
                                double jac = Sketch::SetSketch::jaccardFromCoresBatch(
                                    core_i, c2, m, bip, factor, si, sj, minJac,
                                    ts_i, ts_j, TAIL_STEP);
                                if (jac < minJac) continue;

                                double dist = (jac >= 1.0) ? 0.0
                                    : -std::log(2.0 * jac / (1.0 + jac))
                                      / static_cast<double>(a.kmerSize);
                                if (dist < a.maxDist) {
                                    char line[1024];
                                    int len = snprintf(line, sizeof(line),
                                        "%s\t%s\t%.6f\n",
                                        fileList[i].c_str(), fileList[j].c_str(), dist);
                                    buf.append(line, static_cast<size_t>(len));
                                }
                                if (buf.size() > (1 << 22)) {
                                    fwrite(buf.data(), 1, buf.size(), tf);
                                    buf.clear();
                                }
                            }
                        }
                    }
                    const int iProgress = min(iEnd, N);
                    if (iProgress % progress == 0) {
                        #pragma omp critical
                        cerr << "  dist " << iProgress << " / " << N << "\n";
                    }
                }
                if (!buf.empty()) fwrite(buf.data(), 1, buf.size(), tf);
                fclose(tf);
            }
        }

        // Merge per-thread temp files into final output
        FILE* fout = fopen(finalPath.c_str(), "wb");
        if (!fout) {
            cerr << "ERROR: cannot open " << finalPath << "\n";
            for (int tid = 0; tid < actualThreads; ++tid) remove(partPaths[tid].c_str());
            return;
        }
        setvbuf(fout, nullptr, _IOFBF, 1 << 24);
        vector<char> copyBuf(1 << 20);
        for (int tid = 0; tid < actualThreads; ++tid) {
            FILE* part = fopen(partPaths[tid].c_str(), "rb");
            if (!part) continue;
            while (true) {
                size_t got = fread(copyBuf.data(), 1, copyBuf.size(), part);
                if (got == 0) break;
                fwrite(copyBuf.data(), 1, got, fout);
            }
            fclose(part);
            remove(partPaths[tid].c_str());
        }
        fclose(fout);
    }

    double t6 = get_sec();
    cerr << "dist time: "  << t6 - t5 << " s\n";
    cerr << "total time: " << t6 - t0 << " s\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  §8  Kssd --index  (mirrors test_Kssd.cpp – inlined SIMD encoder + 64-shard
//      merge + flat CSR posting array + stamp-based dist kernel + 5-col output)
// ═══════════════════════════════════════════════════════════════════════════
static void run_index_kssd(const Args& a, const std::vector<std::string>& files)
{
    using namespace std;

    // Build (or fetch from cache) shared kssd parameters.  Because Kssd's
    // raw test does `kssd_parameter_t P;` (uses defaults) we reproduce that
    // by feeding the user-controllable drlevel through but keeping
    // half_k=10, half_subk=6 hard-coded — same defaults as the Sketch
    // library and as test_Kssd.
    const Sketch::kssd_parameter_t& P = get_kssd_params(/*half_k*/ 10,
                                                        /*half_subk*/ 6,
                                                        /*drlevel*/ a.kssdDrlevel);
    const int    kmerSize  = P.kmer_size;
    const bool   use64     = (P.half_k - P.drlevel) > 8;
    const uint64_t tupmask = P.tupmask, domask = P.domask;
    const uint64_t umask0  = P.undomask0, umask1 = P.undomask1;
    const int    revmov    = P.rev_add_move;
    const int    domS      = P.half_outctx_len * 2;
    const int    drS       = P.drlevel * 4;
    const int    undoS     = kmerSize * 2 - P.half_outctx_len * 4;
    const int    dimStart  = P.dim_start;
    const auto&  smap      = P.shuffled_map;

    std::cerr << "half_k=" << P.half_k << " half_subk=" << P.half_subk
              << " drlevel=" << P.drlevel
              << " kmer=" << kmerSize << " use64=" << use64 << "\n";

    // Scalar encode table for non-SIMD tail
    uint8_t BM_U8[256];
    std::memset(BM_U8, 0xFF, sizeof(BM_U8));
    BM_U8['A'] = BM_U8['a'] = 0;
    BM_U8['C'] = BM_U8['c'] = 1;
    BM_U8['G'] = BM_U8['g'] = 2;
    BM_U8['T'] = BM_U8['t'] = 3;

    const int N = static_cast<int>(files.size());
    int progress = N / 20; if (progress < 1) progress = 1;

    struct GSketch {
        std::vector<uint32_t> h32;
        std::vector<uint64_t> h64;
    };
    std::vector<GSketch> sketches(N);

    const int actualThreads = std::min(a.threads, N);
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>>
        threadIdx(actualThreads);

    double t0 = get_sec();

    #pragma omp parallel num_threads(a.threads)
    {
        const int tid = omp_get_thread_num();
        phmap::flat_hash_set<uint32_t> hs32;
        phmap::flat_hash_set<uint64_t> hs64;
        hs32.reserve(1 << 16);
        hs64.reserve(1 << 16);
        auto& localIdx = (tid < actualThreads) ? threadIdx[tid] : threadIdx[0];

#if defined(__AVX512BW__)
        const __m512i vlut512 = _mm512_broadcast_i32x4(_mm_setr_epi8(
            -1, 0, -1, 1, 3, -1, -1, 2,
            -1,-1, -1,-1,-1, -1, -1,-1));
        const __m512i vmask512 = _mm512_set1_epi8(0x0F);
#elif defined(__AVX2__)
        const __m256i vlut256 = _mm256_broadcastsi128_si256(_mm_setr_epi8(
            -1, 0, -1, 1, 3, -1, -1, 2,
            -1,-1, -1,-1,-1, -1, -1,-1));
        const __m256i vmask256 = _mm256_set1_epi8(0x0F);
#endif

        std::vector<char> fileBuf;
        std::vector<char> seqBuf;

        // k-mer chunk processor inlined to keep the hot loop branchless.
        auto processSeq = [&](const char* s, int len) {
            if (len < kmerSize) return;
            alignas(64) uint8_t enc[4096];
            uint64_t fwd = 0, rev = 0;
            int run = 0;
            for (int cs = 0; cs < len; cs += 4096) {
                const int clen = std::min(4096, len - cs);
                const char* cp = s + cs;
                int p = 0;
#if defined(__AVX512BW__)
                for (; p + 64 <= clen; p += 64) {
                    __m512i v = _mm512_loadu_si512(cp + p);
                    v = _mm512_shuffle_epi8(vlut512,
                            _mm512_and_si512(v, vmask512));
                    _mm512_storeu_si512(enc + p, v);
                }
#elif defined(__AVX2__)
                for (; p + 32 <= clen; p += 32) {
                    __m256i v = _mm256_loadu_si256((const __m256i*)(cp + p));
                    v = _mm256_shuffle_epi8(vlut256,
                            _mm256_and_si256(v, vmask256));
                    _mm256_storeu_si256((__m256i*)(enc + p), v);
                }
#endif
                for (; p < clen; ++p)
                    enc[p] = BM_U8[(unsigned char)cp[p]];
                for (int i = 0; i < clen; i++) {
                    const uint8_t b = enc[i];
                    if (__builtin_expect(b <= 3, 1)) {
                        fwd = ((fwd << 2) | b) & tupmask;
                        rev = (rev >> 2) + (((uint64_t)(b ^ 3)) << revmov);
                        if (__builtin_expect(++run >= kmerSize, 1)) {
                            const uint64_t u = (fwd < rev) ? fwd : rev;
                            const uint32_t d = static_cast<uint32_t>(
                                (u & domask) >> domS);
                            auto sit = smap.find(d);
                            if (__builtin_expect(sit == smap.end(), 1)) continue;
                            const uint64_t h =
                                (((u & umask0)
                                  | ((u & umask1) << undoS))
                                 >> drS)
                                | static_cast<uint64_t>(sit->second - dimStart);
                            if (use64) hs64.emplace(h);
                            else       hs32.emplace(static_cast<uint32_t>(h));
                        }
                    } else {
                        run = 0; fwd = 0; rev = 0;
                    }
                }
            }
        };

        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N; t++) {
            hs32.clear(); hs64.clear();

            struct stat st;
            if (stat(files[t].c_str(), &st) != 0) continue;
            const size_t fsize = static_cast<size_t>(st.st_size);
            if (fsize < 2) continue;
            fileBuf.resize(fsize);
            FILE* ff = std::fopen(files[t].c_str(), "rb");
            if (!ff) continue;
            const size_t rd = std::fread(fileBuf.data(), 1, fsize, ff);
            std::fclose(ff);

            const bool isGz =
                (static_cast<unsigned char>(fileBuf[0]) == 0x1f
              && static_cast<unsigned char>(fileBuf[1]) == 0x8b);

            if (isGz) {
                gzFile fp = gzopen(files[t].c_str(), "r");
                if (!fp) continue;
                gzbuffer(fp, 1 << 18);
                kseq_t* ks = kseq_init(fp);
                while (kseq_read(ks) >= 0)
                    processSeq(ks->seq.s, static_cast<int>(ks->seq.l));
                kseq_destroy(ks);
                gzclose(fp);
            } else {
                const char* data = fileBuf.data();
                size_t pos = 0;
                while (pos < rd) {
                    if (data[pos] != '>') { pos++; continue; }
                    while (pos < rd && data[pos] != '\n') pos++;
                    if (pos < rd) pos++;
                    seqBuf.clear();
                    while (pos < rd && data[pos] != '>') {
                        size_t ls = pos;
                        while (pos < rd && data[pos] != '\n'
                               && data[pos] != '\r') pos++;
                        seqBuf.insert(seqBuf.end(), data + ls, data + pos);
                        while (pos < rd && (data[pos] == '\n'
                               || data[pos] == '\r')) pos++;
                    }
                    processSeq(seqBuf.data(), static_cast<int>(seqBuf.size()));
                }
            }

            // hash set → sketch + inverted index in one pass (no sort needed)
            if (use64) {
                sketches[t].h64.reserve(hs64.size());
                for (uint64_t h : hs64) {
                    sketches[t].h64.push_back(h);
                    localIdx[h].push_back(static_cast<uint32_t>(t));
                }
            } else {
                sketches[t].h32.reserve(hs32.size());
                for (uint32_t h : hs32) {
                    sketches[t].h32.push_back(h);
                    localIdx[static_cast<uint64_t>(h)]
                        .push_back(static_cast<uint32_t>(t));
                }
            }
            if (t % progress == 0)
                std::cerr << "  sketch " << t << " / " << N << "\n";
        }
    }
    double t1 = get_sec();
    std::cerr << "sketch + local index: " << t1 - t0 << " s\n";

    // ── Parallel merge into sharded global inverted index ────────────────────
    double t2 = get_sec();
    const int NUM_SHARDS = 64;
    const uint64_t SHARD_MASK = NUM_SHARDS - 1;
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>>
        invShards(NUM_SHARDS);
    {
        omp_lock_t locks[64];
        for (int s = 0; s < NUM_SHARDS; s++) omp_init_lock(&locks[s]);

        #pragma omp parallel num_threads(a.threads)
        {
            int tid = omp_get_thread_num();
            if (tid < actualThreads) {
                for (auto& [key, vec] : threadIdx[tid]) {
                    int shard = static_cast<int>(key & SHARD_MASK);
                    omp_set_lock(&locks[shard]);
                    auto [it, ins] = invShards[shard].try_emplace(key, std::move(vec));
                    if (!ins)
                        it->second.insert(it->second.end(), vec.begin(), vec.end());
                    omp_unset_lock(&locks[shard]);
                }
                phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>()
                    .swap(threadIdx[tid]);
            }
        }
        for (int s = 0; s < NUM_SHARDS; s++) omp_destroy_lock(&locks[s]);
        threadIdx.clear();
    }
    size_t totalUnique = 0;
    for (auto& sh : invShards) totalUnique += sh.size();
    std::cerr << "merge index: " << totalUnique << " unique hashes, "
              << get_sec() - t2 << " s\n";

    // Singleton removal
    size_t totalBefore = totalUnique, totalAfter = 0;
    #pragma omp parallel for num_threads(a.threads) reduction(+:totalAfter)
    for (int s = 0; s < NUM_SHARDS; s++) {
        for (auto it = invShards[s].begin(); it != invShards[s].end(); ) {
            if (it->second.size() <= 1) it = invShards[s].erase(it);
            else { ++it; totalAfter++; }
        }
    }
    std::cerr << "singleton removal: " << totalBefore << " -> " << totalAfter
              << " (" << (totalBefore - totalAfter) << " removed), "
              << get_sec() - t2 << " s\n";

    // ── Flatten posting lists into contiguous CSR array ──────────────────────
    double tFlat = get_sec();
    size_t totalPostings = 0;
    for (int s = 0; s < NUM_SHARDS; s++)
        for (auto& [k, v] : invShards[s])
            totalPostings += v.size();
    std::cerr << "total non-singleton postings: " << totalPostings
              << " (" << totalPostings * 4.0 / (1ULL << 30) << " GB)\n";

    const uint32_t HASH_SPACE = use64 ? 0
        : (1u << std::min(4 * (P.half_k - P.drlevel), 28));
    std::vector<size_t>   csrOff;
    std::vector<uint32_t> csrPosts;

    struct PostRange { size_t off; uint32_t cnt; };
    phmap::flat_hash_map<uint64_t, PostRange> postIdx64;

    csrPosts.resize(totalPostings);

    if (!use64 && HASH_SPACE > 0) {
        csrOff.resize(static_cast<size_t>(HASH_SPACE) + 1, 0);
        for (int s = 0; s < NUM_SHARDS; s++)
            for (auto& [k, v] : invShards[s])
                csrOff[static_cast<uint32_t>(k) + 1] = v.size();
        for (uint32_t h = 0; h < HASH_SPACE; h++)
            csrOff[h + 1] += csrOff[h];
        for (int s = 0; s < NUM_SHARDS; s++) {
            for (auto& [k, v] : invShards[s]) {
                uint32_t h = static_cast<uint32_t>(k);
                std::memcpy(&csrPosts[csrOff[h]], v.data(),
                            v.size() * sizeof(uint32_t));
            }
            invShards[s].clear();
        }
    } else {
        postIdx64.reserve(totalAfter);
        size_t pos = 0;
        for (int s = 0; s < NUM_SHARDS; s++) {
            for (auto& [k, v] : invShards[s]) {
                std::memcpy(&csrPosts[pos], v.data(),
                            v.size() * sizeof(uint32_t));
                postIdx64[k] = {pos, static_cast<uint32_t>(v.size())};
                pos += v.size();
            }
            invShards[s].clear();
        }
    }
    invShards.clear(); invShards.shrink_to_fit();
    std::cerr << "flatten CSR: " << get_sec() - tFlat << " s"
              << (use64 ? " (hash map)" : " (direct array)") << "\n";

    // Sort each genome's sketch for sequential CSR access
    double tSort = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(dynamic)
    for (int i = 0; i < N; i++) {
        if (use64) std::sort(sketches[i].h64.begin(), sketches[i].h64.end());
        else       std::sort(sketches[i].h32.begin(), sketches[i].h32.end());
    }
    std::cerr << "sort sketches: " << get_sec() - tSort << " s\n";

    // ── Phase 3: Distance via flat CSR index ─────────────────────────────────
    double t3 = get_sec();

    const double inv_kmer_size = 1.0 / kmerSize;
    const double p_exp   = std::exp(-static_cast<double>(kmerSize) * a.maxDist);
    const double minJac  = p_exp / (2.0 - p_exp);
    const double radio   = 2.0 * std::exp(a.maxDist * (kmerSize - 1)) - 1.0;
    std::cerr << "pruning: minJac=" << minJac << "  radio=" << radio
              << "  (mashD<" << a.maxDist << ", k=" << kmerSize << ")\n";

    std::vector<int> sketchSz(N);
    for (int i = 0; i < N; i++)
        sketchSz[i] = use64 ? static_cast<int>(sketches[i].h64.size())
                             : static_cast<int>(sketches[i].h32.size());

    const size_t*   csrOffPtr  = csrOff.data();
    const uint32_t* csrPostPtr = csrPosts.data();

    FILE* fout = std::fopen(a.output.c_str(), "w");
    std::setvbuf(fout, nullptr, _IOFBF, 1 << 24);
    std::fprintf(fout, "genome0\tgenome1\tcommon|size0|size1\tjaccard\tmashD\n");

    #pragma omp parallel num_threads(a.threads)
    {
        // Kssd sketches can hold > 65535 unique k-mers per genome, so
        // intersection counts can exceed uint16_t.  Keep int here.
        std::vector<int> isect(N, 0);
        std::vector<int> stamp(N, 0);
        int ep = 0;
        std::vector<int> cand;
        cand.reserve(4096);
        std::string buf;
        buf.reserve(1 << 24);

        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < N; i++) {
            const int s0 = sketchSz[i];
            if (__builtin_expect(s0 == 0, 0)) continue;

            cand.clear();
            ++ep;
            if (__builtin_expect(ep == INT_MAX, 0)) {
                std::memset(stamp.data(), 0, N * sizeof(int));
                ep = 1;
            }

            if (use64) {
                const auto& harr = sketches[i].h64;
                const size_t hsz = harr.size();
                for (size_t idx = 0; idx < hsz; idx++) {
                    if (__builtin_expect(idx + 1 < hsz, 1))
                        __builtin_prefetch(&harr[idx + 1], 0, 1);
                    const uint64_t hv = harr[idx];
                    auto it = postIdx64.find(hv);
                    if (__builtin_expect(it == postIdx64.end(), 0)) continue;
                    const uint32_t* pl = csrPostPtr + it->second.off;
                    const uint32_t plSz = it->second.cnt;
                    for (uint32_t ki = 0; ki < plSz; ki++) {
                        uint32_t j = pl[ki];
                        if (static_cast<int>(j) <= i) continue;
                        if (__builtin_expect(stamp[j] != ep, 1)) {
                            stamp[j] = ep;
                            isect[j] = 1;
                            cand.push_back(j);
                        } else {
                            isect[j]++;
                        }
                    }
                }
            } else {
                const auto& harr = sketches[i].h32;
                const size_t hsz = harr.size();
                for (size_t idx = 0; idx < hsz; idx++) {
                    if (__builtin_expect(idx + 1 < hsz, 1))
                        __builtin_prefetch(&harr[idx + 1], 0, 1);
                    const uint32_t hv = harr[idx];
                    const size_t plStart = csrOffPtr[hv];
                    const size_t plEnd   = csrOffPtr[hv + 1];
                    if (__builtin_expect(plStart == plEnd, 0)) continue;
                    for (size_t ki = plStart; ki < plEnd; ki++) {
                        uint32_t j = csrPostPtr[ki];
                        if (static_cast<int>(j) <= i) continue;
                        if (__builtin_expect(stamp[j] != ep, 1)) {
                            stamp[j] = ep;
                            isect[j] = 1;
                            cand.push_back(j);
                        } else {
                            isect[j]++;
                        }
                    }
                }
            }

            const int minCommon = std::max(1,
                static_cast<int>(std::ceil(minJac * s0)));

            for (int j : cand) {
                const int common = isect[j];
                if (common < minCommon) continue;

                const int s1 = sketchSz[j];
                const int mn = s0 < s1 ? s0 : s1;
                const int mx = s0 > s1 ? s0 : s1;
                if (__builtin_expect(mx > radio * mn, 0)) continue;

                const int denom = s0 + s1 - common;
                const double jac = static_cast<double>(common) / denom;
                const double mashD = (jac >= 1.0) ? 0.0
                    : -inv_kmer_size * std::log(2.0 * jac / (1.0 + jac));
                if (mashD < a.maxDist) {
                    char line[1024];
                    int n = std::snprintf(line, sizeof(line),
                        "%s\t%s\t%d|%d|%d\t%.6f\t%.6f\n",
                        files[j].c_str(), files[i].c_str(),
                        common, s0, s1, jac, mashD);
                    buf.append(line, n);
                }
            }

            if (buf.size() > (1 << 24)) {
                #pragma omp critical
                { std::fwrite(buf.data(), 1, buf.size(), fout); }
                buf.clear();
            }
            if (i % progress == 0)
                std::cerr << "  dist " << i << " / " << N << "\n";
        }
        if (!buf.empty()) {
            #pragma omp critical
            { std::fwrite(buf.data(), 1, buf.size(), fout); }
        }
    }
    std::fclose(fout);

    double t4 = get_sec();
    std::cerr << "dist time: " << t4 - t3 << " s\n";
    std::cerr << "total: "     << t4 - t0 << " s\n";

    // Skip destructors for sketches (208K vectors) and large CSR which would
    // otherwise take 10-30s to free; Sketch::computeDistances does the same.
    std::cerr.flush();
    std::_Exit(0);
}

// ═══════════════════════════════════════════════════════════════════════════
//  §9  BinDash all-pairs  (mirrors test_BinDash.cpp – flat_usigs storage,
//      sort-by-cardinality + cycle-following permutation, 128KB tile-based
//      dist with size-ratio prune, per-thread temp files merged at the end)
// ═══════════════════════════════════════════════════════════════════════════
//
// Note: cleanup-on-signal is intentionally omitted in the CLI version (where
// we expect controlled shutdown).  The test binary registers SIGTERM/SIGINT
// handlers; the CLI relies on the OS to clean its temp files in the worst case.
static void run_allpairs_bindash(const Args& a, const std::vector<std::string>& files)
{
    using namespace std;

    const uint32_t SKETCH64 = a.bdSketch64;
    const int      KSIZE    = a.kmerSize;
    const uint32_t BBITS    = a.bdBbits;
    const uint64_t SEED     = a.seed;
    const uint32_t NBINS    = SKETCH64 * 64;
    const uint32_t NWORDS   = SKETCH64 * BBITS;

    int N = static_cast<int>(files.size());
    std::vector<std::string> fileList = files;  // we mutate during permutation

    if (N == 0) {
        std::cerr << "no input files found\n";
        return;
    }

    const int actualThreads = std::min(a.threads, N);

    std::vector<uint64_t> flat_usigs(static_cast<size_t>(N) * NWORDS, 0);
    std::vector<double>   sizes(N, 0.0);
    std::vector<uint8_t>  valid(N, 0);

    double t0 = get_sec();

    // ── Phase 1: sketch + extract into flat array ──────────────────────────
    #pragma omp parallel for num_threads(actualThreads) schedule(dynamic)
    for (int t = 0; t < N; t++) {
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) continue;
        kseq_t* ks = kseq_init(fp);

        Sketch::BinDash sk(SKETCH64, KSIZE, BBITS, SEED);
        while (kseq_read(ks) >= 0)
            sk.update(ks->seq.s, ks->seq.l);
        sk.finalize();

        kseq_destroy(ks);
        gzclose(fp);

        std::memcpy(&flat_usigs[static_cast<size_t>(t) * NWORDS],
                    sk.getSignatures(), NWORDS * sizeof(uint64_t));

        const double ne = static_cast<double>(sk.getRawNonempty());
        sizes[t] = (ne <= 0.0) ? 0.0
                 : (ne >= NBINS) ? 1e18
                 : -static_cast<double>(NBINS) * std::log1p(-ne / NBINS);
        valid[t] = 1;
    }

    double t1 = get_sec();
    std::cerr << "sketch + finalize: " << t1 - t0 << " s\n";

    // ── Sort by cardinality descending ──────────────────────────────────────
    // After sorting: sizes[i] >= sizes[j] for i < j.  Inner j loop can break
    // when sizes[j] / sizes[i] < minJac (impossible to meet Jaccard threshold
    // given containment bound min(|A|,|B|)/max(|A|,|B|) < minJac).
    std::vector<int> order(N);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a_, int b_) { return sizes[a_] > sizes[b_]; });

    // In-place permutation via cycle following – uses one 2KB sketch buffer
    // instead of an N × NWORDS copy; saves ~428 MB peak RSS at N=208K.
    {
        std::vector<uint64_t> tmp(NWORDS);
        std::vector<bool> done(N, false);
        for (int i = 0; i < N; ++i) {
            if (done[i] || order[i] == i) { done[i] = true; continue; }
            std::memcpy(tmp.data(), &flat_usigs[static_cast<size_t>(i) * NWORDS],
                        NWORDS * 8);
            double      tmp_sz = sizes[i];
            std::string tmp_f  = std::move(fileList[i]);
            uint8_t     tmp_v  = valid[i];
            int j = i;
            int kk = order[j];
            while (kk != i) {
                std::memcpy(&flat_usigs[static_cast<size_t>(j) * NWORDS],
                            &flat_usigs[static_cast<size_t>(kk) * NWORDS],
                            NWORDS * 8);
                sizes[j]    = sizes[kk];
                fileList[j] = std::move(fileList[kk]);
                valid[j]    = valid[kk];
                done[j]     = true;
                j  = kk;
                kk = order[j];
            }
            std::memcpy(&flat_usigs[static_cast<size_t>(j) * NWORDS],
                        tmp.data(), NWORDS * 8);
            sizes[j]    = tmp_sz;
            fileList[j] = std::move(tmp_f);
            valid[j]    = tmp_v;
            done[j]     = true;
        }
    }
    { std::vector<int>().swap(order); }

    // ── Resolve output path (handle directory) ─────────────────────────────
    std::string finalPath = a.output;
    {
        struct stat st;
        if (stat(finalPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            if (finalPath.back() != '/') finalPath += '/';
            finalPath += "rabbitsketch.dist";
        }
    }

    // ── Pruning thresholds ─────────────────────────────────────────────────
    const double p_exp   = std::exp(-static_cast<double>(KSIZE) * a.maxDist);
    const double minJac  = p_exp / (2.0 - p_exp);
    const double p_rand  = 1.0 / static_cast<double>(1ULL << BBITS);
    const uint64_t minSamebits = static_cast<uint64_t>(std::max(1,
        static_cast<int>(std::ceil(NBINS * (minJac * (1.0 - p_rand) + p_rand)))));
    std::cerr << "pruning: minJac=" << minJac
              << "  minSamebits=" << minSamebits << "/" << NBINS
              << "  (mashD<" << a.maxDist << ", k=" << KSIZE << ")\n";

    // ── Per-thread temp files ──────────────────────────────────────────────
    const int procId = static_cast<int>(getpid());
    std::vector<std::string> partPaths(actualThreads);
    for (int tid = 0; tid < actualThreads; ++tid)
        partPaths[tid] = finalPath + ".part." + std::to_string(procId)
                                   + "." + std::to_string(tid);

    double t2 = get_sec();

    static const int TILE = 64;  // 64 sketches × 2KB = 128KB → fits in L2
    const int nTiles = (N + TILE - 1) / TILE;
    int progress = N / 20;
    if (progress < 1) progress = 1;

    // ── Phase 2: tiled distance computation ───────────────────────────────
    #pragma omp parallel num_threads(actualThreads)
    {
        const int tid = omp_get_thread_num();
        FILE* tf = std::fopen(partPaths[tid].c_str(), "wb");
        if (!tf) {
            #pragma omp critical
            std::cerr << "ERROR: cannot open temp file: " << partPaths[tid] << "\n";
        } else {
            std::setvbuf(tf, nullptr, _IOFBF, 1 << 22);
            std::string buf;
            buf.reserve(1 << 22);

            #pragma omp for schedule(dynamic, 1)
            for (int tileI = 0; tileI < nTiles; ++tileI) {
                const int iBeg = tileI * TILE;
                const int iEnd = std::min(N, iBeg + TILE);

                for (int tileJ = tileI; tileJ < nTiles; ++tileJ) {
                    const int jBeg = tileJ * TILE;
                    const int jEnd = std::min(N, jBeg + TILE);

                    if (tileJ > tileI && sizes[iEnd - 1] > 0.0 &&
                        sizes[jBeg] / sizes[iEnd - 1] < minJac)
                        break;

                    for (int i = iBeg; i < iEnd; ++i) {
                        if (!valid[i]) continue;
                        const uint64_t* ui =
                            &flat_usigs[static_cast<size_t>(i) * NWORDS];
                        const double si = sizes[i];

                        const int jStart = (tileI == tileJ) ? std::max(i + 1, jBeg) : jBeg;

                        for (int j = jStart; j < jEnd; ++j) {
                            if (!valid[j]) continue;

                            if (si > 0.0 && sizes[j] / si < minJac) break;

                            if (j + 4 < jEnd)
                                __builtin_prefetch(
                                    &flat_usigs[static_cast<size_t>(j + 4) * NWORDS],
                                    0, 1);

                            const uint64_t* uj =
                                &flat_usigs[static_cast<size_t>(j) * NWORDS];

                            const uint64_t same = Sketch::BinDash::countSameBits(
                                ui, uj, SKETCH64, BBITS);
                            if (same < minSamebits) continue;

                            const double p_match = static_cast<double>(same) / NBINS;
                            const double jac     = (p_match - p_rand) / (1.0 - p_rand);
                            const double dist    = (jac >= 1.0) ? 0.0
                                : -std::log(2.0 * jac / (1.0 + jac))
                                  / static_cast<double>(KSIZE);

                            if (dist < a.maxDist) {
                                char line[1024];
                                int len = std::snprintf(line, sizeof(line),
                                    "%s\t%s\t%.6f\n",
                                    fileList[i].c_str(), fileList[j].c_str(), dist);
                                buf.append(line, static_cast<size_t>(len));
                            }
                        }
                    }

                    if (buf.size() > (1 << 22)) {
                        std::fwrite(buf.data(), 1, buf.size(), tf);
                        buf.clear();
                    }
                }

                const int iProgress = std::min(iEnd, N);
                if (iProgress % progress == 0) {
                    #pragma omp critical
                    std::cerr << "  dist " << iProgress << " / " << N << "\n";
                }
            }
            if (!buf.empty()) std::fwrite(buf.data(), 1, buf.size(), tf);
            std::fclose(tf);
        }
    }

    // ── Merge temp files ───────────────────────────────────────────────────
    auto cleanup_parts = [&]() {
        for (int tid = 0; tid < actualThreads; ++tid)
            std::remove(partPaths[tid].c_str());
    };

    FILE* fout = std::fopen(finalPath.c_str(), "wb");
    if (!fout) {
        std::cerr << "ERROR: cannot open output: " << finalPath << "\n";
        cleanup_parts();
        return;
    }
    std::setvbuf(fout, nullptr, _IOFBF, 1 << 24);
    std::vector<char> copyBuf(1 << 20);
    for (int tid = 0; tid < actualThreads; ++tid) {
        FILE* part = std::fopen(partPaths[tid].c_str(), "rb");
        if (!part) continue;
        while (true) {
            size_t got = std::fread(copyBuf.data(), 1, copyBuf.size(), part);
            if (!got) break;
            std::fwrite(copyBuf.data(), 1, got, fout);
        }
        std::fclose(part);
        std::remove(partPaths[tid].c_str());
    }
    std::fclose(fout);

    double t3 = get_sec();
    std::cerr << "dist time: "  << t3 - t2 << " s\n";
    std::cerr << "total time: " << t3 - t0 << " s\n";
    std::cerr << "output: "     << finalPath << "\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  §9b  HLL --index  (mirrors test_HLL.cpp's NEW inverted-index path)
//
//      Replaces the old LSH band-streaming approach (which was lossy due to
//      MAX_BUCKET skips) with a SetSketch-style witness inverted index:
//
//      1. Each HyperLogLog tracks the 64-bit hash that "won" each register
//         (i.e., the k-mer with the longest leading-zero run that mapped to
//         that bucket). Sharing a register-i witness ⇒ shared k-mer ∈ A∩B.
//      2. Subsample witnesses with WITNESS_STRIDE=4 (one key per 4 registers).
//         For bits=13 this gives 2048 keys per sketch (similar to FastKMV K).
//      3. Sort by cardinality descending so size-ratio pruning can `break`
//         instead of `continue`.
//      4. Build CSR inverted index, singleton-filter keys, then verify
//         candidates exactly with HyperLogLog::distance() (Ertl joint MLE).
//
//      Mode A (expectedOverlap >= 20) → inverted index + exact verify.
//      Mode B fallback → in-memory O(N²) all-pairs over the same HLL vector
//                        (no resketch; cheap because cores are 8KB).
// ═══════════════════════════════════════════════════════════════════════════
static void run_index_hll(const Args& a, const std::vector<std::string>& files_in) {
    using namespace std;

    const int N = static_cast<int>(files_in.size());
    if (N == 0) { cerr << "ERROR: empty file list\n"; return; }

    const int bits = a.hllBits;
    const int m    = 1 << bits;

    // One inverted-index key per occupied register: key = p * 64 + core_[p]
    // P(key_A[p] == key_B[p]) ≈ J_HLL  ← same property as MinHash keys.
    // Eliminates witness hash arrays; no witness subsampling needed.
    cerr << "registers=" << m << "  keyEncoding=(reg_idx*64+core_val)"
         << "  hllBits=" << bits << "\n";

    // ── Phase 1: Sketch construction (with witness tracking) ────────────────
    vector<unique_ptr<Sketch::HyperLogLog>> vhll(N);
    vector<double>           sizes(N, 0.0);
    vector<vector<uint64_t>> skKeys(N);
    vector<string>           fileList = files_in;

    double t0 = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(dynamic)
    for (int t = 0; t < N; ++t) {
        auto sk = std::make_unique<Sketch::HyperLogLog>(bits, /*track_witnesses=*/false, a.kmerSize);
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) { vhll[t] = std::move(sk); continue; }
        kseq_t* ks = kseq_init(fp);
        while (kseq_read(ks) >= 0) sk->update(ks->seq.s);
        kseq_destroy(ks);
        gzclose(fp);

        sizes[t] = sk->cardinality();
        const uint8_t* core = sk->getCore().data();
        skKeys[t].reserve(m);
        for (int p = 0; p < m; ++p)
            if (core[p] != 0)
                skKeys[t].push_back(static_cast<uint64_t>(p) * 64ULL + core[p]);

        vhll[t] = std::move(sk);
    }
    double t1 = get_sec();
    cerr << "sketch time: " << t1 - t0 << " s\n";

    // ── Sort by cardinality descending + reorder all per-genome arrays ──────
    {
        vector<int> order(N);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(),
             [&](int aa, int bb) { return sizes[aa] > sizes[bb]; });

        vector<unique_ptr<Sketch::HyperLogLog>> sH(N);
        vector<string>                          sF(N);
        vector<double>                          sS(N, 0.0);
        vector<vector<uint64_t>>                sK(N);
        for (int ni = 0; ni < N; ++ni) {
            const int oi = order[ni];
            sH[ni] = std::move(vhll[oi]);
            sF[ni] = fileList[oi];
            sS[ni] = sizes[oi];
            sK[ni] = std::move(skKeys[oi]);
        }
        vhll.swap(sH);
        fileList.swap(sF);
        sizes.swap(sS);
        skKeys.swap(sK);
    }

    // ── Decide mode: inverted index needs ≥20 expected register-value overlap ──
    // P(core_A[p] == core_B[p]) ≈ J_HLL, so expected matching keys = m * minJac.
    // expectedOverlap >= 20 ⇒ P(missed valid pair) < e^(-20).
    const double p_exp  = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const double expectedOverlap = static_cast<double>(m) * minJac;
    const bool   useInvIdx = (expectedOverlap >= 20.0);
    cerr << "minJac=" << minJac << "  expectedRegValOverlap=" << expectedOverlap
         << "  mashD<" << a.maxDist << "  k=" << a.kmerSize << "\n";

    double t5 = get_sec();

    if (useInvIdx) {
        // ───────────── MODE A: inverted index (witness hashes) + exact ─────
        cerr << "mode: INVERTED INDEX (HLL register-value keys)\n";

        const int actualThreads = min(a.threads, N);
        vector<phmap::flat_hash_map<uint64_t, vector<uint32_t>>> threadIdx(actualThreads);

        double ti0 = get_sec();
        #pragma omp parallel num_threads(a.threads)
        {
            int tid = omp_get_thread_num();
            if (tid < actualThreads) {
                auto& localIdx = threadIdx[tid];
                #pragma omp for schedule(static)
                for (int t = 0; t < N; ++t)
                    for (uint64_t key : skKeys[t])
                        localIdx[key].push_back(static_cast<uint32_t>(t));
            }
        }
        double ti1 = get_sec();
        cerr << "local index: " << ti1 - ti0 << " s\n";

        auto csrIdx = Sketch::buildCSRIndex<uint64_t>(threadIdx, a.threads);
        double ti2 = get_sec();
        cerr << "CSR build: " << ti2 - ti1 << " s\n";

        // Singleton-filter skKeys: keys present in only one sketch can't
        // produce candidates. Same optimization as run_index_setsketch.
        {
            size_t kbefore = 0, kafter = 0;
            #pragma omp parallel for num_threads(a.threads) schedule(dynamic, 64) \
                                     reduction(+:kbefore,kafter)
            for (int t = 0; t < N; ++t) {
                kbefore += skKeys[t].size();
                vector<uint64_t> kept;
                kept.reserve(skKeys[t].size() / 5);
                for (uint64_t key : skKeys[t])
                    if (csrIdx.postIdx.count(key))
                        kept.push_back(key);
                kafter += kept.size();
                skKeys[t] = std::move(kept);
            }
            cerr << "skKeys singleton-filtered: " << kbefore << " -> " << kafter
                 << " (" << (kbefore > 0 ? 100.0 * kafter / kbefore : 0.0)
                 << "% retained, freed ~"
                 << (kbefore - kafter) * 8 / (1 << 20) << " MB)\n";
        }

        // 8-sigma minCommon (same as SetSketch). Witness sharing follows
        // approximate binomial(witnessesPerSketch, J) with mean expectedOverlap.
        const double sd = std::sqrt(std::max(0.0, expectedOverlap * (1.0 - minJac)));
        const int minCommon = max(1, static_cast<int>(std::floor(expectedOverlap - 8.0 * sd)));
        cerr << "minCommon=" << minCommon << " (expected=" << expectedOverlap
             << ", 8sigma=" << 8.0 * sd << ")\n";

        // Exact verification: HLL::distance() runs Ertl joint MLE on the two
        // register arrays. Size-ratio prune up front saves the MLE call.
        auto exactJaccardFn = [&](int i, int j) -> double {
            const double si = sizes[i], sj = sizes[j];
            if (si > 0.0 && sj / si < minJac) return -1.0;
            const double d = vhll[i]->distance(*vhll[j]);
            return 1.0 - d;  // distance() returns 1 - jaccard
        };
        auto minCommonFn = [minCommon](int) { return minCommon; };

        Sketch::computeDistancesExact<uint64_t>(
            csrIdx, skKeys, fileList,
            N, a.kmerSize, a.maxDist, exactJaccardFn, minCommonFn,
            a.output, a.threads);
    } else {
        // ───────────── MODE B: all-pairs O(N²) over already-built HLLs ─────
        // No re-sketch; reuse vhll. Cheap because HLL distance is one O(m)
        // SIMD pass over two 8KB register arrays. Sorted-by-cardinality
        // already, so size-ratio prune is a `break` not `continue`.
        cerr << "mode: ALL-PAIRS (sorted, size-ratio prune)\n";
        { vector<vector<uint64_t>>().swap(skKeys); }   // free witnesses

        // Resolve output path (handle directory)
        string finalPath = a.output;
        {
            struct stat st;
            if (stat(finalPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
                if (finalPath.back() != '/') finalPath += '/';
                finalPath += "rabbitsketch.dist";
            }
        }
        cerr << "output: " << finalPath << "\n";

        vector<string> bufs(a.threads);
        int progress = N / 20;
        if (progress < 1) progress = 1;

        #pragma omp parallel for num_threads(a.threads) schedule(dynamic, 1)
        for (int i = 0; i < N; ++i) {
            const int tid = omp_get_thread_num();
            string& buf = bufs[tid];
            const double si = sizes[i];
            for (int j = i + 1; j < N; ++j) {
                const double sj = sizes[j];
                if (si > 0.0 && sj / si < minJac) break;  // sorted descending

                const double dist = vhll[i]->distance(*vhll[j]);
                if (dist < a.maxDist) {
                    char line[1024];
                    int len = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                        fileList[i].c_str(), fileList[j].c_str(), dist);
                    buf.append(line, static_cast<size_t>(len));
                }
            }
            if (i % progress == 0) {
                #pragma omp critical
                cerr << "  dist " << i << " / " << N << "\n";
            }
        }

        FILE* fp = fopen(finalPath.c_str(), "w");
        if (!fp) err(errno, "cannot open output: %s", finalPath.c_str());
        setvbuf(fp, nullptr, _IOFBF, 1 << 22);
        for (auto& b : bufs) fwrite(b.data(), 1, b.size(), fp);
        fclose(fp);
    }

    double t6 = get_sec();
    cerr << "dist time: "  << t6 - t5 << " s\n";
    cerr << "total time: " << t6 - t0 << " s\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  §10  CLI dispatch
// ═══════════════════════════════════════════════════════════════════════════
int run_cli(const Args& a) {
    if (!a.listMode) return run_pair(a);

    auto files = load_file_list(a.inputs[0]);
    if (files.empty()) { std::cerr << "ERROR: empty file list\n"; return 1; }

    std::cerr << "===== total files: " << files.size() << "  algo=";
    switch (a.algo) {
        case Algo::MINHASH:     std::cerr << "MinHash";     break;
        case Algo::KSSD:        std::cerr << "Kssd";        break;
        case Algo::HLL:         std::cerr << "HLL";         break;
        case Algo::PROBMINHASH: std::cerr << "ProbMinHash"; break;
        case Algo::BINDASH:     std::cerr << "BinDash";     break;
        case Algo::SETSKETCH:   std::cerr << "SetSketch";   break;
        case Algo::FASTKMV:     std::cerr << "FastKMV";     break;
        default: break;
    }
    std::cerr << (a.useIndex ? "  (inverted index)\n" : "  (all-pairs)\n");

    if (a.useIndex) {
        switch (a.algo) {
            case Algo::MINHASH:     run_index_minhash(a, files);    return 0;
            case Algo::FASTKMV:     run_index_fastkmv(a, files);    return 0;
            case Algo::PROBMINHASH: run_index_probmh(a, files);     return 0;
            case Algo::KSSD:        run_index_kssd(a, files);       return 0;  // never returns (uses _Exit)
            case Algo::SETSKETCH:   run_index_setsketch(a, files);  return 0;
            case Algo::HLL:         run_index_hll(a, files);        return 0;
            default:
                std::cerr << "WARNING: --index not supported for this algorithm; "
                             "falling back to all-pairs\n";
                break;
        }
    }

    // All-pairs (no index) – BinDash uses its specialised fast path; the rest
    // share the generic list_allpairs<T> template.
    if (a.algo == Algo::BINDASH) {
        run_allpairs_bindash(a, files);
        return 0;
    }

    const int kmer = a.kmerSize;
    switch (a.algo) {
    case Algo::MINHASH:
        list_allpairs<Sketch::MinHash>(a, files, build_minhash,
            [](Sketch::MinHash* x, Sketch::MinHash* y) { return x->distance(y); });
        break;

    case Algo::HLL:
        list_allpairs<Sketch::HyperLogLog>(a, files, build_hll,
            [kmer](Sketch::HyperLogLog* x, Sketch::HyperLogLog* y) {
                return mash_distance(x->jaccard_index(*y), kmer);
            });
        break;

    case Algo::KSSD:
        list_allpairs<Sketch::Kssd>(a, files, build_kssd,
            [](Sketch::Kssd* x, Sketch::Kssd* y) { return x->distance(y); });
        break;

    case Algo::PROBMINHASH:
        list_allpairs<Sketch::ProbMinHash4>(a, files, build_probminhash,
            [kmer](Sketch::ProbMinHash4* x, Sketch::ProbMinHash4* y) {
                return mash_distance(x->jaccard_weighted(*y), kmer);
            });
        break;

    case Algo::SETSKETCH:
        list_allpairs<Sketch::SetSketch>(a, files, build_setsketch,
            [kmer](Sketch::SetSketch* x, Sketch::SetSketch* y) {
                return mash_distance(x->jaccard_index(*y), kmer);
            });
        break;

    case Algo::FASTKMV:
        list_allpairs<Sketch::FastKMV>(a, files, build_fastkmv,
            [kmer](Sketch::FastKMV* x, Sketch::FastKMV* y) {
                return mash_distance(x->jaccard(*y), kmer);
            });
        break;

    default: return 1;
    }

    return 0;
}
