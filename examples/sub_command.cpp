/**
 * sub_command.cpp – every rabbitsketch algorithm pipeline lives here.
 *
 * Layout
 * ──────
 *   §1  Helpers (file list IO, mash distance, seq streaming, KSSD param cache)
 *   §2  Per-algorithm sketch builders (build_*)        — pair / all-pairs mode
 *   §3  run_pair                                       — pairwise distance
 *   §4  list_allpairs<T>                               — generic O(N²)
 *   §5  run_index_fastkmv  (bucket CSR + native bottom-k verification)
 *   §6  run_index_probmh   (bucket CSR + frequency-weighted winner keys)
 *   §7  run_index_setsketch (mirrors test_SetSketch.cpp inverted-index mode)
 *   §8  run_index_kssd     (shared with test_Kssd.cpp: bitmap admission,
 *                          SIMD encoding and contiguous CSR index)
 *   §9  run_allpairs_bindash (mirrors test_BinDash.cpp – flat_usigs + sort by
 *                          cardinality + tile-based dist + per-thread temp files)
 *   §10 run_cli            — dispatch table
 */

#include "sub_command.h"

#include "Sketch.h"
#include "KssdBuild.h"
#include "KssdIndex.h"
#include "io/KssdFastx.h"
#include "BinDash.h"
#include "BinDashIndex.h"
#include "fastkmv.h"
#include "fracminhash.h"
#include "probmh.h"
#include "SetSketch.h"
#include "InvertedIndex.h"
#include "shuffle.h"
#include "common.h"
#include "kseq.h"

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

// ── 64-bit mixing + bottom-K (KMV) keys from witnesses ─────────────────────
// Used as a high-entropy LSH layer for register-based sketches (HLL/SetSketch).
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static inline void keys_from_witnesses_kmv(
    const uint64_t* wit,
    int m,
    int K,
    uint64_t seed,
    std::vector<uint64_t>& out)
{
    out.clear();
    out.reserve(static_cast<size_t>(m));
    for (int p = 0; p < m; ++p) {
        const uint64_t w = wit[p];
        if (w != 0) out.push_back(splitmix64(w ^ seed));
    }
    if (K <= 0 || out.empty()) return;
    if (static_cast<int>(out.size()) > K) {
        std::nth_element(out.begin(), out.begin() + K, out.end());
        out.resize(static_cast<size_t>(K));
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
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

// The library owns the shared immutable parameter cache.
static std::shared_ptr<const Sketch::kssd_parameter_t>
get_kssd_params(int half_k, int half_subk, int drlevel) {
    return Sketch::getKssdParameters(half_k, half_subk, drlevel);
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
    const auto parameters = get_kssd_params(a.kmerSize / 2, 6, a.kssdDrlevel);
    auto sketch = std::make_unique<Sketch::Kssd>(parameters, a.seed);
    sketch->fileName = path;
    Sketch::KssdBuild::FastxReader reader;
    Sketch::KssdBuild::Stream stream(*parameters);
    bool okay;
    if (parameters->half_k - parameters->drlevel > 8) {
        Sketch::KssdBuild::Accumulator<uint64_t> hashes(false);
        hashes.reset();
        okay = reader.read(path, stream, [&](uint64_t hash) { hashes.insert(hash); });
        std::vector<uint64_t> values; hashes.finish(values);
        sketch->loadHashes64(std::move(values));
    } else {
        Sketch::KssdBuild::Accumulator<uint32_t> hashes(false);
        hashes.reset();
        okay = reader.read(path, stream, [&](uint64_t hash) { hashes.insert(hash); });
        std::vector<uint32_t> values; hashes.finish(values);
        sketch->loadHashes(std::move(values));
    }
    if (!okay) throw std::runtime_error("KSSD: cannot read FASTX input: " + path);
    sketch->finalize();
    return sketch.release();
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
    // Drop the per-call SIMD encoding scratch buffer (sized to the largest
    // contig seen).  See FastKMV::finalize().
    sk->finalize();
    return sk;
}

static Sketch::FracMinHash* build_fracminhash(const std::string& path, const Args& a) {
    auto* sk = new Sketch::FracMinHash(a.fracScaled, a.kmerSize, a.seed);
    stream_seq(path, [&](char* seq, uint64_t len) { sk->update(seq, len); });
    sk->finalize();
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
    case Algo::FRACMINHASH: {
        std::unique_ptr<Sketch::FracMinHash> sa(build_fracminhash(fA, a));
        std::unique_ptr<Sketch::FracMinHash> sb(build_fracminhash(fB, a));
        d = sa->distance(*sb);
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
    std::vector<std::vector<uint64_t>> skKeys(N);

    const int actualThreads = std::max(1, std::min(a.threads, N));

    double t0 = get_sec();
    #pragma omp parallel for num_threads(actualThreads) schedule(dynamic)
    for (int t = 0; t < N; ++t) {
        gzFile fp = gzopen(files[t].c_str(), "r");
        if (!fp) continue;
        kseq_t* ks = kseq_init(fp);

        Sketch::MinHash sk(a.kmerSize, a.minhashSize,
                           static_cast<uint32_t>(a.seed), /*rc=*/true);
        while (kseq_read(ks) >= 0) sk.update(ks->seq.s);
        kseq_destroy(ks);
        gzclose(fp);

        skKeys[t] = sk.getHashesSorted(); // calls finalize() internally
    }
    double t1 = get_sec();
    std::cerr << "sketch time: " << t1 - t0 << " s\n";

    auto csrIdx = Sketch::buildCSRIndexFromKeys<uint64_t>(skKeys, N, actualThreads);
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
// PostBuildFn: called once on each sketch pointer after all sketches are built,
// before the dist phase.  Use to release construction-only memory (e.g.
// ProbMinHash4::finalize()) without perturbing the parallel build loop.
// Default is a no-op lambda.
template<class T, class BuildFn, class DistFn,
         class PostBuildFn = void(*)(T*)>
static void list_allpairs(const Args& a,
                          const std::vector<std::string>& files,
                          BuildFn build, DistFn dist,
                          PostBuildFn post_build = [](T*){})
{
    const int N = static_cast<int>(files.size());
    std::vector<std::unique_ptr<T>> sks(N);

    double t0 = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(dynamic)
    for (int i = 0; i < N; ++i) sks[i].reset(build(files[i], a));
    std::cerr << "sketch time: " << get_sec() - t0 << " s\n";

    // Post-build pass: free construction-only memory before the dist phase.
    // Running sequentially is fine – these are just pointer resets (O(N) frees).
    for (int i = 0; i < N; ++i) post_build(sks[i].get());

    double t1 = get_sec();

    // Output streamed directly to disk during the dist loop.  Per-thread
    // 16 MiB buffer; flush via `#pragma omp critical` when full.  This mirrors
    // the test_BinDash / test_Kssd / test_MinHash_OLD pairwise pattern:
    // bounded per-thread RSS during dist, regardless of how many close pairs
    // are found.
    FILE* fp = std::fopen(a.output.c_str(), "w");
    if (!fp) err(errno, "cannot open output: %s", a.output.c_str());
    std::setvbuf(fp, nullptr, _IOFBF, 1 << 22);

    // 16 MiB per-thread buffer + critical flush, matching test_Kssd.cpp.
    // At T threads the dist-phase RSS cap is T × 16 MiB (e.g. 32T → 512 MiB).
    // Larger than test_BinDash.cpp's 4 MiB; trades a bit more RSS for fewer
    // fwrite syscalls and lower omp_critical contention.
    static constexpr size_t FLUSH_BYTES = 1u << 24; // 16 MiB

    #pragma omp parallel num_threads(a.threads)
    {
        std::string buf;
        buf.reserve(FLUSH_BYTES);

        #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                const double d = dist(sks[i].get(), sks[j].get());
                if (d < a.maxDist) {
                    char line[2048];
                    int len = std::snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                                            files[i].c_str(), files[j].c_str(), d);
                    buf.append(line, static_cast<size_t>(len));
                }
            }
            if (buf.size() >= FLUSH_BYTES) {
                #pragma omp critical
                { std::fwrite(buf.data(), 1, buf.size(), fp); }
                buf.clear();
            }
        }
        if (!buf.empty()) {
            #pragma omp critical
            { std::fwrite(buf.data(), 1, buf.size(), fp); }
        }
    }

    std::fclose(fp);
    std::cerr << "dist time:   " << get_sec() - t1 << " s\n";
    std::cerr << "output: " << a.output << "\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  §5  FastKMV --index
//
//      1. Extract each genome's bottom-k hash list.
//      2. Build CSR with the shared parallel bucket scatter/sort builder.
//      3. Use total shared keys as a necessary candidate bound, then verify
//         with the native bottom-k UNION-prefix estimator. The intersection
//         of the two complete sketches divided by K is not that estimator:
//         shared keys above the kth union key must not enter its numerator.
// ═══════════════════════════════════════════════════════════════════════════
static void run_index_fastkmv(const Args& a, const std::vector<std::string>& files) {
    const int N = static_cast<int>(files.size());
    std::vector<std::vector<uint64_t>> skKeys(N);
    std::vector<std::unique_ptr<Sketch::FastKMV>> nativeSketches(N);

    const int actualThreads = std::max(1, std::min(a.threads, N));
    double t0 = get_sec();
    #pragma omp parallel for num_threads(actualThreads) schedule(dynamic)
    for (int t = 0; t < N; ++t) {
        std::unique_ptr<Sketch::FastKMV> sk(build_fastkmv(files[t], a));
        const uint64_t* regs = sk->getRegisters();
        skKeys[t].assign(regs, regs + sk->size());
        nativeSketches[t] = std::move(sk);
    }
    double t1 = get_sec();
    std::cerr << "sketch time: " << t1 - t0 << " s\n";

    auto csrIdx = Sketch::buildCSRIndexFromKeys<uint64_t>(skKeys, N, actualThreads);
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

    auto exactJaccard = [&](int i, int j) -> double {
        return nativeSketches[i]->jaccard(*nativeSketches[j], minJac);
    };
    // For partial sketches the denominator is at least their own size;
    // using K unconditionally would incorrectly drop small identical sets.
    auto minCommonFn = [minJac](int size) {
        return std::max(1, static_cast<int>(std::ceil(minJac * size)));
    };

    double t3 = get_sec();
    Sketch::computeDistancesExact<uint64_t>(csrIdx, skKeys, files,
        N, a.kmerSize, a.maxDist, exactJaccard, minCommonFn, a.output, a.threads);
    double t4 = get_sec();
    std::cerr << "dist time: "  << t4 - t3 << " s\n";
    std::cerr << "total time: " << t4 - t0 << " s\n";
}

// FracMinHash: all sketches share the same sampling threshold, so posting
// counts directly give the native sampled-set intersection and union.
static void run_index_fracminhash(const Args& a, const std::vector<std::string>& files) {
    const int N = static_cast<int>(files.size());
    std::vector<int> sizes(N);
    std::vector<std::vector<uint64_t>> keys(N);
    double t0 = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        std::unique_ptr<Sketch::FracMinHash> sk(build_fracminhash(files[i], a));
        keys[i] = sk->hashes();
        // The shared builder mixes only bucket selection, preserving raw keys.
        sizes[i] = static_cast<int>(keys[i].size());
    }
    const double t1 = get_sec();
    std::cerr << "sketch time: " << t1 - t0 << " s\n";
    auto index = Sketch::buildCSRIndexFromKeys<uint64_t>(keys, N, a.threads);
    const double t2 = get_sec();
    std::cerr << "build CSR index: " << t2 - t1 << " s\n";
    const double minJac = 1.0 / (2.0 * std::exp(a.kmerSize * a.maxDist) - 1.0);
    auto jaccard = [](int common, int left, int right) {
        return static_cast<double>(common) / (static_cast<double>(left) + right - common);
    };
    auto minCommon = [minJac](int size) {
        return std::max(1, static_cast<int>(std::ceil(minJac * size)));
    };
    const double t3 = get_sec();
    Sketch::computeDistances<uint64_t, decltype(jaccard), decltype(minCommon), uint32_t>(
        index, keys, sizes, files, N, a.kmerSize, a.maxDist,
        jaccard, minCommon, a.output, a.threads);
    // The native class defines J(empty, empty)=1. These pairs have no posting.
    if (a.maxDist > 0) {
        std::vector<int> empty;
        for (int i = 0; i < N; ++i) if (sizes[i] == 0) empty.push_back(i);
        if (empty.size() > 1) {
            FILE* output = std::fopen(a.output.c_str(), "a");
            if (!output) err(errno, "cannot append empty-sketch pairs");
            for (size_t i = 0; i < empty.size(); ++i)
                for (size_t j = i + 1; j < empty.size(); ++j)
                    std::fprintf(output, "%s\t%s\t0.000000\n",
                        files[empty[i]].c_str(), files[empty[j]].c_str());
            std::fclose(output);
        }
    }
    const double t4 = get_sec();
    std::cerr << "dist time: " << t4 - t3 << " s\n";
    std::cerr << "total time: " << t4 - t0 << " s\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  §6  ProbMinHash --index  (frequency-weighted J_P mode)
//
//  Two-pass frequency-weighted pipeline:
//    Pass-1: count k-mer occurrences (1/PMH_SDIV hash-sampling for memory).
//    Pass-2: insert each unique k-mer once with weight = count/total_kmers.
//  Same-register winner matches estimate the native probability Jaccard.
//  Index key = winners_[k] XOR (k × golden-ratio) via getWinnerIndexKeys().
// ═══════════════════════════════════════════════════════════════════════════

static void run_index_probmh(const Args& a, const std::vector<std::string>& files) {
    const int N = static_cast<int>(files.size());
    const int mSize = static_cast<int>(a.pmhM);
    std::vector<int> sketchSizes(N, mSize);
    std::vector<std::vector<uint64_t>> skKeys(N);

    const int actualThreads = std::max(1, std::min(a.threads, N));

    double t0 = get_sec();

    #pragma omp parallel for num_threads(actualThreads) schedule(dynamic)
    for (int t = 0; t < N; ++t) {
        std::unique_ptr<Sketch::ProbMinHash4> sk(
            pmh_build_weighted(files[t], a.pmhM, a.kmerSize, a.seed, a.pmhMaxL));
        sk->getWinnerIndexKeys(skKeys[t]);
    }

    double t1 = get_sec();
    std::cerr << "sketch time: " << t1 - t0 << " s\n";

    auto csrIdx = Sketch::buildCSRIndexFromKeys<uint64_t>(skKeys, N, actualThreads);
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

    // Candidate keys for inverted index:
    //   Use MinHash/KMV over per-register witnesses (64-bit winning k-mer hashes).
    //   This yields K high-entropy keys per sketch, avoiding the small key space
    //   problem of indexing raw (register,value) pairs.
    const int keysPerSketch = std::max(1, a.minhashSize);
    cerr << "registers=" << m << "  witnessKMV_keysPerSketch=" << keysPerSketch << "\n";

    // ── Phase 1: Sketch construction + key extraction ────────────────────────
    vector<double>  sizes(N, 0.0);
    vector<uint8_t> flat_cores(static_cast<size_t>(N) * m, 0);
    vector<vector<uint64_t>> skKeys(N);
    vector<string>  fileList = files_in;     // mutable copy (sorted in place below)

    double t0 = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(dynamic)
    for (int t = 0; t < N; ++t) {
        Sketch::SetSketch sk(a.ssBits, a.ssBase, a.ssA, a.kmerSize,
                             /*track_witnesses=*/true);
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
        const uint64_t* wit = sk.getWitnesses().data();
        keys_from_witnesses_kmv(wit, m, keysPerSketch, a.seed, skKeys[t]);
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

    // ── Decide mode: inverted index needs ≥20 expected key overlap ──────────
    // expectedOverlap = keysPerSketch * minJac  (P(match per sampled reg) ≈ J).
    const double p_exp  = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const double expectedOverlap = static_cast<double>(keysPerSketch) * minJac;
    bool useInvIdx = (expectedOverlap >= 20.0);
    if (std::getenv("RABBIT_FORCE_ALLPAIRS")) useInvIdx = false;
    if (std::getenv("RABBIT_FORCE_INDEX")) useInvIdx = true;
    cerr << "minJac=" << minJac << "  expectedRegValOverlap=" << expectedOverlap
         << "  (keysPerSketch=" << keysPerSketch << ")\n";

    double t5 = get_sec();

    if (useInvIdx) {
        // ───────────── MODE A: witness-KMV inverted index + exact verify ─────
        // Approximate candidate generation (may miss); exact verification is unchanged.
        cerr << "mode: INVERTED INDEX (witness-KMV keys)\n";

        double ti0 = get_sec();
        auto csrIdx = Sketch::buildCSRIndexFromKeys<uint64_t>(skKeys, N, a.threads);
        double ti2 = get_sec();
        cerr << "CSR build (radix): " << ti2 - ti0 << " s\n";

        // ── σ-calibrated minCommon (witness-KMV is a Bernoulli LSH) ───────────
        // For a true pair with Jaccard J, # KMV matches ~ Binomial(K, J).
        // At J = minJac:    E = K·minJac,    σ = √(K·minJac·(1-minJac))
        // Use minCommon = max(1, ⌊E − 6σ⌋) → P(false negative | J ≥ minJac) ≤ Φ(−6) ≈ 1e-9.
        // Auto-clamps to 1 when E is small (low-d / small-N), so it falls back
        // to the previous "any-shared-key" behaviour and preserves recall.
        const double sigma     = std::sqrt(expectedOverlap * (1.0 - minJac));
        const int    minCommon = std::max(1,
            static_cast<int>(std::floor(expectedOverlap - 6.0 * sigma)));
        cerr << "minCommon=" << minCommon << "/" << keysPerSketch
             << "  (E=" << expectedOverlap << ", 6sigma=" << 6.0*sigma << ")\n";

        // Two-stage Jaccard: incl-excl SIMD screen (+ early abort), then MLE
        // refinement for surviving candidates. The screen drops ~99% of pairs
        // in typical bio workloads; the MLE step then upgrades the kept pairs
        // to the joint-MLE accuracy (matches `jaccard_index()` semantics).
        const double ssBase = a.ssBase;
        const uint32_t ssQ  = static_cast<uint32_t>(proto.getQ());
        auto exactJaccardFn = [&, ssBase, ssQ](int i, int j) -> double {
            const double si = sizes[i], sj = sizes[j];
            if (si > 0.0 && sj / si < minJac) return -1.0;
            const uint8_t* c1 = &flat_cores[static_cast<size_t>(i) * m];
            const uint8_t* c2 = &flat_cores[static_cast<size_t>(j) * m];
            // Candidates already passed the witness-key inverted index; go
            // straight to MLE.  Skipping the incl-excl pre-screen avoids
            // false negatives at the threshold boundary where incl-excl has a
            // small downward bias that could incorrectly reject valid pairs.
            return Sketch::SetSketch::jaccardFromCoresMLE(
                c1, c2, m, bip, factor, si, sj, ssBase, ssQ);
        };
        // Partial witness sets must not be screened with the full-K count:
        // two identical short inputs can have fewer than minCommon keys.
        auto minCommonFn = [minCommon, keysPerSketch](int actualKeys) {
            return actualKeys < keysPerSketch ? 1 : minCommon;
        };

        Sketch::computeDistancesExact<uint64_t>(
            csrIdx, skKeys, fileList,
            N, a.kmerSize, a.maxDist, exactJaccardFn, minCommonFn,
            a.output, a.threads);

#if 0
        // ─── MODE A: direct-address CSR (register-value keys) + exact verify ──
        //
        // Key space is bounded: key = p*64 + core[p] < m*64 = KEY_SPACE.
        // Build a flat uint32_t CSR array indexed directly by key value instead
        // of a hash map.  Direct array access is ~5-10× faster than phmap
        // probing for this workload (L3 hit vs hash + probe overhead).
        //
        // Build in 3 parallel passes:
        //   P1: atomic count per key
        //   P2: prefix-sum (serial, tiny)
        //   P3: atomic fill into posting array
        cerr << "mode: INVERTED INDEX (direct-address CSR)\n";

        // P1: count
        vector<uint32_t> dcnt(KEY_SPACE, 0);
        {
            // Use per-thread local counts to avoid atomic contention in P1.
            const int nt = min(a.threads, N);
            vector<vector<uint32_t>> lcnt(nt, vector<uint32_t>(KEY_SPACE, 0));
            #pragma omp parallel for num_threads(nt) schedule(dynamic, 4)
            for (int t = 0; t < N; ++t) {
                const int tid = omp_get_thread_num();
                auto& lc = lcnt[tid];
                for (uint64_t key : skKeys[t]) lc[key]++;
            }
            // reduce
            #pragma omp parallel for num_threads(a.threads) schedule(static)
            for (size_t k = 0; k < KEY_SPACE; ++k)
                for (int tid = 0; tid < nt; ++tid) dcnt[k] += lcnt[tid][k];
        }

        // P2: prefix-sum + singleton removal
        vector<size_t> doff(KEY_SPACE + 1, 0);
        size_t totalPosts = 0, uniqueKeys = 0;
        for (size_t k = 0; k < KEY_SPACE; ++k) {
            doff[k] = totalPosts;
            if (dcnt[k] > 1) { totalPosts += dcnt[k]; uniqueKeys++; }
        }
        doff[KEY_SPACE] = totalPosts;
        cerr << "direct CSR: " << uniqueKeys << " unique keys, "
             << totalPosts << " postings\n";

        // P3: fill posting lists (atomic cursor per key)
        vector<uint32_t> dposts(totalPosts);
        {
            // Reuse dcnt as atomic cursor (reset to doff[k] for each key)
            vector<uint32_t> cur(KEY_SPACE);
            for (size_t k = 0; k < KEY_SPACE; ++k)
                cur[k] = static_cast<uint32_t>(doff[k]);

            #pragma omp parallel for num_threads(a.threads) schedule(dynamic, 4)
            for (int t = 0; t < N; ++t) {
                for (uint64_t key : skKeys[t]) {
                    if (doff[key+1] - doff[key] <= 1) continue; // singleton
                    uint32_t pos;
                    #pragma omp atomic capture
                    pos = cur[key]++;
                    dposts[pos] = static_cast<uint32_t>(t);
                }
            }
        }
        { vector<uint32_t>().swap(dcnt); } // free count array

        // IDF prune + non-singleton filter.
        // SetSketch register values cluster around a narrow range → average
        // posting list >> 1 for large N. Keys appearing in more than
        // max(50, expectedOverlap×4) sequences are "stop words": they add
        // candidate explosion without helping discrimination.
        // When N is small (< ~800) expectedOverlap×4 > N → no pruning,
        // which is correct since all-similar small datasets have no "common"
        // register value that can be safely skipped.
        const uint32_t maxPost = static_cast<uint32_t>(
            std::max(50.0, expectedOverlap * 4.0));
        {
            size_t kbefore = 0, kafter = 0;
            #pragma omp parallel for num_threads(a.threads) schedule(dynamic, 4) \
                                     reduction(+:kbefore,kafter)
            for (int t = 0; t < N; ++t) {
                kbefore += skKeys[t].size();
                vector<uint64_t> kept;
                kept.reserve(skKeys[t].size() / 4);
                for (uint64_t key : skKeys[t]) {
                    const uint32_t plN = static_cast<uint32_t>(doff[key+1] - doff[key]);
                    if (plN > 1 && plN <= maxPost) kept.push_back(key);
                }
                kafter += kept.size();
                skKeys[t] = std::move(kept);
            }
            cerr << "skKeys after IDF prune (maxPost=" << maxPost << "): "
                 << kbefore << " -> " << kafter << "\n";
        }

        // Recompute minCommon based on effective (retained) keys per sequence.
        // Use 3σ instead of 8σ: after IDF pruning effExpected is smaller, and
        // for random (unrelated) pairs E[rare_key_match] ≈ 0 so 3σ gives
        // strong precision without sacrificing recall.
        const double effKeys = static_cast<double>(
            [&]{ size_t s=0; for(int t=0;t<N;++t) s+=skKeys[t].size(); return s; }()
        ) / std::max(1, N);
        const double effExpected = effKeys * minJac;
        const double sd = std::sqrt(std::max(0.0, effExpected * (1.0 - minJac)));
        const int minCommon = max(1, static_cast<int>(std::floor(effExpected - 3.0 * sd)));
        cerr << "effKeysPerSeq=" << effKeys
             << "  effExpected=" << effExpected
             << "  minCommon=" << minCommon << "\n";

        // Inline distance computation: direct-address CSR + temp-file output
        {
            const bool useMashDist = (a.kmerSize > 0);
            const double p_exp2 = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
            const double minJacByDist = useMashDist ? p_exp2/(2.0-p_exp2) : 1.0-a.maxDist;

            string finalPath = a.output;
            { struct stat st; if (stat(finalPath.c_str(),&st)==0 && S_ISDIR(st.st_mode)){
                if (finalPath.back()!='/') finalPath+='/';
                finalPath+="rabbitsketch.dist";}}

            const int nt = min(a.threads, N);
            const int procId = static_cast<int>(getpid());
            vector<string> parts(nt);
            for (int tid=0;tid<nt;++tid)
                parts[tid] = finalPath+".part."+to_string(procId)+"."+to_string(tid);

            int progress2 = N/20; if(progress2<1) progress2=1;

            #pragma omp parallel num_threads(nt)
            {
                const int tid = omp_get_thread_num();
                FILE* tf = fopen(parts[tid].c_str(),"wb");
                if(!tf) { cerr<<"ERROR: cannot open "<<parts[tid]<<"\n"; }
                else {
                    setvbuf(tf,nullptr,_IOFBF,1<<22);
                    string buf; buf.reserve(1<<22);
                    vector<uint16_t> isect(N,0);
                    vector<int>      stamp(N,0);
                    int ep=0;
                    vector<int> cand; cand.reserve(4096);

                    #pragma omp for schedule(dynamic,4)
                    for (int i=0;i<N;++i) {
                        const double si=sizes[i];
                        cand.clear(); ++ep;
                        if(__builtin_expect(ep==INT_MAX,0)){
                            memset(stamp.data(),0,N*sizeof(int));ep=1;}

                        for (uint64_t key : skKeys[i]) {
                            const uint32_t* pl  = dposts.data()+doff[key];
                            const uint32_t  plN = static_cast<uint32_t>(doff[key+1]-doff[key]);
                            for (uint32_t pi=0;pi<plN;++pi){
                                int j=static_cast<int>(pl[pi]);
                                if(j<=i) continue;
                                if(__builtin_expect(stamp[j]!=ep,1)){
                                    stamp[j]=ep;isect[j]=1;cand.push_back(j);
                                } else { ++isect[j]; }
                            }
                        }

                        for (int j:cand){
                            if(isect[j]<static_cast<uint16_t>(minCommon)) continue;
                            const double sj=sizes[j];
                            if(si>0.0 && sj/si<minJac) continue;
                            const double jac=Sketch::SetSketch::jaccardFromCoresBatch(
                                &flat_cores[static_cast<size_t>(i)*m],
                                &flat_cores[static_cast<size_t>(j)*m],
                                m,bip,factor,si,sj,minJac,
                                &tailSums[static_cast<size_t>(i)*nCP],
                                &tailSums[static_cast<size_t>(j)*nCP],
                                TAIL_STEP);
                            if(jac<minJacByDist) continue;
                            const double dist=(jac>=1.0)?0.0
                                :-std::log(2.0*jac/(1.0+jac))/static_cast<double>(a.kmerSize);
                            if(dist<a.maxDist){
                                char line[1024];
                                int len=snprintf(line,sizeof(line),"%s\t%s\t%.6f\n",
                                    fileList[i].c_str(),fileList[j].c_str(),dist);
                                buf.append(line,static_cast<size_t>(len));
                            }
                            if(buf.size()>(1<<22)){fwrite(buf.data(),1,buf.size(),tf);buf.clear();}
                        }
                        if(i%progress2==0){
                            #pragma omp critical
                            cerr<<"  dist "<<i<<" / "<<N<<"\n";
                        }
                    }
                    if(!buf.empty()) fwrite(buf.data(),1,buf.size(),tf);
                    fclose(tf);
                }
            }
            FILE* fout=fopen(finalPath.c_str(),"wb");
            if(fout){
                setvbuf(fout,nullptr,_IOFBF,1<<24);
                vector<char> cbuf(1<<20);
                for(int tid=0;tid<nt;++tid){
                    FILE* pf=fopen(parts[tid].c_str(),"rb");
                    if(!pf) continue;
                    while(true){size_t got=fread(cbuf.data(),1,cbuf.size(),pf);if(!got)break;fwrite(cbuf.data(),1,got,fout);}
                    fclose(pf); remove(parts[tid].c_str());
                }
                fclose(fout);
            }
            cerr<<"output: "<<finalPath<<"\n";
        }
#endif
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

        // SetSketch base / q for the joint-MLE refinement step.
        const double   ssBase = a.ssBase;
        const uint32_t ssQ    = static_cast<uint32_t>(proto.getQ());

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

                                // Stage 1: SIMD incl-excl screen with early abort.
                                double jac = Sketch::SetSketch::jaccardFromCoresBatch(
                                    core_i, c2, m, bip, factor, si, sj, minJac,
                                    ts_i, ts_j, TAIL_STEP);
                                if (jac < minJac) continue;

                                // Stage 2: joint-MLE refinement on the survivor.
                                jac = Sketch::SetSketch::jaccardFromCoresMLE(
                                    core_i, c2, m, bip, factor, si, sj, ssBase, ssQ);
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
//  §8  KSSD indexed comparison (shared SIMD/bitmap construction + contiguous
//      merge + flat CSR posting array + stamp-based dist kernel + 5-col output)
// ═══════════════════════════════════════════════════════════════════════════
static void run_index_kssd(const Args& a, const std::vector<std::string>& files)
{
    using namespace std;
    // Same effective k as the native builder: 2 * floor(k / 2), default 20.
    const auto parameters = get_kssd_params(a.kmerSize / 2, 6, a.kssdDrlevel);
    const auto& P = *parameters;
    const int kmerSize = P.kmer_size;
    const bool use64 = (P.half_k - P.drlevel) > 8;
    const int N = static_cast<int>(files.size());
    const int progress = std::max(1, N / 20);
    const int actualThreads = std::max(1, std::min(a.threads, N));
    std::cerr << "half_k=" << P.half_k << " half_subk=" << P.half_subk
              << " drlevel=" << P.drlevel << " kmer=" << kmerSize
              << " use64=" << use64 << "\n";
    struct GSketch { std::vector<uint32_t> h32; std::vector<uint64_t> h64; };
    std::vector<GSketch> sketches(N);
    std::atomic<bool> failed(false);
    const double t0 = get_sec();
    #pragma omp parallel num_threads(actualThreads)
    {
        Sketch::KssdBuild::FastxReader reader;
        Sketch::KssdBuild::Stream stream(P);
        // Only the active width allocates storage. Workspaces survive files.
        Sketch::KssdBuild::Accumulator<uint32_t> hashes32(false);
        Sketch::KssdBuild::Accumulator<uint64_t> hashes64(false);
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < N; ++i) {
            bool okay;
            if (use64) {
                hashes64.reset();
                okay = reader.read(files[i], stream, [&](uint64_t h) { hashes64.insert(h); });
                hashes64.finish(sketches[i].h64);
            } else {
                hashes32.reset();
                okay = reader.read(files[i], stream, [&](uint64_t h) { hashes32.insert(h); });
                hashes32.finish(sketches[i].h32);
            }
            if (!okay) {
                failed.store(true, std::memory_order_relaxed);
                #pragma omp critical
                std::cerr << "KSSD: cannot read FASTX input: " << files[i] << "\n";
            }
            if (i % progress == 0)
                std::cerr << "  sketch " << i << " / " << N << "\n";
        }
    }
    if (failed.load()) throw std::runtime_error("KSSD input read failed");
    const double t1 = get_sec();
    std::cerr << "sketch time: " << t1 - t0 << " s\n";
    const unsigned keyBits = 4 * (P.half_k - P.drlevel);
    auto index = use64
        ? Sketch::KssdBuild::buildIndex<uint64_t>(N, actualThreads, keyBits,
            [&](int i) -> const std::vector<uint64_t>& { return sketches[i].h64; })
        : Sketch::KssdBuild::buildIndex<uint32_t>(N, actualThreads, keyBits,
            [&](int i) -> const std::vector<uint32_t>& { return sketches[i].h32; });
    auto& csrOff = index.offsets;
    auto& csrPosts = index.postings;
    auto& postIdx64 = index.sparse;
    const double t2 = get_sec();
    std::cerr << "build CSR index: " << t2 - t1 << " s\n"
              << "unique hashes: " << index.unique_before << " -> " << index.unique_after << "\n"
              << "total non-singleton postings: " << csrPosts.size() << "\n"
              << "index offsets: " << (csrOff.empty() ? "sparse" : "direct array") << "\n";

    // ── Phase 3: Distance via flat CSR index ─────────────────────────────────
    double t3 = get_sec();

    const double inv_kmer_size = 1.0 / kmerSize;
    const double p_exp   = std::exp(-static_cast<double>(kmerSize) * a.maxDist);
    const double minJac  = p_exp / (2.0 - p_exp);
    const double radio   = 2.0 * std::exp(a.maxDist * kmerSize) - 1.0;
    std::cerr << "pruning: minJac=" << minJac << "  radio=" << radio
              << "  (mashD<" << a.maxDist << ", k=" << kmerSize << ")\n";

    std::vector<int> sketchSz(N);
    for (int i = 0; i < N; i++)
        sketchSz[i] = use64 ? static_cast<int>(sketches[i].h64.size())
                             : static_cast<int>(sketches[i].h32.size());

    const size_t*   csrOffPtr  = csrOff.data();
    const uint32_t* csrPostPtr = csrPosts.data();

    FILE* fout = std::fopen(a.output.c_str(), "w");
    if (!fout) throw std::runtime_error("KSSD: cannot open output");
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

            if (csrOff.empty()) {
                const size_t hsz = use64 ? sketches[i].h64.size() : sketches[i].h32.size();
                for (size_t idx = 0; idx < hsz; idx++) {
                    const uint64_t hv = use64 ? sketches[i].h64[idx] : sketches[i].h32[idx];
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
// BinDash index: exact (bin position, b-bit value) keys. Posting counts equal
// the native packed comparison's matching-bin count, including densified bins.
template<class Key>
static int run_index_bindash_typed(const Args& a, const std::vector<std::string>& files)
{
    namespace BI = Sketch::BinDashIndex;
    const int n = static_cast<int>(files.size());
    const int threads = std::max(1, std::min(a.threads, n));
    const uint32_t bins = a.bdSketch64 * 64;
    std::vector<std::vector<Key>> keys(n);
    std::vector<int> sizes(n, static_cast<int>(bins));
    std::atomic<bool> failed(false);
    const double t0 = get_sec();
    #pragma omp parallel for num_threads(threads) schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        gzFile fp = gzopen(files[i].c_str(), "r");
        if (!fp) { failed.store(true); continue; }
        kseq_t* reader = kseq_init(fp);
        Sketch::BinDash sketch(a.bdSketch64, a.kmerSize, a.bdBbits, a.seed);
        int status;
        while ((status = kseq_read(reader)) >= 0)
            sketch.update(reader->seq.s, reader->seq.l);
        kseq_destroy(reader);
        gzclose(fp);
        if (status < -1) { failed.store(true); continue; }
        keys[i] = BI::keys<Key>(sketch);
    }
    if (failed.load()) {
        std::cerr << "ERROR: BinDash could not read one or more FASTX inputs\n";
        return 1;
    }
    const double t1 = get_sec();
    std::cerr << "sketch + index keys: " << t1 - t0 << " s\n";
    auto index = Sketch::buildCSRIndexFromKeys<Key>(keys, n, threads, BI::Bucket{});
    const double t2 = get_sec();
    std::cerr << "build CSR index: " << t2 - t1 << " s\n";
    const int minimum = static_cast<int>(BI::minimumMatches(
        bins, a.bdBbits, a.kmerSize, a.maxDist));
    std::cerr << "BinDash index: bins=" << bins << " bbits=" << a.bdBbits
              << " key_bytes=" << sizeof(Key) << " minSamebits=" << minimum
              << " (mashD<" << a.maxDist << ")\n";
    auto jaccard = [bins, bits = a.bdBbits](int common, int, int) {
        return BI::jaccard(static_cast<uint32_t>(common), bins, bits);
    };
    auto minCommon = [minimum](int) { return minimum; };
    std::string output = a.output;
    struct stat st;
    if (stat(output.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
        if (output.back() != '/') output += '/';
        output += "rabbitsketch.dist";
    }
    // 32-bit counters also support configurations exceeding 65535 bins.
    Sketch::computeDistances<Key, decltype(jaccard), decltype(minCommon), uint32_t>(
        index, keys, sizes, files, n, a.kmerSize, a.maxDist,
        jaccard, minCommon, output, threads);
    const double t3 = get_sec();
    std::cerr << "dist time: " << t3 - t2 << " s\n"
              << "total time: " << t3 - t0 << " s\n";
    return 0;
}

static int run_index_bindash(const Args& a, const std::vector<std::string>& files)
{
    const uint64_t bins = uint64_t(a.bdSketch64) * 64;
    if (bins == 0 || bins >= INT_MAX || a.bdBbits < 1 || a.bdBbits > 64 ||
        a.kmerSize < 1 || a.kmerSize > 32 || files.size() > INT_MAX) {
        std::cerr << "ERROR: invalid BinDash sketch size, bit width, k-mer size or input count\n";
        return 1;
    }
    const unsigned bits = Sketch::BinDashIndex::keyBits(static_cast<uint32_t>(bins), a.bdBbits);
    if (bits <= 32) return run_index_bindash_typed<uint32_t>(a, files);
    if (bits <= 64) return run_index_bindash_typed<uint64_t>(a, files);
    return run_index_bindash_typed<Sketch::BinDashIndex::WideKey>(a, files);
}

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

    // Candidate keys for inverted index:
    //   Use MinHash/KMV over per-register witnesses (64-bit winning k-mer hashes).
    const int keysPerSketch = std::max(1, a.minhashSize);
    cerr << "registers=" << m << "  witnessKMV_keysPerSketch=" << keysPerSketch
         << "  hllBits=" << bits << "\n";

    // ── Phase 1: Sketch construction ────────────────────────────────────────
    vector<unique_ptr<Sketch::HyperLogLog>> vhll(N);
    vector<double>           sizes(N, 0.0);
    vector<vector<uint64_t>> skKeys(N);
    vector<string>           fileList = files_in;

    double t0 = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(dynamic)
    for (int t = 0; t < N; ++t) {
        auto sk = std::make_unique<Sketch::HyperLogLog>(bits, /*track_witnesses=*/true, a.kmerSize);
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) { vhll[t] = std::move(sk); continue; }
        kseq_t* ks = kseq_init(fp);
        while (kseq_read(ks) >= 0) sk->update(ks->seq.s);
        kseq_destroy(ks);
        gzclose(fp);

        sizes[t] = sk->cardinality();
        const uint64_t* wit = sk->getWitnesses().data();
        keys_from_witnesses_kmv(wit, m, keysPerSketch, a.seed, skKeys[t]);
        // Witnesses no longer needed (verify uses HLL register array only).
        sk->clearWitnesses();

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

    // ── Decide mode: inverted index needs ≥20 expected key overlap ──────────
    // expectedOverlap = keysPerSketch * minJac  (P(match per sampled reg) ≈ J).
    const double p_exp  = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const double expectedOverlap = static_cast<double>(keysPerSketch) * minJac;
    bool useInvIdx = (expectedOverlap >= 20.0);
    if (std::getenv("RABBIT_FORCE_ALLPAIRS")) useInvIdx = false;
    if (std::getenv("RABBIT_FORCE_INDEX")) useInvIdx = true;
    cerr << "minJac=" << minJac << "  expectedRegValOverlap=" << expectedOverlap
         << "  (keysPerSketch=" << keysPerSketch << ")"
         << "  mashD<" << a.maxDist << "  k=" << a.kmerSize << "\n";

    double t5 = get_sec();

    if (useInvIdx) {
        // ───────────── MODE A: witness-KMV inverted index + exact verify ─────
        // Approximate candidate generation (may miss); exact verification uses HLL.
        cerr << "mode: INVERTED INDEX (witness-KMV keys, HLL)\n";

        double ti0 = get_sec();
        auto csrIdx = Sketch::buildCSRIndexFromKeys<uint64_t>(skKeys, N, a.threads);
        double ti2 = get_sec();
        cerr << "CSR build (radix): " << ti2 - ti0 << " s\n";

        // ── σ-calibrated minCommon (see SetSketch path for derivation) ────
        const double sigma     = std::sqrt(expectedOverlap * (1.0 - minJac));
        const int    minCommon = std::max(1,
            static_cast<int>(std::floor(expectedOverlap - 6.0 * sigma)));
        cerr << "minCommon=" << minCommon << "/" << keysPerSketch
             << "  (E=" << expectedOverlap << ", 6sigma=" << 6.0*sigma << ")\n";

        auto exactJaccardFn = [&](int i, int j) -> double {
            const double si = sizes[i], sj = sizes[j];
            if (si > 0.0 && sj / si < minJac) return -1.0;
            return vhll[i]->jaccard_index(*vhll[j]);
        };
        // Partial witness sets must not be screened with the full-K count:
        // two identical short inputs can have fewer than minCommon keys.
        auto minCommonFn = [minCommon, keysPerSketch](int actualKeys) {
            return actualKeys < keysPerSketch ? 1 : minCommon;
        };

        Sketch::computeDistancesExact<uint64_t>(
            csrIdx, skKeys, fileList,
            N, a.kmerSize, a.maxDist, exactJaccardFn, minCommonFn,
            a.output, a.threads);

#if 0
        // ─── MODE A: direct-address CSR (register-value keys) + HLL exact ─────
        cerr << "mode: INVERTED INDEX (direct-address CSR, HLL)\n";

        // P1: per-thread local count to avoid atomic contention
        vector<uint32_t> dcnt(KEY_SPACE, 0);
        {
            const int nt = min(a.threads, N);
            vector<vector<uint32_t>> lcnt(nt, vector<uint32_t>(KEY_SPACE, 0));
            #pragma omp parallel for num_threads(nt) schedule(dynamic, 4)
            for (int t = 0; t < N; ++t) {
                const int tid = omp_get_thread_num();
                auto& lc = lcnt[tid];
                for (uint64_t key : skKeys[t]) lc[key]++;
            }
            #pragma omp parallel for num_threads(a.threads) schedule(static)
            for (size_t k = 0; k < KEY_SPACE; ++k)
                for (int tid = 0; tid < nt; ++tid) dcnt[k] += lcnt[tid][k];
        }

        // P2: prefix-sum + singleton removal
        vector<size_t> doff(KEY_SPACE + 1, 0);
        size_t totalPosts = 0, uniqueKeys = 0;
        for (size_t k = 0; k < KEY_SPACE; ++k) {
            doff[k] = totalPosts;
            if (dcnt[k] > 1) { totalPosts += dcnt[k]; uniqueKeys++; }
        }
        doff[KEY_SPACE] = totalPosts;
        cerr << "direct CSR: " << uniqueKeys << " unique keys, "
             << totalPosts << " postings\n";

        // P3: fill posting lists (atomic cursor per key)
        vector<uint32_t> dposts(totalPosts);
        {
            vector<uint32_t> cur(KEY_SPACE);
            for (size_t k = 0; k < KEY_SPACE; ++k)
                cur[k] = static_cast<uint32_t>(doff[k]);

            #pragma omp parallel for num_threads(a.threads) schedule(dynamic, 4)
            for (int t = 0; t < N; ++t) {
                for (uint64_t key : skKeys[t]) {
                    if (doff[key+1] - doff[key] <= 1) continue;
                    uint32_t pos;
                    #pragma omp atomic capture
                    pos = cur[key]++;
                    dposts[pos] = static_cast<uint32_t>(t);
                }
            }
        }
        { vector<uint32_t>().swap(dcnt); }

        // IDF prune + non-singleton filter on skKeys.
        // HLL register values follow a geometric distribution and cluster around
        // log2(n_p) (n_p = items per bucket). For large N, many (p,v) keys are
        // extremely common (posting lists of ~N/10 – N/5), making candidate
        // generation degenerate to O(N²). Drop keys in > max(50, expectedOverlap×4)
        // sequences (IDF "stop word" elimination). When N is small,
        // expectedOverlap×4 > N → maxPost > N → no pruning (correct for
        // all-similar small datasets where every key is "common").
        const uint32_t maxPost = static_cast<uint32_t>(
            std::max(50.0, expectedOverlap * 4.0));
        {
            size_t kbefore = 0, kafter = 0;
            #pragma omp parallel for num_threads(a.threads) schedule(dynamic, 4) \
                                     reduction(+:kbefore,kafter)
            for (int t = 0; t < N; ++t) {
                kbefore += skKeys[t].size();
                vector<uint64_t> kept;
                kept.reserve(skKeys[t].size() / 4);
                for (uint64_t key : skKeys[t]) {
                    const uint32_t plN = static_cast<uint32_t>(doff[key+1] - doff[key]);
                    if (plN > 1 && plN <= maxPost) kept.push_back(key);
                }
                kafter += kept.size();
                skKeys[t] = std::move(kept);
            }
            cerr << "skKeys after IDF prune (maxPost=" << maxPost << "): "
                 << kbefore << " -> " << kafter << "\n";
        }

        // Recompute effective expectedOverlap and minCommon based on retained keys.
        // effectiveKeysPerSeq = average number of retained keys (approx kafter/N).
        const double effKeys = static_cast<double>(
            [&]{ size_t s=0; for(int t=0;t<N;++t) s+=skKeys[t].size(); return s; }()
        ) / std::max(1, N);
        const double effExpected = effKeys * minJac;
        const double sd = std::sqrt(std::max(0.0, effExpected * (1.0 - minJac)));
        // 3σ: for IDF-pruned rare keys, unrelated pairs have E[match]≈0, so 3σ
        // gives strong precision without sacrificing recall.
        const int minCommon = max(1, static_cast<int>(std::floor(effExpected - 3.0 * sd)));
        cerr << "effKeysPerSeq=" << effKeys
             << "  effExpected=" << effExpected
             << "  minCommon=" << minCommon << "\n";

        // Inline distance computation with direct-address CSR
        {
            string finalPath = a.output;
            { struct stat st; if (stat(finalPath.c_str(),&st)==0 && S_ISDIR(st.st_mode)){
                if (finalPath.back()!='/') finalPath+='/';
                finalPath+="rabbitsketch.dist";}}

            const int nt = min(a.threads, N);
            const int procId = static_cast<int>(getpid());
            vector<string> parts(nt);
            for (int tid=0;tid<nt;++tid)
                parts[tid] = finalPath+".part."+to_string(procId)+"."+to_string(tid);

            int progress = N/20; if(progress<1) progress=1;

            #pragma omp parallel num_threads(nt)
            {
                const int tid = omp_get_thread_num();
                FILE* tf = fopen(parts[tid].c_str(),"wb");
                if(!tf) { cerr<<"ERROR: cannot open "<<parts[tid]<<"\n"; }
                else {
                    setvbuf(tf,nullptr,_IOFBF,1<<22);
                    string buf; buf.reserve(1<<22);
                    vector<uint16_t> isect(N,0);
                    vector<int>      stamp(N,0);
                    int ep=0;
                    vector<int> cand; cand.reserve(4096);

                    #pragma omp for schedule(dynamic,4)
                    for (int i=0;i<N;++i) {
                        const double si=sizes[i];
                        cand.clear(); ++ep;
                        if(__builtin_expect(ep==INT_MAX,0)){
                            memset(stamp.data(),0,N*sizeof(int));ep=1;}

                        for (uint64_t key : skKeys[i]) {
                            const uint32_t* pl  = dposts.data()+doff[key];
                            const uint32_t  plN = static_cast<uint32_t>(doff[key+1]-doff[key]);
                            for (uint32_t pi=0;pi<plN;++pi){
                                int j=static_cast<int>(pl[pi]);
                                if(j<=i) continue;
                                if(__builtin_expect(stamp[j]!=ep,1)){
                                    stamp[j]=ep;isect[j]=1;cand.push_back(j);
                                } else { ++isect[j]; }
                            }
                        }

                        for (int j:cand){
                            if(isect[j]<static_cast<uint16_t>(minCommon)) continue;
                            const double sj=sizes[j];
                            if(si>0.0 && sj/si<minJac) continue;
                            // Exact verify: Ertl joint MLE Jaccard → Mash distance.
                            // After IDF pruning + minCommon filter, only truly similar
                            // pairs reach here, so ERTL calls are infrequent.
                            const double jac_d = vhll[i]->distance(*vhll[j]); // 1-J
                            const double jac = 1.0 - jac_d;
                            if(jac < minJac) continue;
                            const double d = (jac >= 1.0) ? 0.0
                                : -std::log(2.0*jac/(1.0+jac))
                                  / static_cast<double>(a.kmerSize);
                            if(d < a.maxDist){
                                char line[1024];
                                int len=snprintf(line,sizeof(line),"%s\t%s\t%.6f\n",
                                    fileList[i].c_str(),fileList[j].c_str(),d);
                                buf.append(line,static_cast<size_t>(len));
                            }
                            if(buf.size()>(1<<22)){fwrite(buf.data(),1,buf.size(),tf);buf.clear();}
                        }
                        if(i%progress==0){
                            #pragma omp critical
                            cerr<<"  dist "<<i<<" / "<<N<<"\n";
                        }
                    }
                    if(!buf.empty()) fwrite(buf.data(),1,buf.size(),tf);
                    fclose(tf);
                }
            }
            FILE* fout=fopen(finalPath.c_str(),"wb");
            if(fout){
                setvbuf(fout,nullptr,_IOFBF,1<<24);
                vector<char> cbuf(1<<20);
                for(int tid=0;tid<nt;++tid){
                    FILE* pf=fopen(parts[tid].c_str(),"rb");
                    if(!pf) continue;
                    while(true){size_t got=fread(cbuf.data(),1,cbuf.size(),pf);if(!got)break;fwrite(cbuf.data(),1,got,fout);}
                    fclose(pf); remove(parts[tid].c_str());
                }
                fclose(fout);
            }
            cerr<<"output: "<<finalPath<<"\n";
        }
#endif
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

                const double jac = vhll[i]->jaccard_index(*vhll[j]);
                if (jac < minJac) continue;
                const double dist = (jac >= 1.0) ? 0.0
                    : -std::log(2.0 * jac / (1.0 + jac))
                      / static_cast<double>(a.kmerSize);
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
    if (std::isnan(a.maxDist) || a.maxDist < 0.0) {
        std::cerr << "ERROR: -d must be nonnegative or +infinity\n";
        return 1;
    }

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
        case Algo::FRACMINHASH: std::cerr << "FracMinHash"; break;
        default: break;
    }
    // KSSD list comparisons use the optimized index by default. The explicit
    // diagnostic override preserves the native estimator as a validation oracle.
    const bool useIndex = a.useIndex ||
        (a.algo == Algo::KSSD && !std::getenv("RABBIT_FORCE_ALLPAIRS"));
    std::cerr << (useIndex ? "  (inverted index)\n" : "  (all-pairs)\n");

    if (useIndex) {
        switch (a.algo) {
            case Algo::MINHASH:     run_index_minhash(a, files);    return 0;
            case Algo::FASTKMV:     run_index_fastkmv(a, files);    return 0;
            case Algo::FRACMINHASH: run_index_fracminhash(a, files); return 0;
            case Algo::PROBMINHASH: run_index_probmh(a, files);     return 0;
            case Algo::BINDASH:     return run_index_bindash(a, files);
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
            [cutoff = a.maxDist](Sketch::Kssd* x, Sketch::Kssd* y) {
                return x->distance(y, cutoff);
            });
        break;

    case Algo::PROBMINHASH:
        list_allpairs<Sketch::ProbMinHash4>(a, files, build_probminhash,
            [kmer](Sketch::ProbMinHash4* x, Sketch::ProbMinHash4* y) {
                return mash_distance(x->jaccard_weighted(*y), kmer);
            },
            // Post-build: free ted_params_ (~40 KB) + perm_ (~8 KB) per sketch.
            // Runs after all sketches are built; no malloc-contention during build.
            [](Sketch::ProbMinHash4* sk) { sk->finalize(); });
        break;

    case Algo::SETSKETCH:
        list_allpairs<Sketch::SetSketch>(a, files, build_setsketch,
            [kmer](Sketch::SetSketch* x, Sketch::SetSketch* y) {
                return mash_distance(x->jaccard_index(*y), kmer);
            });
        break;

    case Algo::FASTKMV:
        list_allpairs<Sketch::FastKMV>(a, files, build_fastkmv,
            [cutoff = a.maxDist](Sketch::FastKMV* x, Sketch::FastKMV* y) {
                // The CLI historically reports distance 1 for zero overlap.
                // Broad cutoffs must keep those rows, including -d inf.
                if (cutoff >= 1.0)
                    return mash_distance(x->jaccard(*y), x->getKmerSize());
                return x->distance(*y, cutoff);
            });
        break;

    case Algo::FRACMINHASH:
        list_allpairs<Sketch::FracMinHash>(a, files, build_fracminhash,
            [](Sketch::FracMinHash* x, Sketch::FracMinHash* y) {
                return x->distance(*y);
            });
        break;

    default: return 1;
    }

    return 0;
}
