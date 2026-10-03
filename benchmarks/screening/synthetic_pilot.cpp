// One isolated native-backend measurement. See PROTOCOL.md for measurement scope.
#include "Sketch.h"
#include "fastkmv.h"
#include "fracminhash.h"
#include "BinDash.h"
#include "probmh.h"
#include "api/RuntimeInfo.h"
#include "real_input.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <sys/resource.h>

using Clock = std::chrono::steady_clock;
static double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
static uint64_t peakRss() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage)) throw std::runtime_error("getrusage failed");
    return static_cast<uint64_t>(usage.ru_maxrss) * 1024; // Linux KiB
}
static std::string quote(const std::string& value) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u" << std::hex << std::setw(4)
                             << std::setfill('0') << static_cast<unsigned>(c);
        else out << c;
    }
    return out.str() + '"';
}
struct Json {
    std::map<std::string, std::string> fields;
    void put(const std::string& key, const std::string& value) { fields[key] = quote(value); }
    void put(const std::string& key, double value) {
        std::ostringstream out;
        out << std::setprecision(17) << value;
        fields[key] = std::isfinite(value) ? out.str() : "null";
    }
    void integer(const std::string& key, uint64_t value) { fields[key] = std::to_string(value); }
    void print() const {
        std::ostringstream out;
        out << '{';
        bool first = true;
        for (const auto& field : fields) {
            if (!first) out << ',';
            first = false;
            out << quote(field.first) << ':' << field.second;
        }
        out << "}\n";
        static std::mutex output_mutex;
        const std::lock_guard<std::mutex> lock(output_mutex);
        std::cout << out.str() << std::flush;
    }
};

// Match main/sub_command's shared KSSD dictionary policy. Construction calls
// libc srand/rand, so it must happen once, before the OpenMP sample loop.
static std::shared_ptr<const Sketch::kssd_parameter_t> sharedKssdParameters(
    int half_k, int half_subk, int drlevel) {
    static std::mutex mutex;
    static std::map<std::tuple<int, int, int>,
                    std::shared_ptr<const Sketch::kssd_parameter_t>> cache;
    const std::lock_guard<std::mutex> lock(mutex);
    const auto key = std::make_tuple(half_k, half_subk, drlevel);
    auto found = cache.find(key);
    if (found == cache.end())
        found = cache.emplace(key, std::make_shared<Sketch::kssd_parameter_t>(
            half_k, half_subk, drlevel)).first;
    return found->second;
}

void screeningWarmKssd(int k, int half_subk, int drlevel) {
    (void)sharedKssdParameters(k / 2, half_subk, drlevel);
}

static uint64_t reverseComplement(uint64_t word, int k) {
    uint64_t result = 0;
    for (int i = 0; i < k; ++i) {
        result = (result << 2) | (3 - (word & 3));
        word >>= 2;
    }
    return result;
}
static void appendWord(std::string& result, uint64_t word, int k) {
    for (int i = k - 1; i >= 0; --i) result.push_back("ACGT"[(word >> (2 * i)) & 3]);
    result.push_back('N'); // No k-mer can cross a set-element boundary.
}
struct Dataset { std::string a, b; };
static Dataset generate(uint64_t na, uint64_t nb, uint64_t common, int k, uint64_t seed) {
    const uint64_t total = na + nb - common;
    const uint64_t mask = k == 32 ? UINT64_MAX : (UINT64_C(1) << (2 * k)) - 1;
    // Stay far below saturation, including reverse-complement equivalence.
    if (total > mask / 4) throw std::invalid_argument("set too large for chosen k");
    std::mt19937_64 rng(seed);
    std::unordered_set<uint64_t> seen;
    seen.reserve(total);
    std::vector<uint64_t> a, b;
    a.reserve(na); b.reserve(nb);
    for (uint64_t index = 0; index < total;) {
        uint64_t word = rng() & mask;
        word = std::min(word, reverseComplement(word, k));
        if (!seen.insert(word).second) continue;
        if (index < common) { a.push_back(word); b.push_back(word); }
        else if (index < na) a.push_back(word);
        else b.push_back(word);
        ++index;
    }
    std::shuffle(a.begin(), a.end(), rng);
    std::shuffle(b.begin(), b.end(), rng);
    Dataset result;
    result.a.reserve(na * (k + 1)); result.b.reserve(nb * (k + 1));
    for (uint64_t word : a) appendWord(result.a, word, k);
    for (uint64_t word : b) appendWord(result.b, word, k);
    return result;
}
static uint64_t digest(const std::string& value) {
    uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned char c : value) { hash ^= c; hash *= UINT64_C(1099511628211); }
    return hash;
}
// Independent enumeration of carrier strings, used by --verify-data checks.
static std::unordered_set<uint64_t> enumerate(const std::string& s, int k) {
    std::unordered_set<uint64_t> result;
    for (size_t i = 0; i + k <= s.size(); ++i) {
        uint64_t word = 0;
        bool valid = true;
        for (int j = 0; j < k; ++j) {
            const auto base = std::string("ACGT").find(s[i + j]);
            if (base == std::string::npos) { valid = false; break; }
            word = (word << 2) | base;
        }
        if (valid) result.insert(std::min(word, reverseComplement(word, k)));
    }
    return result;
}

struct MeasuredPair {
    std::function<double()> query;
    std::function<double()> cardinality_a, cardinality_b;
    uint64_t bytes_a = 0, bytes_b = 0;
    bool empty_a = false, empty_b = false;
    double setup_seconds = 0, build_a_seconds = 0, build_b_seconds = 0;
};

template<class T, class Make, class Update, class Finalize, class Bytes, class Empty, class Query>
static MeasuredPair buildPair(Dataset& data, Make make, Update update, Finalize finalize,
                      Bytes bytes, Empty empty, Query query) {
    MeasuredPair result;
    auto start = Clock::now();
    auto a = make(); update(*a, data.a); finalize(*a);
    result.build_a_seconds = seconds(start);
    start = Clock::now();
    auto b = make(); update(*b, data.b); finalize(*b);
    result.build_b_seconds = seconds(start);
    result.bytes_a = bytes(*a); result.bytes_b = bytes(*b);
    result.empty_a = empty(*a); result.empty_b = empty(*b);
    result.query = [a, b, query]() { return query(*a, *b); };
    return result;
}

int screeningMeasureOne(int argc, char** argv) {
    Json out;
    try {
        std::map<std::string, std::string> options;
        for (int i = 1; i < argc; i += 2) {
            if (i + 1 == argc) throw std::invalid_argument("options require values");
            if (!options.emplace(argv[i], argv[i + 1]).second)
                throw std::invalid_argument("duplicate option");
        }
        const std::unordered_set<std::string> allowed = {
            "--algorithm", "--parameter", "--k", "--n-a", "--n-b", "--intersection",
            "--data-seed", "--hash-seed", "--query-min-ms", "--half-subk", "--bbits", "--verify-data",
            "--job-id", "--batch-threads", "--worker-thread", "--fasta-a", "--fasta-b"};
        if (options.count("--job-id")) out.put("job_id", options.at("--job-id"));
        for (const auto& kv : options)
            if (!allowed.count(kv.first)) throw std::invalid_argument("unknown option: " + kv.first);
        auto number = [&](const char* key, uint64_t fallback) {
            const auto it = options.find(key);
            if (it == options.end()) return fallback;
            if (it->second.empty() || it->second[0] == '-') throw std::invalid_argument(key);
            size_t end = 0;
            const uint64_t value = std::stoull(it->second, &end);
            if (end != it->second.size()) throw std::invalid_argument(key);
            return value;
        };
        const std::string algorithm = options.at("--algorithm");
        const auto raw_k = number("--k", 20), raw_p = number("--parameter", 0);
        if (raw_k < 1 || raw_k > 32 || raw_p > UINT32_MAX)
            throw std::invalid_argument("invalid k or parameter");
        const int k = static_cast<int>(raw_k);
        const uint32_t p = static_cast<uint32_t>(raw_p);
        const uint64_t na = number("--n-a", 10000), nb = number("--n-b", 10000);
        const uint64_t common = number("--intersection", 0);
        const uint64_t data_seed = number("--data-seed", 42), seed = number("--hash-seed", 42);
        const uint64_t min_ms = number("--query-min-ms", 5);
        const uint64_t half_subk = number("--half-subk", 6), bbits = number("--bbits", 16);
        const uint64_t batch_threads = number("--batch-threads", 1);
        const bool shared_setup = options.count("--job-id") != 0;
        const bool real_input = options.count("--fasta-a") != 0;
        if (real_input != bool(options.count("--fasta-b")))
            throw std::invalid_argument("both FASTA paths required");
        out.integer("batch_threads", batch_threads);
        out.integer("worker_thread", number("--worker-thread", 0));
        out.put("measurement_profile", real_input ? "main_openmp_real_pair" :
                shared_setup ? "main_openmp_batch" : "serial_pair");
        out.put("rss_scope", shared_setup ? "shared_batch_process" : "isolated_pair_process");
        out.put("setup_scope", shared_setup ? "batch_shared" : "pair_private");
        if (!na || !nb || common > std::min(na, nb) || na > UINT64_MAX - nb ||
            na + nb > SIZE_MAX / static_cast<size_t>(k + 1) || min_ms < 1 || min_ms > 1000)
            throw std::invalid_argument("invalid cardinalities, intersection, or timing budget");
        if ((algorithm == "hll" && (p < 4 || p > 20)) ||
            (algorithm == "setsketch" && (p < 4 || p > 16)) ||
            (algorithm == "bindash" && (!p || p % 64 || bbits < 1 || bbits > 16)) ||
            (algorithm == "probminhash" && p < 2) ||
            ((algorithm == "minhash" || algorithm == "fracminhash") && !p) ||
            (algorithm == "kssd" && (k % 2 || half_subk < 3 || half_subk > 7 ||
                                     half_subk > raw_k / 2 || p + 3 > half_subk)))
            throw std::invalid_argument("illegal native configuration");
        auto gen_start = Clock::now();
        auto data = real_input ? Dataset{readScreeningFasta(options.at("--fasta-a")),
                                         readScreeningFasta(options.at("--fasta-b"))}
                               : generate(na, nb, common, k, data_seed);
        out.put("generation_seconds", seconds(gen_start));
        out.integer("input_digest_a", digest(data.a)); out.integer("input_digest_b", digest(data.b));
        if (algorithm == "exact") {
            if (!real_input) throw std::invalid_argument("exact mode requires FASTA");
            auto a=exactScreeningWords(data.a,k), b=exactScreeningWords(data.b,k);
            const auto intersection=exactScreeningIntersection(a,b);
            if (a.empty() || b.empty()) throw std::invalid_argument("empty exact k-mer set");
            out.integer("n_a",a.size()); out.integer("n_b",b.size());
            out.integer("intersection",intersection);
            out.put("ground_truth",double(intersection)/(a.size()+b.size()-intersection));
            out.integer("input_bytes",data.a.size()+data.b.size());
            out.put("status","ok"); out.print(); return 0;
        }
        if (number("--verify-data", 0)) {
            auto a = enumerate(data.a, k), b = enumerate(data.b, k);
            uint64_t intersection = 0;
            for (uint64_t value : a) intersection += b.count(value);
            if (a.size() != na || b.size() != nb || intersection != common)
                throw std::runtime_error("independent truth verification failed");
            out.integer("verified_distinct_a", a.size()); out.integer("verified_distinct_b", b.size());
            out.integer("verified_intersection", intersection);
        }
        out.integer("process_peak_rss_before_build_bytes", peakRss());
        out.integer("input_bytes", data.a.size() + data.b.size());
        out.put("ground_truth", double(common) / double(na + nb - common));
        out.put("ground_truth_containment_a", double(common) / na);
        out.put("ground_truth_containment_b", double(common) / nb);
        MeasuredPair pair;
        auto no_finalize = [](auto&) {};
        auto not_empty = [](auto&) { return false; };
        auto update_length = [](auto& sketch, std::string& s) { sketch.update(s.data(), s.size()); };
        auto update_string = [](auto& sketch, std::string& s) { sketch.update(s.data()); };
        auto finalize = [](auto& sketch) { sketch.finalize(); };
        auto native_j = [](auto& a, auto& b) { return a.jaccard(b); };
        if (algorithm == "minhash") {
            // User-facing MinHash is FastKMV. Legacy MinHash is intentionally absent.
            pair = buildPair<Sketch::FastKMV>(data,
                [&]() { return std::make_shared<Sketch::FastKMV>(p, k, seed); },
                update_length, finalize, [](auto& s) { return s.size() * sizeof(uint64_t); },
                [](auto& s) { return s.size() == 0; }, native_j);
            out.put("backend_class", "Sketch::FastKMV");
        } else if (algorithm == "fracminhash") {
            pair = buildPair<Sketch::FracMinHash>(data,
                [&]() { return std::make_shared<Sketch::FracMinHash>(p, k, seed); },
                update_length, finalize, [](auto& s) { return s.size() * sizeof(uint64_t); },
                [](auto& s) { return s.size() == 0; }, native_j);
            out.put("backend_class", "Sketch::FracMinHash");
        } else if (algorithm == "hll" || algorithm == "setsketch") {
            // Capture cardinality access separately; querying may fill native caches.
            auto register_pair = [&](auto make) {
                using T = typename decltype(make())::element_type;
                std::vector<std::shared_ptr<T>> sketches;
                auto tracked_make = [&]() { auto s = make(); sketches.push_back(s); return s; };
                auto r = buildPair<T>(data, tracked_make, update_string, no_finalize,
                    [](auto& s) { return s.getCore().size() * sizeof(uint8_t); }, not_empty,
                    [](auto& a, auto& b) { return a.jaccard_index(b); });
                r.cardinality_a = [a = sketches[0]]() { return a->cardinality(); };
                r.cardinality_b = [b = sketches[1]]() { return b->cardinality(); };
                return r;
            };
            if (algorithm == "hll") {
                pair = register_pair([&]() { return std::make_shared<Sketch::HyperLogLog>(p, k, seed, false); });
                out.put("backend_class", "Sketch::HyperLogLog");
            } else {
                pair = register_pair([&]() { return std::make_shared<Sketch::SetSketch>(p, 2.0, 5.0, k, false, seed); });
                out.put("backend_class", "Sketch::SetSketch");
            }
        } else if (algorithm == "probminhash") {
            pair = buildPair<Sketch::ProbMinHash4>(data,
                [&]() { return std::make_shared<Sketch::ProbMinHash4>(p, k, seed, 0); },
                update_length, finalize, [](auto& s) { return s.getM() * sizeof(double); },
                not_empty, native_j);
            out.put("backend_class", "Sketch::ProbMinHash4_unweighted");
        } else if (algorithm == "bindash") {
            pair = buildPair<Sketch::BinDash>(data,
                [&]() { return std::make_shared<Sketch::BinDash>(p / 64, k, bbits, seed); },
                update_length, finalize, [](auto& s) { return s.memoryBytes(); },
                [](auto& s) { return s.getRawNonempty() == 0; }, native_j);
            out.put("backend_class", "Sketch::BinDash");
        } else if (algorithm == "kssd") {
            const auto start = Clock::now();
            auto params = shared_setup ? sharedKssdParameters(k / 2, half_subk, p)
                : std::make_shared<Sketch::kssd_parameter_t>(k / 2, half_subk, p);
            const double setup_seconds = seconds(start);
            pair = buildPair<Sketch::Kssd>(data,
                [&]() { return std::make_shared<Sketch::Kssd>(params, seed); },
                update_string, finalize,
                [](auto& s) { return s.getHashes().size() * sizeof(uint32_t) +
                                    s.getHashes64().size() * sizeof(uint64_t); },
                [](auto& s) { return s.getHashes().empty() && s.getHashes64().empty(); },
                [](auto& a, auto& b) { return a.jaccard(&b); });
            pair.setup_seconds = setup_seconds;
            out.put("backend_class", "Sketch::Kssd");
        } else throw std::invalid_argument("unsupported algorithm: " + algorithm);
        out.integer("process_peak_rss_after_build_bytes", peakRss());
        const auto cold_start = Clock::now();
        const double estimate = pair.query();
        out.put("query_cold_seconds", seconds(cold_start));
        // Separate first call (including native cache fill) from repeated warm comparisons.
        uint64_t queries = 0;
        volatile double sink = 0;
        auto warm_start = Clock::now();
        do {
            for (unsigned i = 0; i < 32; ++i) sink = pair.query();
            queries += 32;
        } while (seconds(warm_start) < double(min_ms) / 1000);
        const double elapsed = seconds(warm_start);
        (void)sink;
        out.put("estimate", estimate);
        out.put("status", std::isfinite(estimate) && estimate >= 0 && estimate <= 1 ? "ok" : "invalid_estimate");
        out.integer("empty_sample_a", pair.empty_a); out.integer("empty_sample_b", pair.empty_b);
        out.integer("payload_bytes_a", pair.bytes_a); out.integer("payload_bytes_b", pair.bytes_b);
        out.integer("sketch_payload_bytes", pair.bytes_a + pair.bytes_b);
        out.put("setup_seconds", pair.setup_seconds);
        out.put("build_a_seconds", pair.build_a_seconds); out.put("build_b_seconds", pair.build_b_seconds);
        out.put("build_seconds", pair.setup_seconds + pair.build_a_seconds + pair.build_b_seconds);
        out.put("query_seconds", elapsed / queries); out.integer("query_iterations", queries);
        out.put("cardinality_estimate_a", pair.cardinality_a ? pair.cardinality_a() : NAN);
        out.put("cardinality_estimate_b", pair.cardinality_b ? pair.cardinality_b() : NAN);
        out.integer("process_peak_rss_bytes", peakRss());
        const auto& runtime = Sketch::Runtime::runtimeInfo();
        out.put("compiler", runtime.compiler); out.put("byte_path", runtime.selected_byte_path);
        out.put("rank_path", runtime.selected_rank_path); out.put("u64_path", runtime.selected_u64_path);
        out.print();
        return 0;
    } catch (const std::exception& e) {
        out.put("status", "error"); out.put("reason", e.what()); out.print();
        return 2;
    }
}

#ifndef RABBITSKETCH_SCREENING_EMBEDDED
int main(int argc, char** argv) { return screeningMeasureOne(argc, argv); }
#endif
