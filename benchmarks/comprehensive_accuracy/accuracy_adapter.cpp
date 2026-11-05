// Accuracy-only adapter; every estimate is produced by a native backend.
#define RABBITSKETCH_SCREENING_EMBEDDED
#include "../screening/synthetic_pilot.cpp"
#include "sub_command.h"
#include <omp.h>
#include <fstream>

int run_regimes_manifest(const std::string&, int);

namespace comprehensive {
struct Input {
    std::string id;
    int k, ratio, replicate;
    uint64_t na, nb, common, data_seed, hash_seed;
};

template<class Make, class Update, class Finish, class Bytes, class Empty,
         class Jaccard, class Cardinality, class Containment>
void evaluate(const Input& x, Dataset& data, const std::string& algorithm,
              int parameter, int budget, int half_subk, uint64_t da, uint64_t db,
              Make make, Update update, Finish finish, Bytes bytes, Empty empty,
              Jaccard jaccard, Cardinality cardinality, Containment containment) {
    auto a = make(), b = make();
    update(*a, data.a); update(*b, data.b); finish(*a); finish(*b);
    const auto ba = bytes(*a), bb = bytes(*b);
    const bool ea = empty(*a), eb = empty(*b);
    const double j = jaccard(*a, *b), ca = cardinality(*a), cb = cardinality(*b);
    double cab = NAN, cba = NAN;
    if (std::isfinite(j) && ca > 0 && cb > 0 && std::isfinite(ca) && std::isfinite(cb)) {
        double intersection = j * (ca + cb) / (1 + j);
        intersection = std::max(0., std::min(intersection, std::min(ca, cb)));
        cab = intersection / ca; cba = intersection / cb;
    }
    Json row;
    row.put("unit_id", x.id); row.put("algorithm", algorithm);
    row.integer("parameter", parameter); row.integer("requested_budget", budget);
    row.integer("half_subk", half_subk); row.integer("k", x.k);
    row.integer("replicate", x.replicate); row.integer("data_seed", x.data_seed);
    row.integer("hash_seed", algorithm == "kssd" ? 42 : x.hash_seed);
    row.integer("n_a", x.na); row.integer("n_b", x.nb); row.integer("intersection", x.common);
    row.integer("input_digest_a", da); row.integer("input_digest_b", db);
    row.integer("payload_bytes_a", ba); row.integer("payload_bytes_b", bb);
    row.integer("empty_a", ea); row.integer("empty_b", eb);
    row.put("jaccard", j); row.put("cardinality_a", ca); row.put("cardinality_b", cb);
    row.put("containment_ab", cab); row.put("containment_ba", cba);
    row.put("native_containment_ab", containment(*a, *b));
    row.put("status", std::isfinite(j) && j >= 0 && j <= 1 ? "ok" : "invalid_estimate");
    row.print();
}

void one(const Input& x, bool verify) {
    auto data = generate(x.na, x.nb, x.common, x.k, x.data_seed);
    if (verify) {
        const auto a = exactScreeningWords(data.a, x.k), b = exactScreeningWords(data.b, x.k);
        if (a.size() != x.na || b.size() != x.nb || exactScreeningIntersection(a, b) != x.common)
            throw std::runtime_error("independent exact truth verification failed");
    }
    const auto da = digest(data.a), db = digest(data.b);
    auto up = [](auto& s, std::string& data) { s.update(data.data(), data.size()); };
    auto upstr = [](auto& s, std::string& data) { s.update(data.data()); };
    auto finish = [](auto& s) { s.finalize(); };
    auto noop = [](auto&) {};
    auto empty = [](auto& s) { return s.size() == 0; };
    auto notempty = [](auto&) { return false; };
    auto hashbytes = [](auto& s) { return s.size() * 8; };
    auto registerbytes = [](auto& s) { return s.getCore().size(); };
    auto card = [](auto& s) { return s.cardinality(); };
    auto jaccard = [](auto& a, auto& b) { return a.jaccard(b); };
    auto rjaccard = [](auto& a, auto& b) { return a.jaccard_index(b); };
    auto contain = [](auto& a, auto& b) { return a.containment(b); };
    auto nocard = [](auto&) { return NAN; };
    auto nocontain = [](auto&, auto&) { return NAN; };
    for (int budget : {512, 2048, 8192, 32768}) {
        const int m = budget / 16, bins = budget / 4;
        const int p = static_cast<int>(std::log2(budget / 2));
        const int scaled = std::max(1, static_cast<int>(std::ceil(8. * (x.na + x.nb) / (.8 * budget))));
        evaluate(x, data, "minhash", m, budget, 6, da, db,
            [&] { return std::make_unique<Sketch::FastKMV>(m, x.k, x.hash_seed); },
            up, finish, hashbytes, empty, jaccard, card, contain);
        evaluate(x, data, "fracminhash", scaled, budget, 6, da, db,
            [&] { return std::make_unique<Sketch::FracMinHash>(scaled, x.k, x.hash_seed); },
            up, finish, hashbytes, empty, jaccard, card, contain);
        evaluate(x, data, "hll", p, budget, 6, da, db,
            [&] { return std::make_unique<Sketch::HyperLogLog>(p, x.k, x.hash_seed, false); },
            upstr, noop, registerbytes, notempty, rjaccard, card, contain);
        evaluate(x, data, "setsketch", p, budget, 6, da, db,
            [&] { return std::make_unique<Sketch::SetSketch>(p, 2., 5., x.k, false, x.hash_seed); },
            upstr, noop, registerbytes, notempty, rjaccard, card, contain);
        evaluate(x, data, "probminhash", m, budget, 6, da, db,
            [&] { return std::make_unique<Sketch::ProbMinHash4>(m, x.k, x.hash_seed, 0); },
            up, finish, [](auto& s) { return s.getM() * 8; }, notempty, jaccard, nocard, nocontain);
        evaluate(x, data, "bindash", bins, budget, 6, da, db,
            [&] { return std::make_unique<Sketch::BinDash>(bins / 64, x.k, 16, x.hash_seed); },
            up, finish, [](auto& s) { return s.memoryBytes(); },
            [](auto& s) { return s.getRawNonempty() == 0; }, jaccard, card, contain);
    }
    for (int dr : {0, 1, 2, 3, 4}) {
        const int half = dr == 4 ? 7 : 6;
        const auto params = sharedKssdParameters(x.k / 2, half, dr);
        evaluate(x, data, "kssd", dr, 0, half, da, db,
            [&] { return std::make_unique<Sketch::Kssd>(params, 42); }, upstr, finish,
            [](auto& s) { return s.getHashes().size() * 4 + s.getHashes64().size() * 8; },
            [](auto& s) { return s.getHashes().empty() && s.getHashes64().empty(); },
            [](auto& a, auto& b) { return a.jaccard(&b); },
            [](auto& s) { return (s.getHashes().size() + s.getHashes64().size()) / s.admissionProbability(); },
            nocontain);
    }
}
}

int run_screening_manifest(const Args& args) {
    try {
        std::ifstream in(args.screeningManifest);
        if (!in) throw std::runtime_error("missing accuracy manifest");
        std::string line; std::getline(in, line);
        if (line == "#regimes-v1") return run_regimes_manifest(args.screeningManifest, args.threads);
        const bool verify = line == "#comprehensive-verify-v1";
        if (!verify && line != "#comprehensive-v1") throw std::runtime_error("bad accuracy marker");
        std::vector<comprehensive::Input> inputs;
        std::unordered_set<std::string> ids;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream stream(line); comprehensive::Input x;
            if (!(stream >> x.id >> x.k >> x.na >> x.nb >> x.common >> x.replicate >> x.data_seed >> x.hash_seed))
                throw std::runtime_error("bad accuracy input");
            if (!ids.insert(x.id).second || !x.na || !x.nb || x.common > std::min(x.na, x.nb) ||
                x.k < 20 || x.k > 32 || x.k % 2 || args.threads < 1)
                throw std::runtime_error("invalid accuracy input");
            if (!inputs.empty() && x.k != inputs[0].k) throw std::runtime_error("mixed k in batch");
            inputs.push_back(x);
        }
        if (inputs.empty()) throw std::runtime_error("empty accuracy manifest");
        for (int dr : {0, 1, 2, 3, 4}) screeningWarmKssd(inputs[0].k, dr == 4 ? 7 : 6, dr);
        omp_set_dynamic(0);
        int failed = 0;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(args.threads) reduction(+:failed)
        for (size_t i = 0; i < inputs.size(); ++i) {
            try { comprehensive::one(inputs[i], verify); }
            catch (const std::exception& e) {
                Json row; row.put("unit_id", inputs[i].id); row.put("status", "error");
                row.put("reason", e.what()); row.print(); ++failed;
            }
        }
        return failed ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << "comprehensive accuracy: " << e.what() << '\n'; return 2;
    }
}
