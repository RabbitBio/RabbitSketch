// Reuse the actual CLI FASTX/sketch constructors, including PMH frequency
// normalization. Only the benchmark driver and key snapshot format are new.
#include "../../examples/sub_command.cpp"
#include <chrono>
#include <iomanip>
#include <sys/resource.h>

namespace {
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point begin) {
    return std::chrono::duration<double>(Clock::now() - begin).count();
}
long peakRSS() {
    struct rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage)) throw std::runtime_error("getrusage failed");
    return usage.ru_maxrss;
}
template<class T> void writeValue(std::ostream& out, T value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}
template<class T> T readValue(std::istream& in) {
    T value{};
    in.read(reinterpret_cast<char*>(&value), sizeof(value));
    if (!in) throw std::runtime_error("truncated snapshot");
    return value;
}
void writeSnapshot(const std::string& path, const Args& args,
                   const std::vector<std::vector<uint64_t>>& keys) {
    const std::string temporary = path + ".tmp";
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("cannot create snapshot");
    writeValue<uint64_t>(out, UINT64_C(0x5253494e44580001));
    writeValue<uint32_t>(out, args.algo == Algo::FASTKMV ? 1 : 2);
    writeValue<uint32_t>(out, args.kmerSize);
    writeValue<uint32_t>(out, args.fkmvK);
    writeValue<uint64_t>(out, args.seed);
    writeValue<uint64_t>(out, keys.size());
    for (const auto& values : keys) {
        writeValue<uint32_t>(out, values.size());
        out.write(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(uint64_t));
    }
    out.close();
    if (!out || std::rename(temporary.c_str(), path.c_str()))
        throw std::runtime_error("cannot finish snapshot");
}
std::vector<std::vector<uint64_t>> readSnapshot(const std::string& path,
                                               const Args& args, size_t genomes) {
    std::ifstream in(path, std::ios::binary);
    if (!in || readValue<uint64_t>(in) != UINT64_C(0x5253494e44580001))
        throw std::runtime_error("invalid snapshot magic/version");
    if (readValue<uint32_t>(in) != (args.algo == Algo::FASTKMV ? 1u : 2u) ||
        readValue<uint32_t>(in) != static_cast<uint32_t>(args.kmerSize) ||
        readValue<uint32_t>(in) != args.fkmvK || readValue<uint64_t>(in) != args.seed ||
        readValue<uint64_t>(in) != genomes)
        throw std::runtime_error("snapshot configuration differs");
    std::vector<std::vector<uint64_t>> keys(genomes);
    for (auto& values : keys) {
        const auto count = readValue<uint32_t>(in);
        if (count > args.fkmvK) throw std::runtime_error("invalid snapshot key count");
        values.resize(count);
        in.read(reinterpret_cast<char*>(values.data()), count * sizeof(uint64_t));
        if (!in) throw std::runtime_error("truncated snapshot keys");
    }
    if (in.peek() != std::char_traits<char>::eof()) throw std::runtime_error("trailing snapshot data");
    return keys;
}

int prepare(const Args& args, const std::vector<std::string>& files,
            const std::string& snapshot) {
    const auto begin = Clock::now();
    const int n = files.size();
    std::vector<std::vector<uint64_t>> keys(n);
    std::atomic<int> completed{0};
    #pragma omp parallel for num_threads(args.threads) schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        if (args.algo == Algo::FASTKMV) {
            std::unique_ptr<Sketch::FastKMV> sk(build_fastkmv(files[i], args));
            const uint64_t* values = sk->getRegisters();
            keys[i].assign(values, values + sk->size());
        } else {
            std::unique_ptr<Sketch::ProbMinHash4> sk(build_probminhash(files[i], args));
            sk->getWinnerIndexKeys(keys[i]);
        }
        const int done = ++completed;
        if (done % 4096 == 0 || done == n) {
            #pragma omp critical
            std::cerr << "prepared " << done << " / " << n << '\n';
        }
    }
    const double buildSeconds = seconds(begin);
    const auto saveBegin = Clock::now();
    writeSnapshot(snapshot, args, keys);
    const double saveSeconds = seconds(saveBegin);
    uint64_t total = 0, empty = 0;
    for (const auto& values : keys) { total += values.size(); empty += values.empty(); }
    std::cout << std::setprecision(12) << "{\"files\":" << n << ",\"keys\":" << total
              << ",\"empty_sketches\":" << empty << ",\"sketch_seconds\":" << buildSeconds
              << ",\"snapshot_write_seconds\":" << saveSeconds
              << ",\"peak_rss_kib\":" << peakRSS() << "}\n";
    return 0;
}

int replay(const Args& args, const std::vector<std::string>& files,
           const std::string& snapshot, const std::string& mode) {
    const auto begin = Clock::now();
    auto keys = readSnapshot(snapshot, args, files.size());
    const double loadSeconds = seconds(begin);
    const int n = files.size();
    std::vector<std::unique_ptr<Sketch::FastKMV>> native;
    const auto nativeBegin = Clock::now();
    if (args.algo == Algo::FASTKMV) {
        native.resize(n);
        #pragma omp parallel for num_threads(args.threads) schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            native[i] = std::make_unique<Sketch::FastKMV>(Sketch::FastKMV::fromRegisters(
                args.fkmvK, args.kmerSize, args.seed, keys[i]));
            native[i]->finalize();
        }
    }
    const double restoreSeconds = seconds(nativeBegin);
    std::cerr << "snapshot loaded; starting " << mode << " index\n";
    const auto indexBegin = Clock::now();
    Sketch::InvertedIndex<uint64_t> index;
    if (mode == "map") {
        std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>> local(args.threads);
        #pragma omp parallel num_threads(args.threads)
        {
            auto& table = local[omp_get_thread_num()];
            #pragma omp for schedule(dynamic)
            for (int i = 0; i < n; ++i)
                for (uint64_t key : keys[i]) table[key].push_back(i);
        }
        index = Sketch::buildCSRIndex(local, args.threads);
    } else {
        index = Sketch::buildCSRIndexFromKeys(keys, n, args.threads);
    }
    const double indexSeconds = seconds(indexBegin);
    const long indexRSS = peakRSS();
    std::cerr << "complete index time: " << indexSeconds << " s\n";
    const double p = std::exp(-args.kmerSize * args.maxDist);
    const double minJac = p / (2.0 - p);
    const auto queryBegin = Clock::now();
    if (args.algo == Algo::FASTKMV) {
        auto minimum = [minJac](int size) { return std::max(1, int(std::ceil(minJac * size))); };
        auto jaccard = [&](int i, int j) { return native[i]->jaccard(*native[j], minJac); };
        Sketch::computeDistancesExact(index, keys, files, n, args.kmerSize,
            args.maxDist, jaccard, minimum, args.output, args.threads);
    } else {
        const int required = std::max(1, int(std::ceil(minJac * args.pmhM)));
        auto minimum = [required](int) { return required; };
        auto jaccard = [&](int common, int, int) {
            return Sketch::ProbMinHash4::jaccardFromCommon(common, args.pmhM);
        };
        const std::vector<int> sizes(n, args.pmhM);
        Sketch::computeDistances(index, keys, sizes, files, n, args.kmerSize,
            args.maxDist, jaccard, minimum, args.output, args.threads);
    }
    const double querySeconds = seconds(queryBegin);
    std::cout << std::setprecision(12) << "{\"files\":" << n << ",\"load_seconds\":" << loadSeconds
              << ",\"restore_seconds\":" << restoreSeconds << ",\"index_seconds\":" << indexSeconds
              << ",\"index_peak_rss_kib\":" << indexRSS << ",\"index_keys\":" << index.postIdx.size()
              << ",\"postings\":" << index.csrPosts.size() << ",\"query_seconds\":" << querySeconds
              << ",\"index_query_seconds\":" << indexSeconds + querySeconds
              << ",\"internal_total_seconds\":" << seconds(begin) << "}\n";
    return 0;
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 7 && argc != 8)
            throw std::runtime_error("usage: index_build_233 prepare|map|bucket fastkmv|probminhash files.list snapshot threads output [K]");
        const std::string mode = argv[1], algorithm = argv[2];
        if (mode != "prepare" && mode != "map" && mode != "bucket") throw std::runtime_error("unknown mode");
        Args args;
        if (algorithm == "fastkmv") args.algo = Algo::FASTKMV;
        else if (algorithm == "probminhash") args.algo = Algo::PROBMINHASH;
        else throw std::runtime_error("unknown algorithm");
        args.threads = std::stoi(argv[5]);
        if (args.threads < 1) throw std::runtime_error("invalid thread count");
        args.output = argv[6];
        if (argc == 8) args.fkmvK = args.pmhM = std::stoul(argv[7]);
        const auto files = load_file_list(argv[3]);
        if (files.empty() || files.size() > INT_MAX) throw std::runtime_error("invalid file count");
        return mode == "prepare" ? prepare(args, files, argv[4]) : replay(args, files, argv[4], mode);
    } catch (const std::exception& error) {
        std::cerr << "index_build_233: " << error.what() << '\n';
        return 1;
    }
}
