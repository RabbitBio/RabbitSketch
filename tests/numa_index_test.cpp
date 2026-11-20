#include "InvertedIndex.h"
#include "BinDashCompressedIndex.h"
#include "api/SketchIndex.h"
#include <map>
#include <set>

using namespace Sketch;
void check(bool value, const char* why) { if (!value) throw std::runtime_error(why); }
constexpr int workers = 16;
NumaIndexOptions localOptions() { NumaIndexOptions o; o.replicate = true; o.max_replica_bytes = 32u << 20; return o; }

template<class Index, class Key>
void queryAll(const Index& index, const std::map<Key, std::vector<uint32_t>>& oracle, size_t nodes) {
    std::atomic<bool> okay{true};
    #pragma omp parallel num_threads(workers)
    {
#ifdef __linux__
        const auto before = Numa::Mask::affinity();
        Numa::Mask oldMemory;
        int oldPolicy = -1;
        const bool havePolicy = syscall(SYS_get_mempolicy, &oldPolicy, oldMemory.words.data(), oldMemory.bits(), nullptr, 0) == 0;
#endif
        {
            auto local = index.local();
            if (nodes > 1 && !local.locally_bound) okay = false;
            for (const auto& item : oracle) {
                size_t off = 0; uint32_t count = 0;
                bool found = local.findPostingRange(item.first, off, count);
                if (item.second.size() < 2) { if (found) okay = false; continue; }
                if (!found) { okay = false; continue; }
                std::vector<uint32_t> ids;
                local.forEachPosting(off, count, [&](uint32_t id) { ids.push_back(id); });
                if (ids != item.second) okay = false;
            }
        }
#ifdef __linux__
        if (Numa::Mask::affinity().words != before.words) okay = false;
        if (havePolicy) {
            Numa::Mask afterMemory; int afterPolicy = -1;
            if (syscall(SYS_get_mempolicy, &afterPolicy, afterMemory.words.data(), afterMemory.bits(), nullptr, 0) ||
                afterPolicy != oldPolicy || afterMemory.words != oldMemory.words) okay = false;
        }
#endif
    }
    check(okay, "posting multiset, local binding, or caller policy restoration failed");
    check(index.bindingFailures() == 0, "query worker fell back to remote replica");
    check(index.replicaCount() == (nodes > 1 ? nodes : 0), "wrong node count");
    check(index.placementMatches(), "replica backing pages are on the wrong node");
}

void generic(size_t nodes) {
    std::vector<std::vector<uint64_t>> keys(2048);
    std::map<uint64_t, std::vector<uint32_t>> oracle;
    for (size_t i = 0; i < keys.size(); ++i) {
        keys[i] = {UINT64_MAX, 0, 7, 7, i * 107 + 25};
        if (i % 3 == 0) keys[i].push_back(17);
        for (auto k : keys[i]) oracle[k].push_back(i);
    }
    auto index = makeNumaIndex(buildCSRIndexFromKeys<uint64_t>(keys, keys.size(), workers), workers, localOptions());
    queryAll(index, oracle, nodes);
    auto shared = makeNumaIndex(buildCSRIndexFromKeys<uint64_t>(keys, keys.size(), workers), workers);
    queryAll(shared, oracle, 0);
    if (nodes > 1) {
        auto limited = localOptions(); limited.max_replica_bytes = 1;
        bool rejected = false;
        try { auto tooBig = makeNumaIndex(buildCSRIndexFromKeys<uint64_t>(keys, keys.size(), workers), workers, limited); }
        catch (const std::length_error&) { rejected = true; }
        check(rejected, "replica memory cap ignored");
    }
    auto empty = makeNumaIndex(InvertedIndex<uint64_t>{}, workers, localOptions());
    check(empty.replicaCount() == 0 && empty.replicaBytes() == 0, "empty index allocated replicas");
}

void direct(size_t nodes) {
    namespace BI = BinDashIndex;
    for (auto mode : {BI::LookupMode::SortedValues, BI::LookupMode::Bitmap}) {
        constexpr size_t n = 512, bins = 96;
        std::vector<std::vector<uint16_t>> values(n, std::vector<uint16_t>(bins));
        std::vector<BI::ValueKeyView<uint64_t, uint16_t>> views(n);
        std::map<uint64_t, std::vector<uint32_t>> oracle;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < bins; ++j) {
                values[i][j] = (i % 13 == 0) ? uint16_t(i) : uint16_t((i + j) % 23);
                oracle[(j << 16) | values[i][j]].push_back(i);
            }
            views[i] = {values[i].data(), bins, 16};
        }
        BI::DirectIndexBuildOptions options; options.lookup_mode = mode;
        auto index = makeNumaIndex(BI::buildCompressedDirectIndex(views, n, workers, options), workers, localOptions());
        queryAll(index, oracle, nodes);
    }
}

void api(size_t nodes) {
    using namespace API;
    SketchConfig config; config.algorithm = Algorithm::FastKMV; config.resolution = 128; config.kmer_size = 21;
    MultiSketchBuilder builder({config});
    IO::FastxRecord record; record.name = "test";
    uint64_t x = 11;
    for (int i = 0; i < 5000; ++i) { x = x * 6364136223846793005ULL + 1; record.sequence += "ACGT"[x >> 62]; }
    builder.update(record);
    auto sketch = builder.finalize("test").front().sketch;
    std::vector<BuiltSketch> sketches(workers, sketch);
    SketchIndex shared(sketches, workers), local(sketches, workers, localOptions());
    check(local.numaReplicaCount() == (nodes > 1 ? nodes : 0), "public API did not enable NUMA");
    auto expected = shared.pairs(0.1, workers), observed = local.pairs(0.1, workers);
    check(expected.size() == observed.size() && observed.size() == 120, "API pairs missing");
    for (size_t i = 0; i < observed.size(); ++i)
        check(observed[i].left == expected[i].left && observed[i].right == expected[i].right &&
              observed[i].estimate.value == expected[i].estimate.value, "API pair result changed");
    auto keys = sketch.indexKeys();
    check(local.postings(keys[0].value, keys[0].position) == shared.postings(keys[0].value, keys[0].position), "API postings changed");
    check(local.candidates(0).size() == shared.candidates(0).size(), "API candidates changed");
    check(local.search(0, 0.1).size() == shared.search(0, 0.1).size(), "API search changed");
}

int main() {
    try {
        const auto topology = Numa::Topology::detect(workers);
        const size_t nodes = topology.nodes.size();
        if (std::getenv("RABBIT_REQUIRE_NUMA")) check(nodes >= 2, "test requires two usable nodes");
        generic(nodes); direct(nodes); api(nodes);
        std::cout << "numa_index_test: PASS nodes=" << nodes << '\n';
    } catch (const std::exception& error) {
        std::cerr << "numa_index_test: FAIL: " << error.what() << '\n'; return 1;
    }
}
