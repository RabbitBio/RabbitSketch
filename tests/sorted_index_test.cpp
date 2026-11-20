#include "SortedInvertedIndex.h"
#include <map>
#include <random>
#include <stdexcept>

namespace {
void check(bool okay, const char* message) {
    if (!okay) throw std::runtime_error(message);
}

template<class Key>
void checkIndex(const std::vector<std::vector<Key>>& keys, int threads,
                size_t workspace = size_t(8) << 30) {
    std::map<Key, std::vector<uint32_t>> oracle;
    std::vector<Sketch::IndexKeyView<Key>> views;
    for (size_t i = 0; i < keys.size(); ++i) {
        check(std::is_sorted(keys[i].begin(), keys[i].end()), "test keys are not sorted");
        views.push_back({keys[i].data(), keys[i].size()});
        for (Key key : keys[i]) oracle[key].push_back(static_cast<uint32_t>(i));
    }
    Sketch::SortedIndexBuildStats stats;
    const auto index = Sketch::buildCSRIndexFromSortedKeys<Key>(
        views, static_cast<int>(keys.size()), threads, {1024, workspace}, &stats);
    size_t expected_keys = 0, expected_posts = 0;
    for (const auto& entry : oracle) {
        const auto found = index.postIdx.find(entry.first);
        if (entry.second.size() == 1) {
            check(found == index.postIdx.end(), "singleton retained");
            continue;
        }
        ++expected_keys;
        expected_posts += entry.second.size();
        check(found != index.postIdx.end(), "key missing or transformed");
        const auto range = found->second;
        check(range.off < index.encodedPosts.size(), "invalid encoded range");
        std::vector<uint32_t> ids;
        index.forEachPosting(range.off, range.cnt, [&](uint32_t sid) { ids.push_back(sid); });
        check(ids == entry.second, "posting IDs differ from independent oracle");
    }
    check(index.postIdx.size() == expected_keys, "unexpected index key");
    check(index.totalPostings == expected_posts,
          "incorrect posting count");
    check(stats.sort_workspace_bound_bytes <= workspace,
          "temporary sorting or merging exceeded the workspace budget");
}

template<class Key>
void cases() {
    for (int threads : {1, 4, 16}) {
        checkIndex<Key>({}, threads);
        checkIndex<Key>({{}, {}, {}}, threads);
        checkIndex<Key>({{0, 1, 2}}, threads);
        checkIndex<Key>({{0, 0, 1}, {0, 0, 0, 1}, {2, 2}, {}}, threads);
        checkIndex<Key>({{0, 1, std::numeric_limits<Key>::max()},
                         {0, 2, std::numeric_limits<Key>::max()}, {}, {3}}, threads);
    }
    std::mt19937_64 random(20261120);
    std::vector<std::vector<Key>> keys(11);
    for (uint64_t i = 0; i < 240000; ++i) {
        const Key key = static_cast<Key>(random());
        keys[i % 7].push_back(key);
        if (i % 3) keys[(i + 2) % 7].push_back(key);
    }
    keys[7] = keys[0];
    keys[8] = {0};
    keys[9] = {0}; // One empty sketch and a very skewed length distribution.
    for (auto& values : keys) {
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());
    }
    for (int threads : {1, 4, 16}) checkIndex(keys, threads, 4 << 20);
    // Dense low keys plus repeated cut values exercise bottom-k ranges and
    // ensure equal keys are never split between partitions.
    keys.assign(9, {});
    for (auto& values : keys)
        for (uint64_t i = 0; i < 80000; ++i) values.push_back(static_cast<Key>(i));
    checkIndex(keys, 16, 1 << 20);
    // Skewed partitions exceed the budget but must still produce all repeated
    // occurrences, including zero deltas, without an oversized pair buffer.
    keys.assign(11, std::vector<Key>(15000, Key(7)));
    checkIndex(keys, 16, 4096); // bounded k-way merge
    checkIndex(keys, 16, 1);    // no heap fits: constant-scratch rescan
    keys[3].push_back(Key(9)); keys[8].push_back(Key(9));
    checkIndex(keys, 16, 4096);
}
} // namespace

int main() {
    try {
        cases<uint32_t>();
        cases<uint64_t>();
        std::cout << "sorted_index_test: PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "sorted_index_test: FAIL: " << e.what() << '\n';
        return 1;
    }
}
