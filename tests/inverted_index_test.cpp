#include "InvertedIndex.h"
#include "probmh.h"
#include <map>
#include <random>
#include <stdexcept>

namespace {
void check(bool okay, const char* message) {
    if (!okay) throw std::runtime_error(message);
}

// Independent reference: group original keys with std::map, then compare
// every retained posting and pair count. No bucket/layout assumptions.
template<class Key>
std::vector<std::vector<uint32_t>> checkIndex(
    const std::vector<std::vector<Key>>& keys, int threads) {
    const int n = static_cast<int>(keys.size());
    std::map<Key, std::vector<uint32_t>> reference;
    std::vector<std::vector<Key>> sorted = keys;
    for (int i = 0; i < n; ++i) {
        std::sort(sorted[i].begin(), sorted[i].end());
        check(std::adjacent_find(sorted[i].begin(), sorted[i].end()) == sorted[i].end(),
              "test input contains duplicate keys");
        for (Key key : keys[i]) reference[key].push_back(i);
    }
    const auto index = Sketch::buildCSRIndexFromKeys<Key>(keys, n, threads);
    size_t expectedKeys = 0, expectedPosts = 0;
    std::vector<bool> visited(index.csrPosts.size());
    for (const auto& item : reference) {
        const auto found = index.postIdx.find(item.first);
        if (item.second.size() == 1) {
            check(found == index.postIdx.end(), "singleton retained");
            continue;
        }
        ++expectedKeys;
        expectedPosts += item.second.size();
        check(found != index.postIdx.end(), "original key missing or transformed");
        const auto range = found->second;
        check(range.off <= index.csrPosts.size() &&
              range.cnt <= index.csrPosts.size() - range.off, "invalid CSR range");
        std::vector<uint32_t> observed(index.csrPosts.begin() + range.off,
            index.csrPosts.begin() + range.off + range.cnt);
        std::sort(observed.begin(), observed.end());
        check(observed == item.second, "posting IDs differ from reference");
        for (size_t p = range.off; p < range.off + range.cnt; ++p) {
            check(!visited[p], "overlapping CSR ranges");
            visited[p] = true;
        }
    }
    check(index.postIdx.size() == expectedKeys, "unexpected index key");
    check(index.csrPosts.size() == expectedPosts && index.totalPostings == expectedPosts,
          "incorrect posting total");
    std::vector<std::vector<uint32_t>> counts(n, std::vector<uint32_t>(n));
    for (int i = 0; i < n; ++i) {
        for (Key key : keys[i]) {
            const auto found = index.postIdx.find(key);
            if (found == index.postIdx.end()) continue;
            const auto range = found->second;
            for (size_t p = range.off; p < range.off + range.cnt; ++p)
                ++counts[i][index.csrPosts[p]];
        }
        for (int j = i + 1; j < n; ++j) {
            std::vector<Key> common;
            std::set_intersection(sorted[i].begin(), sorted[i].end(),
                sorted[j].begin(), sorted[j].end(), std::back_inserter(common));
            check(counts[i][j] == common.size(), "pair count differs from set intersection");
        }
    }
    return counts;
}

template<class Key>
void integerCases() {
    for (int threads : {1, 4, 16}) {
        checkIndex<Key>({}, threads);
        checkIndex<Key>({{}, {}, {}}, threads);
        checkIndex<Key>({{0, 1, 2}}, threads);
        checkIndex<Key>({{0, 1, std::numeric_limits<Key>::max()},
                         {0, 2, std::numeric_limits<Key>::max()}, {}, {3}}, threads);
    }
    std::mt19937_64 random(20261112);
    for (bool lowRange : {true, false}) {
        std::vector<std::vector<Key>> keys(9);
        for (uint64_t i = 0; i < 12000; ++i) {
            // Bottom-k/threshold-style low keys, or unrestricted integer keys.
            const Key key = lowRange ? static_cast<Key>(i * 17) : static_cast<Key>(random());
            keys[i % 7].push_back(key);
            if (i % 3) keys[(i + 2) % 7].push_back(key);
        }
        keys[7] = keys[0]; // Identical sketches; last sketch remains empty.
        for (auto& values : keys) {
            std::sort(values.begin(), values.end());
            values.erase(std::unique(values.begin(), values.end()), values.end());
            std::shuffle(values.begin(), values.end(), random);
        }
        for (int threads : {1, 4, 16}) checkIndex(keys, threads);
    }
}

void winnerCases() {
    constexpr uint32_t registers = 64;
    std::vector<Sketch::ProbMinHash4> sketches;
    std::vector<std::vector<uint64_t>> keys(5);
    for (int i = 0; i < 5; ++i) {
        sketches.emplace_back(registers, 21, 42);
        for (uint64_t item = 1; item <= 200; ++item) {
            if (i == 4 || (i == 2 && item % 3 == 0)) continue;
            const uint64_t hash = item + (i == 3 ? 1000 : 0);
            sketches.back().addHash(hash, 1.0 + (item % 11));
        }
        sketches.back().getWinnerIndexKeys(keys[i]);
    }
    for (int threads : {1, 4, 16}) {
        const auto counts = checkIndex(keys, threads);
        for (int i = 0; i < 4; ++i) for (int j = i + 1; j < 4; ++j)
            check(Sketch::ProbMinHash4::jaccardFromCommon(counts[i][j], registers) ==
                  sketches[i].jaccard_weighted(sketches[j]),
                  "winner posting counts changed the native estimator");
    }
}
} // namespace

int main() {
    try {
        integerCases<uint32_t>();
        integerCases<uint64_t>();
        winnerCases();
        std::cout << "inverted_index_test: PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "inverted_index_test: FAIL: " << e.what() << '\n';
        return 1;
    }
}
