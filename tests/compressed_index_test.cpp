#include "InvertedIndex.h"
#include "api/SketchBuilder.h"
#include <map>
#include <random>
#include <stdexcept>

void check(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

template<class Key>
void compare(const Sketch::InvertedIndex<Key>& index,
             const std::vector<std::vector<Key>>& input) {
    std::map<Key, std::vector<uint32_t>> oracle;
    for (size_t i = 0; i < input.size(); ++i)
        for (const auto& key : input[i]) oracle[key].push_back(i);
    size_t unique = 0, posts = 0;
    for (const auto& item : oracle) {
        size_t off = 0; uint32_t count = 0;
        bool found = index.findPostingRange(item.first, off, count);
        if (item.second.size() == 1) { check(!found, "singleton retained"); continue; }
        check(found, "key missing");
        std::vector<uint32_t> observed;
        index.forEachPosting(off, count, [&](uint32_t sid) { observed.push_back(sid); });
        check(observed == item.second, "posting multiset or order changed");
        ++unique; posts += observed.size();
    }
    check(index.postIdx.size() == unique && index.totalPostings == posts, "index counts changed");
}

void codec() {
    const std::vector<uint32_t> deltas = {0, 1, 127, 128, 16383, 16384, 2097151,
                                        2097152, 268435455, 268435456, UINT32_MAX};
    const std::vector<std::vector<uint8_t>> bytes = {{0}, {1}, {127}, {128,1},
        {255,127}, {128,128,1}, {255,255,127}, {128,128,128,1},
        {255,255,255,127}, {128,128,128,128,1}, {255,255,255,255,15}};
    for (size_t i = 0; i < deltas.size(); ++i) {
        uint8_t actual[5];
        const auto length = Sketch::PostingCodec::encode(actual, deltas[i]);
        check(length == bytes[i].size() && std::equal(bytes[i].begin(), bytes[i].end(), actual),
              "noncanonical varint or lost high bits");
        const uint8_t* p = bytes[i].data();
        check(Sketch::PostingCodec::decode(p) == deltas[i] && p == bytes[i].data() + bytes[i].size(),
              "varint boundary decode failed");
    }
    Sketch::InvertedIndex<uint64_t> index;
    index.encodedPosts = {0, 0, 255, 255, 255, 255, 15, 0};
    std::vector<uint32_t> ids;
    index.forEachPosting(0, 4, [&](uint32_t id) { ids.push_back(id); });
    check(ids == std::vector<uint32_t>({0,0,UINT32_MAX,UINT32_MAX}), "sample IDs or duplicates truncated");
}

template<class Key>
void integers() {
    std::vector<std::vector<Key>> input(10000);
    for (size_t i = 0; i < input.size(); ++i) {
        input[i] = {Key(7), Key(0), Key(7), std::numeric_limits<Key>::max()};
        if (i % 7 == 0) input[i].push_back(Key(17));
    }
    for (size_t budget : {size_t(1), size_t(4096), size_t(1 << 20)}) {
        const auto original = input;
        Sketch::IndexBuildStats stats;
        auto idx = Sketch::buildCSRIndexFromKeys<Key>(input, input.size(), 16,
            Sketch::MixedIndexBucket{}, {1024, budget}, &stats);
        compare(idx, original);
        check(stats.sort_workspace_bound_bytes <= budget, "workspace cap exceeded");
        check(idx.encodedPosts.size() == idx.totalPostings, "dense IDs did not compress to one byte");
        // Const unsorted entry point is bounded too, and cannot mutate input.
        idx = Sketch::buildCSRIndexFromKeys<Key>(original, original.size(), 16,
            Sketch::MixedIndexBucket{}, {1024, budget}, &stats);
        compare(idx, original);
        check(stats.sort_workspace_bound_bytes <= budget, "const workspace cap exceeded");
        input = original;
    }
    std::vector<phmap::flat_hash_map<Key,std::vector<uint32_t>>> maps(3);
    maps[0][Key(7)] = {9, 0, 9}; maps[1][Key(7)] = {1, 9};
    maps[2][Key(11)] = {5};
    auto index = Sketch::buildCSRIndex(maps, 4);
    std::vector<uint32_t> ids;
    size_t off = 0; uint32_t count = 0;
    check(index.findPostingRange(Key(7), off, count), "legacy merged key lost");
    index.forEachPosting(off, count, [&](uint32_t sid) { ids.push_back(sid); });
    check(ids == std::vector<uint32_t>({0,1,9,9,9}), "legacy duplicates lost");
    check(index.totalPostings == 5 && index.postIdx.size() == 1, "legacy singleton or count changed");
}

int main() {
    try {
        codec(); integers<uint32_t>(); integers<uint64_t>();
        using Key = Sketch::API::IndexKey;
        std::vector<std::vector<Key>> keys = {
            {{UINT64_MAX, UINT32_MAX}, {0,0}, {UINT64_MAX,0}, {0,0}},
            {{0,1}, {UINT64_MAX, UINT32_MAX}, {0,0}}, {}, {{0,1}}};
        auto original = keys;
        const auto index = Sketch::buildCSRIndexFromKeys<Key>(keys, keys.size(), 16);
        compare(index, original);
        std::cout << "compressed_index_test: PASS\n";
    } catch (const std::exception& error) {
        std::cerr << "compressed_index_test: FAIL: " << error.what() << '\n';
        return 1;
    }
}
