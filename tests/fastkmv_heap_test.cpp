#include "fastkmv.h"
#include "rank/CanonicalKmer.h"
#include "rank/RankStream.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void check(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

std::string dna(size_t n, uint64_t seed) {
    std::mt19937_64 random(seed);
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[random() & 3];
    return s;
}

std::vector<uint64_t> keys(const std::vector<std::string>& records, int kmer, uint64_t seed) {
    std::vector<uint64_t> result;
    for (const auto& s : records) {
        Sketch::Rank::CanonicalKmerIterator it(s.data(), s.size(), kmer);
        uint64_t code;
        while (it.next(code)) {
            uint64_t h = code;
            if (Sketch::FastKMV::runtimeHashProfile() != Sketch::Rank::HashProfile::CanonicalBits) {
                h = Sketch::Rank::RankStream::fmix64(h, seed);
                if (Sketch::FastKMV::runtimeHashProfile() == Sketch::Rank::HashProfile::Murmur3Fmix64Twice)
                    h = Sketch::Rank::RankStream::fmix64(h, seed);
            }
            result.push_back(h >> 11);
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

void matches(const Sketch::FastKMV& s, const std::vector<uint64_t>& expected) {
    const size_t n = std::min(size_t(s.getK()), expected.size());
    check(s.size() == n, "unique retained count differs from oracle");
    // Cardinality must work before requesting any sorted view.
    const double cardinality = n < s.getK() || n < 2 || expected[n-1] == 0
        ? double(n) : double(n-1)*double(Sketch::FastKMV::KEY_MAX)/double(expected[n-1]);
    check(std::abs(s.cardinality()-cardinality) <= 1e-12*std::max(1.0,cardinality),
          "heap threshold cardinality differs from oracle");
    const uint64_t* registers = s.getRegisters();
    check(std::equal(registers, registers+n, expected.begin()), "bottom-k differs from full sort oracle");
    check(s.getRegisters() == registers, "repeated ordered read changed storage");
    for (size_t i=n; i<s.getK(); ++i)
        check(registers[i] == UINT64_MAX, "partial sketch lost padding sentinel");
}

void streamingAndQueries() {
    for (int kmer : {1, 21, 32}) for (uint32_t k : {2u, 8u, 64u, 1024u, 4096u}) {
        Sketch::FastKMV s(k, kmer, 42);
        std::vector<std::string> records;
        matches(s, keys(records, kmer, 42));
        std::string ambiguous = dna(12000, 17);
        ambiguous.replace(80, 6, "NnSDwU");
        ambiguous.replace(900, 4, "acgt");
        const std::vector<std::string> additions = {
            "", "ACG", std::string(90, 'A'), dna(150, 11),
            ambiguous, ambiguous, dna(18000, 23), std::string(600, 'A')};
        for (const auto& added : additions) {
            s.update(added.data(), added.size());
            records.push_back(added);
            matches(s, keys(records, kmer, 42)); // Query followed by more writes.
        }
        s.finalize(); s.finalize();
        matches(s, keys(records, kmer, 42));
        bool rejected = false;
        try { s.update("ACGT", 4); } catch (const std::logic_error&) { rejected = true; }
        check(rejected, "sealed sketch accepted an update");
        auto copy = s;
        matches(copy, keys(records, kmer, 42));
    }
}

void copiesRestoreProjectMerge() {
    const auto a = dna(20000, 71), b = dna(15000, 91), c = dna(12000, 113);
    Sketch::FastKMV source(256, 21, 42);
    source.update(a.data(), a.size()); // Copy while still in heap order.
    Sketch::FastKMV copy(source), assigned(8, 21, 42);
    assigned = source;
    copy.update(b.data(), b.size());
    assigned.update(b.data(), b.size());
    matches(copy, keys({a,b},21,42));
    matches(assigned, keys({a,b},21,42));
    matches(source, keys({a},21,42));

    Sketch::FastKMV sorted_copy(source); // Copy after an ordered query.
    Sketch::FastKMV moved(std::move(sorted_copy));
    moved.update(b.data(), b.size());
    matches(moved, keys({a,b},21,42));

    std::vector<uint64_t> saved(source.getRegisters(),source.getRegisters()+source.size());
    auto restored = Sketch::FastKMV::fromRegisters(256,21,42,saved);
    restored.update(a.data(),a.size()); // All existing retained keys are duplicates.
    matches(restored,keys({a},21,42));
    restored.update(b.data(),b.size());
    matches(restored,keys({a,b},21,42));

    auto projected = source.project(64);
    projected.update(b.data(),b.size());
    matches(projected,keys({a,b},21,42));
    Sketch::FastKMV other(64,21,42);
    other.update(b.data(),b.size());
    auto merged = source.merge(other);
    matches(merged,keys({a,b},21,42));
    merged.update(c.data(),c.size());
    matches(merged,keys({a,b,c},21,42));

    auto partial = Sketch::FastKMV::fromRegisters(256,21,42,{});
    partial.update("AAAAAAAAAAAAAAAAAAAAAAAAA",25);
    matches(partial,keys({std::string(25,'A')},21,42));
}

double referenceJaccard(std::vector<uint64_t> a,std::vector<uint64_t> b,size_t k) {
    if(a.empty() || b.empty()) return 0;
    if(a.size()>k) a.resize(k);
    if(b.size()>k) b.resize(k);
    std::vector<uint64_t> joined;
    std::set_union(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(joined));
    if(joined.size()>k) joined.resize(k);
    size_t common=0;
    for(uint64_t h:joined)
        common+=std::binary_search(a.begin(),a.end(),h) && std::binary_search(b.begin(),b.end(),h);
    return double(common)/joined.size();
}

void jaccardAfterLazyReads() {
    const auto shared = dna(10000,109), left = dna(12000,131), right = dna(14000,157);
    for(uint32_t k:{8u,64u,1024u}) {
        Sketch::FastKMV a(k,21,42), b(k*2,21,42);
        a.update(shared.data(),shared.size()); b.update(shared.data(),shared.size());
        check(a.jaccard(b)==1.0,"initial shared sketches do not match");
        a.update(left.data(),left.size()); b.update(right.data(),right.size());
        const double expected=referenceJaccard(keys({shared,left},21,42),keys({shared,right},21,42),k);
        check(std::abs(a.jaccard(b)-expected)<1e-12,"Jaccard after resumed build differs from oracle");
        check(std::abs(a.jaccard(b)-expected)<1e-12,"repeated Jaccard changed result");
        a.finalize(); b.finalize();
        check(std::abs(a.jaccard(b)-expected)<1e-12,"finalization changed Jaccard");
        check(std::abs(a.jaccard(b,expected*.9)-expected)<1e-12,"threshold pruned a qualifying pair");
    }
}
}

int main() {
    try {
        streamingAndQueries();
        copiesRestoreProjectMerge();
        jaccardAfterLazyReads();
        std::cout << "fastkmv_heap_test: PASS\n";
    } catch(const std::exception& error) {
        std::cerr << "fastkmv_heap_test: FAIL: " << error.what() << '\n';
        return 1;
    }
}
