// Smoke test: SetSketch inclusion-exclusion vs joint-MLE estimator
//
// Constructs two random DNA strings with a controlled overlap pattern,
// then compares
//   J_true  : exact k-mer Jaccard (from set intersection)
//   J_incl  : SetSketch inclusion-exclusion estimate  (jaccard_index_inclexcl)
//   J_mle   : SetSketch joint-MLE estimate            (jaccard_index_mle, the
//                                                       new default of
//                                                       jaccard_index())
// across a grid of (true J, bits, a).

#include "SetSketch.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <cmath>

static const char BASES[4] = {'A','C','G','T'};

// Pair generation: build a single ancestral random sequence, then create
// A by substituting each base with prob p_a, B similarly. Result is two
// sequences with controlled divergence d; expected J(k) ≈ (1-d)^k for k-mers.
struct Pair { std::string a, b; };
static Pair gen_pair(size_t L, double d, std::mt19937_64& rng) {
    std::uniform_int_distribution<int> base4(0, 3);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::string ancestor(L, 'A');
    for (size_t i = 0; i < L; ++i) ancestor[i] = BASES[base4(rng)];
    auto mut = [&](const std::string& s) {
        std::string out = s;
        for (size_t i = 0; i < L; ++i) {
            if (u(rng) < d) {
                int x = base4(rng);
                while (BASES[x] == out[i]) x = base4(rng);
                out[i] = BASES[x];
            }
        }
        return out;
    };
    return { mut(ancestor), mut(ancestor) };
}

// Canonical k-mer hash by encoding 2 bits per base (A=0,C=1,G=2,T=3).
static uint64_t encode(const std::string& s, size_t pos, int k) {
    uint64_t x = 0;
    for (int i = 0; i < k; ++i) {
        int c = -1;
        switch (s[pos + i]) {
            case 'A': c = 0; break;
            case 'C': c = 1; break;
            case 'G': c = 2; break;
            case 'T': c = 3; break;
            default: return UINT64_MAX;
        }
        x = (x << 2) | (uint64_t)c;
    }
    return x;
}
static uint64_t rc(uint64_t x, int k) {
    uint64_t y = 0;
    for (int i = 0; i < k; ++i) {
        y = (y << 2) | ((~x) & 3);
        x >>= 2;
    }
    return y;
}
static double exact_kmer_jaccard(const std::string& A, const std::string& B, int k) {
    std::unordered_set<uint64_t> SA, SB;
    auto fill = [&](const std::string& s, std::unordered_set<uint64_t>& S) {
        if ((int)s.size() < k) return;
        for (size_t i = 0; i + (size_t)k <= s.size(); ++i) {
            uint64_t h = encode(s, i, k);
            if (h == UINT64_MAX) continue;
            uint64_t hr = rc(h, k);
            S.insert(h < hr ? h : hr);
        }
    };
    fill(A, SA); fill(B, SB);
    size_t inter = 0;
    if (SA.size() > SB.size()) std::swap(SA, SB);
    for (auto x : SA) if (SB.count(x)) inter++;
    size_t un = SA.size() + SB.size() - inter;
    return un ? (double)inter / (double)un : 0.0;
}

int main() {
    const int k = 21;
    const size_t L = 200000;          // 200 kb
    const int n_trials = 20;
    std::mt19937_64 rng(42);

    // (bits, a)
    std::vector<std::pair<int,double>> configs = {
        {13,  5.0},    // RabbitSketch default
        {13, 3.0},     // alisim-recommended a
        {13, 10.0},
        {11,  20.0},
        {10,  20.0},
        {10,   3.0},
        {15,  20.0},
        {16,  20.0},
    };

    // divergence d → expected k-mer Jaccard ≈ (1-d)^k
    std::vector<double> divs = { 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15 };

    printf("%6s  %6s  %5s  %6s  %8s  %8s  %8s  %9s  %9s\n",
           "bits", "a", "trial", "div", "J_true", "J_incl", "J_mle", "err_incl", "err_mle");

    for (auto [b, av] : configs) {
        for (double d : divs) {
            double sum_true = 0, sum_incl = 0, sum_mle = 0;
            double err_incl = 0, err_mle = 0;
            for (int t = 0; t < n_trials; ++t) {
                auto pr = gen_pair(L, d, rng);
                Sketch::SetSketch s1(b, 2.0, av, k), s2(b, 2.0, av, k);
                s1.update(const_cast<char*>(pr.a.c_str()));
                s2.update(const_cast<char*>(pr.b.c_str()));
                double j_true = exact_kmer_jaccard(pr.a, pr.b, k);
                double j_incl = s1.jaccard_index_inclexcl(s2);
                double j_mle  = s1.jaccard_index_mle(s2);
                sum_true += j_true; sum_incl += j_incl; sum_mle += j_mle;
                err_incl += std::fabs(j_incl - j_true);
                err_mle  += std::fabs(j_mle  - j_true);
            }
            double Jt = sum_true / n_trials;
            double Ji = sum_incl / n_trials;
            double Jm = sum_mle  / n_trials;
            err_incl /= n_trials; err_mle /= n_trials;
            printf("%6d  %6.2f  %5d  %6.3f  %8.5f  %8.5f  %8.5f  %9.5f  %9.5f\n",
                   b, av, n_trials, d, Jt, Ji, Jm, err_incl, err_mle);
        }
        printf("\n");
    }
    return 0;
}
