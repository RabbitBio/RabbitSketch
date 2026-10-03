#pragma once
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

// A separator at each record boundary prevents artificial inter-contig k-mers.
inline std::string readScreeningFasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read FASTA: " + path);
    std::string sequence, line;
    bool record = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            sequence.push_back('N'); record = true; continue;
        }
        if (!record && !line.empty()) throw std::runtime_error("FASTA header missing");
        for (unsigned char c : line) {
            if (std::isspace(c)) continue;
            c = std::toupper(c);
            sequence.push_back(c=='A'||c=='C'||c=='G'||c=='T' ? c : 'N');
        }
    }
    if (!record || sequence.empty()) throw std::runtime_error("empty FASTA: " + path);
    return sequence;
}

// Exact 2-bit canonical words, not probabilistic hashes. k <= 32.
inline std::vector<uint64_t> exactScreeningWords(const std::string& s, int k) {
    if (k < 1 || k > 32) throw std::invalid_argument("invalid exact k");
    const uint64_t mask = k==32 ? UINT64_MAX : (UINT64_C(1)<<(2*k))-1;
    uint64_t forward=0, reverse=0;
    int valid=0;
    std::vector<uint64_t> words;
    words.reserve(s.size());
    for (char c : s) {
        const int v=c=='A'?0:c=='C'?1:c=='G'?2:c=='T'?3:-1;
        if (v<0) { forward=reverse=0; valid=0; continue; }
        forward=((forward<<2)|v)&mask;
        reverse=(reverse>>2)|(uint64_t(3-v)<<(2*(k-1)));
        if (++valid>=k) words.push_back(std::min(forward,reverse));
    }
    std::sort(words.begin(),words.end());
    words.erase(std::unique(words.begin(),words.end()),words.end());
    return words;
}

inline uint64_t exactScreeningIntersection(const std::vector<uint64_t>& a,
                                           const std::vector<uint64_t>& b) {
    size_t i=0,j=0; uint64_t common=0;
    while (i<a.size() && j<b.size()) {
        if (a[i]<b[j]) ++i;
        else if (a[i]>b[j]) ++j;
        else { ++common; ++i; ++j; }
    }
    return common;
}
