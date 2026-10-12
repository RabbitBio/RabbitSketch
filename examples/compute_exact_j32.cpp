/*
 * compute_exact_j32.cpp
 *
 * Reads existing FASTA pairs from bench_data/Lvar or Lasym and computes
 * exact Binary Jaccard with k=32 (matching SetSketch/HLL's internal k-mer size).
 *
 * Usage:
 *   ./compute_exact_j32 <data_dir> <mode>
 *     mode = "lvar"  → reads pairs_meta.tsv with col: pair_id, rate, len_bp
 *     mode = "lasym" → reads pairs_meta.tsv with col: pair_id, rate, len_A, len_B
 *
 * Outputs: <data_dir>/pairs_exact_j32.tsv
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>

static const int K = 32;
static const int COMP[4] = {3, 2, 1, 0};

// Read a plain-text FASTA file (no gz) into a string
static std::string read_fasta(const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) return "";
    std::string seq;
    char buf[65536];
    bool in_seq = false;
    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == '>') { in_seq = true; continue; }
        if (!in_seq) continue;
        int n = (int)strlen(buf);
        while (n > 0 && (buf[n-1]=='\n'||buf[n-1]=='\r')) --n;
        seq.append(buf, n);
    }
    fclose(f);
    return seq;
}

static double exact_jaccard_k32(const std::string& a, const std::string& b) {
    auto collect = [&](const std::string& seq) {
        std::vector<uint64_t> v;
        v.reserve(seq.size());
        const uint64_t mask  = (K >= 32) ? ~0ULL : (1ULL<<(2*K))-1;
        const int      hi    = 2*(K-1);
        uint64_t fwd=0, rev=0; int valid=0;
        for (unsigned char ch : seq) {
            int b2;
            switch(ch){
                case 'A':case 'a':b2=0;break;
                case 'C':case 'c':b2=1;break;
                case 'G':case 'g':b2=2;break;
                case 'T':case 't':b2=3;break;
                default:valid=0;fwd=rev=0;continue;
            }
            fwd=((fwd<<2)|b2)&mask;
            rev=(rev>>2)|((uint64_t)COMP[b2]<<hi);
            if(++valid>=K) v.push_back(fwd<=rev?fwd:rev);
        }
        std::sort(v.begin(),v.end());
        v.erase(std::unique(v.begin(),v.end()),v.end());
        return v;
    };
    auto va = collect(a), vb = collect(b);
    size_t i=0,j=0,inter=0;
    while(i<va.size()&&j<vb.size()){
        if(va[i]<vb[j])++i;
        else if(va[i]>vb[j])++j;
        else{++inter;++i;++j;}
    }
    size_t uni=va.size()+vb.size()-inter;
    return uni?(double)inter/uni:0.0;
}

struct PairInfo {
    int    pid;
    double rate;
    int    la, lb;   // lb only used in lasym
};

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <data_dir> <lvar|lasym>\n", argv[0]);
        return 1;
    }
    const std::string dir(argv[1]);
    const bool is_lvar = (strcmp(argv[2], "lvar") == 0);

    // Read pairs_meta.tsv
    std::vector<PairInfo> pairs;
    {
        FILE* f = fopen((dir + "/pairs_meta.tsv").c_str(), "r");
        if (!f) { fprintf(stderr, "Cannot open pairs_meta.tsv\n"); return 1; }
        char line[256]; fgets(line, sizeof(line), f); // skip header
        while (fgets(line, sizeof(line), f)) {
            PairInfo p{};
            if (is_lvar) {
                double theo;
                if (sscanf(line, "%d\t%lf\t%d\t%lf", &p.pid, &p.rate, &p.la, &theo) >= 3) {
                    p.lb = p.la;
                    pairs.push_back(p);
                }
            } else {
                double ratio, theo;
                if (sscanf(line, "%d\t%lf\t%d\t%d\t%lf\t%lf",
                           &p.pid, &p.rate, &p.la, &p.lb, &ratio, &theo) >= 4) {
                    pairs.push_back(p);
                }
            }
        }
        fclose(f);
    }
    fprintf(stderr, "Loaded %zu pairs from %s\n", pairs.size(), dir.c_str());

    // Compute exact Jaccard k=32 in parallel
    std::vector<double> results(pairs.size(), 0.0);

    #pragma omp parallel for schedule(dynamic,1) num_threads(32)
    for (int i = 0; i < (int)pairs.size(); i++) {
        char fa[512], fb[512];
        snprintf(fa, sizeof(fa), "%s/pair%06d_A.fa", dir.c_str(), pairs[i].pid);
        snprintf(fb, sizeof(fb), "%s/pair%06d_B.fa", dir.c_str(), pairs[i].pid);
        std::string seqA = read_fasta(fa);
        std::string seqB = read_fasta(fb);
        if (seqA.empty() || seqB.empty()) {
            fprintf(stderr, "Warning: empty sequence for pair %d\n", pairs[i].pid);
            continue;
        }
        results[i] = exact_jaccard_k32(seqA, seqB);
        if (i % 1000 == 0)
            fprintf(stderr, "  %d/%zu done\n", i, pairs.size());
    }

    // Write output
    FILE* fout = fopen((dir + "/pairs_exact_j32.tsv").c_str(), "w");
    if (is_lvar)
        fprintf(fout, "pair_id\trate\tlen_bp\texact_j_k32\n");
    else
        fprintf(fout, "pair_id\trate\tlen_A\tlen_B\texact_j_k32\n");

    for (int i = 0; i < (int)pairs.size(); i++) {
        if (is_lvar)
            fprintf(fout, "%d\t%.4f\t%d\t%.8f\n",
                    pairs[i].pid, pairs[i].rate, pairs[i].la, results[i]);
        else
            fprintf(fout, "%d\t%.4f\t%d\t%d\t%.8f\n",
                    pairs[i].pid, pairs[i].rate, pairs[i].la, pairs[i].lb, results[i]);
    }
    fclose(fout);
    fprintf(stderr, "Written: %s/pairs_exact_j32.tsv\n", dir.c_str());
    return 0;
}
