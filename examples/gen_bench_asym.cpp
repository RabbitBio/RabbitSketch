/*
 * gen_bench_asym.cpp  –  Generate asymmetric-length simulated sequence pairs
 *
 * A is derived from the first L_A bases of the reference (mutated),
 * B is derived from the full L_B bases of the reference (mutated).
 * This models the scenario where one sequence is much shorter than the other
 * (e.g. a short query vs a long reference).
 *
 * Exact Jaccard is computed by canonical k-mer set intersection/union.
 *
 * Usage:  ./gen_bench_asym <out_dir> [seed]
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <omp.h>

static const int    K       = 21;
static const int    N_PAIRS = 500;
static const int    WRAP    = 80;

// Asymmetric length pairs: (L_A, L_B) – A is shorter, B is longer
struct LenPair { int la; int lb; };
static const std::vector<LenPair> LEN_PAIRS = {
    {100000,   1000000},   // 100k vs 1M  (ratio 10×)
    {100000,   8000000},   // 100k vs 8M  (ratio 80×)
    {500000,   8000000},   // 500k vs 8M  (ratio 16×)
    {1000000,  8000000},   // 1M   vs 8M  (ratio 8×)
    {500000,   4000000},   // 500k vs 4M  (ratio 8×)
    {200000,   2000000},   // 200k vs 2M  (ratio 10×)
};
static const std::vector<double> RATES = {0.001, 0.010, 0.050, 0.100, 0.200};

static const char BASES[4] = {'A','C','G','T'};
static const int  COMP[4]  = { 3, 2, 1, 0 };

static void write_fasta(const char* path, const char* name, const std::string& seq) {
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); exit(1); }
    fprintf(f, ">%s\n", name);
    int n = (int)seq.size();
    for (int i = 0; i < n; i += WRAP) {
        int end = std::min(i+WRAP, n);
        fwrite(seq.data()+i, 1, end-i, f);
        fputc('\n', f);
    }
    fclose(f);
}

static std::string mutate(const std::string& seq, double rate, std::mt19937_64& rng) {
    std::bernoulli_distribution     mut_d(rate);
    std::uniform_int_distribution<int> alt_d(0, 2);
    std::string out = seq;
    for (char& c : out) {
        if (mut_d(rng)) {
            int orig;
            switch(c){case 'A':orig=0;break;case 'C':orig=1;break;
                      case 'G':orig=2;break;default:orig=3;break;}
            c = BASES[(orig+1+alt_d(rng))%4];
        }
    }
    return out;
}

static double exact_jaccard(const std::string& a, const std::string& b) {
    auto collect = [&](const std::string& seq) {
        std::vector<uint64_t> v;
        v.reserve(seq.size());
        const uint64_t mask  = (1ULL<<(2*K))-1;
        const int hi         = 2*(K-1);
        uint64_t fwd=0, rev=0; int valid=0;
        for (unsigned char ch : seq) {
            int b2;
            switch(ch){case 'A':case 'a':b2=0;break;
                       case 'C':case 'c':b2=1;break;
                       case 'G':case 'g':b2=2;break;
                       case 'T':case 't':b2=3;break;
                       default:valid=0;fwd=rev=0;continue;}
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

int main(int argc, char** argv) {
    if(argc<2){fprintf(stderr,"Usage: %s <out_dir> [seed]\n",argv[0]);return 1;}
    const std::string out_dir(argv[1]);
    uint64_t seed = argc>=3?(uint64_t)atoll(argv[2]):99999ULL;
    mkdir(out_dir.c_str(), 0755);

    FILE* fmeta  = fopen((out_dir+"/pairs_meta.tsv").c_str(),"w");
    FILE* fexact = fopen((out_dir+"/pairs_exact_j.tsv").c_str(),"w");
    FILE* falist = fopen((out_dir+"/A_list.txt").c_str(),"w");
    FILE* fblist = fopen((out_dir+"/B_list.txt").c_str(),"w");
    fprintf(fmeta, "pair_id\trate\tlen_A\tlen_B\tratio\ttheo_j_k21\n");
    fprintf(fexact,"pair_id\trate\tlen_A\tlen_B\texact_j_k21\n");

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> base_rand(0,3);

    int pair_id = 0;
    int total   = (int)(LEN_PAIRS.size()*RATES.size());
    int combo   = 0;

    for (auto& lp : LEN_PAIRS) {
        int LA = lp.la, LB = lp.lb;

        // Generate one reference of length LB
        std::string ref(LB,'A');
        for (char& c : ref) c = BASES[base_rand(rng)];
        const std::string ref_A = ref.substr(0, LA);  // first LA bases

        for (double rate : RATES) {
            ++combo;
            // Theoretical Jaccard for asymmetric case:
            // Expected shared k-mers from overlap ≈ (1-r)^k * LA
            // Expected unique k-mers in B beyond overlap ≈ (1-r)^k * (LB - LA)
            // E[J] = LA_kmers / (LA_kmers + (LB-LA)_kmers) × (1-r)^k contributions cancel
            // Simplified: E[J] ≈ (1-r)^k * LA / ((1-r)^k * LB)  (for large sequences)
            double pk = std::pow(1.0-rate, K);
            double theo_j = pk * LA / (pk * (LA + LB) - pk*pk*LA);
            // More precisely: E[inter]=pk*LA, E[A]=pk*LA, E[B]=pk*LB → E[J]=LA/(LA+LB)·(pk/(2-pk))
            // Use simpler bound: theo ≈ (LA/LB) * (pk/(2-pk))
            theo_j = (double)LA / LB * pk / (2.0 - pk);

            fprintf(stderr,"[%d/%d] LA=%d LB=%d rate=%.3f\n",
                    combo, total, LA, LB, rate);

            std::vector<uint64_t> seeds_A(N_PAIRS), seeds_B(N_PAIRS);
            for(int p=0;p<N_PAIRS;p++){seeds_A[p]=rng();seeds_B[p]=rng();}

            struct PR{std::string seqA,seqB;double ej;};
            std::vector<PR> res(N_PAIRS);

            #pragma omp parallel for schedule(dynamic,1) num_threads(32)
            for(int p=0;p<N_PAIRS;p++){
                std::mt19937_64 ra(seeds_A[p]),rb(seeds_B[p]);
                res[p].seqA = mutate(ref_A, rate, ra);
                res[p].seqB = mutate(ref,   rate, rb);
                res[p].ej   = exact_jaccard(res[p].seqA, res[p].seqB);
            }

            for(int p=0;p<N_PAIRS;p++,pair_id++){
                char fa[512],fb[512],na[80],nb[80];
                snprintf(fa,sizeof(fa),"%s/pair%06d_A.fa",out_dir.c_str(),pair_id);
                snprintf(fb,sizeof(fb),"%s/pair%06d_B.fa",out_dir.c_str(),pair_id);
                snprintf(na,sizeof(na),"pair%06d_A_L%d_r%.3f",pair_id,LA,rate);
                snprintf(nb,sizeof(nb),"pair%06d_B_L%d_r%.3f",pair_id,LB,rate);
                write_fasta(fa,na,res[p].seqA);
                write_fasta(fb,nb,res[p].seqB);
                fprintf(fmeta, "%d\t%.4f\t%d\t%d\t%.2f\t%.6f\n",
                        pair_id,rate,LA,LB,(double)LB/LA,theo_j);
                fprintf(fexact,"%d\t%.4f\t%d\t%d\t%.8f\n",
                        pair_id,rate,LA,LB,res[p].ej);
                fprintf(falist,"%s\n",fa);
                fprintf(fblist,"%s\n",fb);
            }
        }
    }

    fclose(fmeta);fclose(fexact);fclose(falist);fclose(fblist);
    fprintf(stderr,"Done. %d pairs → %s\n",pair_id,out_dir.c_str());
    return 0;
}
