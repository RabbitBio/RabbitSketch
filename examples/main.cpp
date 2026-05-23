/**
 * rabbitsketch – unified CLI front-end for the RabbitSketch library.
 *
 * This file is intentionally thin.  It only:
 *
 *   1. Defines the usage text.
 *   2. Parses command-line arguments into an Args struct.
 *   3. Hands the parsed arguments to run_cli() in sub_command.cpp,
 *      which contains every algorithm pipeline.
 *
 * Algorithms (mutually exclusive):
 *   --minhash      classical bottom-k MinHash
 *   --kssd         dimensionality-reduced k-mer sketch
 *   --hll          HyperLogLog cardinality / Jaccard sketch
 *   --probminhash  frequency-weighted ProbMinHash4 (PMH-norm, J_P estimator)
 *   --bindash      b-bit one-permutation hashing
 *   --setsketch    register-based set sketch
 *   --fastkmv      FastKMV bottom-k sketch
 *
 * Modes:
 *   pair mode      -i fileA fileB                 (two FASTA/FASTQ inputs)
 *   list mode      -i list.txt -l                 (file list, all-pairs)
 *   --index        in list mode, run the inverted-index pipeline
 *                  (currently supported for: --fastkmv --probminhash
 *                                            --kssd     --setsketch)
 *
 * Common options:
 *   -o <path>      output file (default: rabbitsketch.dist)
 *   -k <int>       k-mer size (default per algorithm)
 *   -t <int>       OpenMP threads     (default: omp_get_max_threads())
 *   -d <float>     max Mash distance for filter (default 0.05; use 1.0 for no filter)
 *   --seed <u64>   hash seed          (default 42)
 *
 * Algorithm-specific options:
 *   --minhash       -s   <int>     sketch size                    (default 1000)
 *   --kssd          --drlevel <int>                                (default 3)
 *   --hll           -b   <int>     log2(register count)            (default 13)
 *   --probminhash   -m   <int>     registers                       (default 1024)
 *                   -L   <int>     max-L truncation (0 = full)     (default 0)
 *   --bindash       --bbits <int>  bits per bin                    (default 16)
 *                   --sketch64 <int> 64-bin groups                 (default 16; NBINS = sketch64 * 64)
 *   --setsketch     -b   <int>     log2(register count)            (default 13)
 *                   -a   <float>   parameter a                     (default 5.0)
 *                   -B   <float>   base                            (default 2.0)
 *   --fastkmv       -K   <int>     sketch size (k smallest)        (default 1024)
 *
 * Output (TSV):  <fileA>\t<fileB>\t<distance>
 *   (Kssd --index uses an extended 5-column format identical to test_Kssd.cpp.)
 * Distance is the Mash distance D = -ln(2J/(1+J))/k for all algorithms.
 */

#include "sub_command.h"

#include <omp.h>

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>

// ─────────────────────────────────────────────────────────────────────────────
//  Usage
// ─────────────────────────────────────────────────────────────────────────────
static void print_usage(const char* prog) {
    std::cerr <<
"Usage: " << prog << " --<algo> [--index] -i <file>... [-l] [options]\n\n"
"Algorithms (choose one):\n"
"  --minhash --kssd --hll --probminhash --bindash --setsketch --fastkmv\n\n"
"Modes:\n"
"  -i fileA fileB         pairwise distance between two FASTA/FASTQ files\n"
"  -i list.txt -l         all-pairs over the file list\n"
"  --index                use inverted-index pipeline (list mode only;\n"
"                         supported for minhash/fastkmv/probminhash/kssd/setsketch/hll)\n\n"
"Common options:\n"
"  -o <path>              output file (default: rabbitsketch.dist)\n"
"  -k <int>               k-mer size (default 21; hll & setsketch use 32)\n"
"  -t <int>               threads (default: omp_get_max_threads())\n"
"  -d <float>             max Mash distance filter (default 0.05; use 1.0 for no filter)\n"
"  --seed <u64>           hash seed (default 42)\n\n"
"Algorithm-specific:\n"
"  --minhash      -s <int>             sketch size      (1000)\n"
"  --kssd         --drlevel <int>      (3)\n"
"  --hll          -b <int>             log2(registers)  (13)\n"
"  --probminhash  -m <int> -L <int>   sketch size / max-L (1024, 0)\n"
"                  freq-weighted PMH-norm (J_P estimator, coverage-robust)\n"
"  --bindash      --bbits <int> --sketch64 <int>  (16, 16)\n"
"  --setsketch    -b <int> -a <float> -B <float>  (13, 5.0, 2.0)\n"
"  --fastkmv      -K <int>             sketch size      (1024)\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  Argument parsing.  Populates the shared Args struct (defined in
//  sub_command.h) and validates mutual exclusion / required flags.
// ─────────────────────────────────────────────────────────────────────────────
static int parse_args(int argc, char* argv[], Args& a) {
    auto need = [&](int i, const char* opt) {
        if (i + 1 >= argc) {
            std::cerr << "ERROR: option " << opt << " requires an argument\n";
            std::exit(1);
        }
    };
    auto setAlgo = [&](Algo x) {
        if (a.algo != Algo::NONE) {
            std::cerr << "ERROR: more than one --<algo> flag specified\n";
            std::exit(1);
        }
        a.algo = x;
    };

    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if      (s == "--minhash")      setAlgo(Algo::MINHASH);
        else if (s == "--kssd")         setAlgo(Algo::KSSD);
        else if (s == "--hll")          setAlgo(Algo::HLL);
        else if (s == "--probminhash")  setAlgo(Algo::PROBMINHASH);
        else if (s == "--bindash")      setAlgo(Algo::BINDASH);
        else if (s == "--setsketch")    setAlgo(Algo::SETSKETCH);
        else if (s == "--fastkmv")      setAlgo(Algo::FASTKMV);
        else if (s == "--index")        a.useIndex = true;
        else if (s == "-l")             a.listMode = true;
        else if (s == "-i") {
            need(i, "-i");
            const size_t before = a.inputs.size();
            while (i + 1 < argc && argv[i + 1][0] != '-')
                a.inputs.emplace_back(argv[++i]);
            if (a.inputs.size() == before) {
                std::cerr << "ERROR: -i must be followed by at least one file "
                             "(got next token: '" << argv[i + 1] << "')\n";
                std::exit(1);
            }
        }
        else if (s == "-o")             { need(i, "-o"); a.output    = argv[++i]; }
        else if (s == "-k")             { need(i, "-k"); a.kmerSize  = std::stoi(argv[++i]); }
        else if (s == "-t")             { need(i, "-t"); a.threads   = std::max(1, std::stoi(argv[++i])); }
        else if (s == "-d")             { need(i, "-d"); a.maxDist   = std::stod(argv[++i]); }
        else if (s == "--seed")         { need(i, "--seed"); a.seed  = std::stoull(argv[++i]); }
        else if (s == "-s")             { need(i, "-s"); a.minhashSize = std::stoi(argv[++i]); }
        else if (s == "--drlevel")      { need(i, "--drlevel"); a.kssdDrlevel = std::stoi(argv[++i]); }
        else if (s == "-b")             { need(i, "-b"); a.hllBits = a.ssBits = std::stoi(argv[++i]); }
        else if (s == "-m")             { need(i, "-m"); a.pmhM = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (s == "-L")             { need(i, "-L"); a.pmhMaxL = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (s == "--bbits")        { need(i, "--bbits");    a.bdBbits    = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (s == "--sketch64")     { need(i, "--sketch64"); a.bdSketch64 = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (s == "-a")             { need(i, "-a"); a.ssA    = std::stod(argv[++i]); }
        else if (s == "-B")             { need(i, "-B"); a.ssBase = std::stod(argv[++i]); }
        else if (s == "-K")             { need(i, "-K"); a.fkmvK  = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (s == "-h" || s == "--help") { print_usage(argv[0]); std::exit(0); }
        else {
            std::cerr << "ERROR: unknown option: " << s << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (a.algo == Algo::NONE) {
        std::cerr << "ERROR: no --<algo> selected\n";
        print_usage(argv[0]);
        return 1;
    }
    if (a.inputs.empty()) {
        std::cerr << "ERROR: -i is required\n";
        return 1;
    }
    if (a.listMode) {
        if (a.inputs.size() != 1) {
            std::cerr << "ERROR: -l mode expects exactly one file (the list)\n";
            return 1;
        }
    } else {
        if (a.inputs.size() != 2) {
            std::cerr << "ERROR: pair mode requires exactly two files (or use -l)\n";
            return 1;
        }
    }
    if (a.useIndex && !a.listMode) {
        std::cerr << "WARNING: --index is only meaningful in list (-l) mode; ignored\n";
        a.useIndex = false;
    }
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Entry point.  All actual work happens in run_cli() in sub_command.cpp.
// ─────────────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    Args a;
    a.threads = omp_get_max_threads();
    if (parse_args(argc, argv, a) != 0) return 1;
    return run_cli(a);
}
