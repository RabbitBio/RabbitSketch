/**
 * sub_command.h – API surface for the rabbitsketch CLI sub-commands.
 *
 *   main.cpp            : argument parsing + dispatch only
 *   sub_command.{h,cpp} : Args struct definition + every algorithm's
 *                         actual sketch / distance pipeline.  Each
 *                         indexed runner mirrors its corresponding
 *                         test_*.cpp byte-for-byte.
 */

#ifndef RABBITSKETCH_SUB_COMMAND_H
#define RABBITSKETCH_SUB_COMMAND_H

#include <cstdint>
#include <string>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
//  Algorithm selector (mutually exclusive --<algo> flag in CLI)
// ─────────────────────────────────────────────────────────────────────────────
enum class Algo {
    NONE,
    MINHASH,
    KSSD,
    HLL,
    PROBMINHASH,
    BINDASH,
    SETSKETCH,
    FASTKMV,
};

// ─────────────────────────────────────────────────────────────────────────────
//  Parsed CLI arguments – shared between main.cpp (parser) and
//  sub_command.cpp (executor).
// ─────────────────────────────────────────────────────────────────────────────
struct Args {
    Algo        algo        = Algo::NONE;
    bool        useIndex    = false;
    bool        listMode    = false;
    std::vector<std::string> inputs;
    std::string output      = "rabbitsketch.dist";

    int         kmerSize    = 21;
    int         threads     = 1;     // overwritten by main() with omp_get_max_threads()
    double      maxDist     = 0.05;  // typical genomic Mash-distance cutoff
    uint64_t    seed        = 42;

    int         minhashSize = 1000;
    int         kssdDrlevel = 3;
    int         hllBits     = 13;

    uint32_t    pmhM        = 1024;
    uint32_t    pmhMaxL     = 0;

    uint32_t    bdBbits     = 16;
    uint32_t    bdSketch64  = 16;   // matches test_BinDash.cpp (NBINS = 16 * 64 = 1024)

    int         ssBits      = 13;
    double      ssA         = 5.0;
    double      ssBase      = 2.0;

    uint32_t    fkmvK       = 1024;
};

// ─────────────────────────────────────────────────────────────────────────────
//  Single CLI dispatch entry point.
//
//  Picks the right runner based on (a.algo, a.listMode, a.useIndex) and
//  returns a process exit code (0 = success, non-zero = error).
// ─────────────────────────────────────────────────────────────────────────────
int run_cli(const Args& a);

#endif  // RABBITSKETCH_SUB_COMMAND_H
