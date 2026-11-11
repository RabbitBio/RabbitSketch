// Keep the historical positional interface while sharing the optimized CLI.
#include "sub_command.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <file.list> <threshold> <threads> [output_file]\n";
        return 1;
    }
    try {
        Args args;
        args.algo = Algo::KSSD;
        args.listMode = true;
        args.useIndex = true;
        args.inputs = {argv[1]};
        args.maxDist = std::stod(argv[2]);
        args.threads = std::stoi(argv[3]);
        args.output = argc == 5 ? argv[4] : "result.sketch.dist";
        if (args.threads < 1 || !std::isfinite(args.maxDist) || args.maxDist < 0)
            throw std::invalid_argument("threads must be positive and threshold finite and nonnegative");
        return run_cli(args);
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
