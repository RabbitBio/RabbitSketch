// Streaming, order-independent validation of the actual CLI's distance output.
// XXH3 is an integrity fingerprint, not a cryptographic proof of equality.
#define XXH_INLINE_ALL
#include "xxhash.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static std::string hex128(unsigned __int128 value) {
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16)
        << static_cast<uint64_t>(value >> 64) << std::setw(16)
        << static_cast<uint64_t>(value);
    return out.str();
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("usage: output_digest distance.tsv");
        std::vector<char> buffer(1 << 22);
        std::ifstream input;
        input.rdbuf()->pubsetbuf(buffer.data(), buffer.size());
        input.open(argv[1], std::ios::binary);
        if (!input) throw std::runtime_error("cannot open distance output");
        uint64_t rows = 0, bytes = 0;
        unsigned __int128 sum = 0, xors = 0, field_sum = 0, field_xors = 0;
        std::string line, canonical;
        while (std::getline(input, line)) {
            if (input.eof()) throw std::runtime_error("truncated final line");
            bytes += line.size() + 1;
            if (rows == 0 && line == "genome0\tgenome1\tcommon|size0|size1\tjaccard\tmashD") continue;
            const size_t first = line.find('\t');
            const size_t second = line.find('\t', first == std::string::npos ? first : first + 1);
            const size_t last = line.rfind('\t');
            const size_t columns = 1 + std::count(line.begin(), line.end(), '\t');
            if (first == std::string::npos || second == std::string::npos || first == 0 ||
                second == first + 1 || last + 1 == line.size() || (columns != 3 && columns != 5))
                throw std::runtime_error("malformed three/five-column TSV");
            const std::string a = line.substr(0, first), b = line.substr(first + 1, second - first - 1);
            if (a == b) throw std::runtime_error("unexpected self pair");
            char* end = nullptr;
            const double distance = std::strtod(line.c_str() + last + 1, &end);
            if (*end || !std::isfinite(distance) || distance < 0 || distance > 0.0500005)
                throw std::runtime_error("distance outside the default output filter");
            if (a < b) canonical = a + '\t' + b + line.substr(last);
            else canonical = b + '\t' + a + line.substr(last);
            const XXH128_hash_t h = XXH3_128bits(canonical.data(), canonical.size());
            const unsigned __int128 h128 = (static_cast<unsigned __int128>(h.high64) << 64) | h.low64;
            sum += h128; xors ^= h128;
            // Both indexed implementations emit the same file-list orientation.
            // Also fingerprint all estimator columns, not just distance.
            const auto fields_hash = XXH3_128bits(line.data(), line.size());
            const unsigned __int128 fields = (static_cast<unsigned __int128>(fields_hash.high64) << 64) | fields_hash.low64;
            field_sum += fields; field_xors ^= fields;
            ++rows;
        }
        if (input.bad()) throw std::runtime_error("error reading distance output");
        std::cout << "{\"rows\":" << rows << ",\"bytes\":" << bytes
                  << ",\"canonical_xxh3_128_sum\":\"" << hex128(sum)
                  << "\",\"canonical_xxh3_128_xor\":\"" << hex128(xors) << "\",\"all_fields_sum\":\"" << hex128(field_sum)
                  << "\",\"all_fields_xor\":\"" << hex128(field_xors) << "\"}\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
