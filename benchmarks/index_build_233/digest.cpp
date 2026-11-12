// Parallel equivalent of the existing three-column output fingerprint.
// Fingerprints check integrity; they are not a mathematical equality proof.
#define XXH_INLINE_ALL
#include "xxhash.h"
#include <omp.h>
#include <algorithm>
#include <cstring>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace {
struct Digest {
    uint64_t rows = 0, bytes = 0;
    unsigned __int128 sum = 0, xors = 0;
};
std::string hex128(unsigned __int128 value) {
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << uint64_t(value >> 64)
        << std::setw(16) << uint64_t(value);
    return out.str();
}
Digest scan(const char* data, size_t size, size_t lo, size_t hi) {
    Digest result;
    std::string reversed;
    while (lo > 0 && lo < size && data[lo - 1] != '\n') ++lo;
    while (lo < hi && lo < size) {
        const char* end = static_cast<const char*>(std::memchr(data + lo, '\n', size - lo));
        if (!end) throw std::runtime_error("truncated final line");
        const std::string_view line(data + lo, end - data - lo);
        const auto first = line.find('\t'), last = line.rfind('\t');
        if (first == std::string_view::npos || last == first || first == 0 || last == first + 1 ||
            line.find('\t', first + 1) != last)
            throw std::runtime_error("malformed three-column output");
        const auto a = line.substr(0, first), b = line.substr(first + 1, last - first - 1);
        const auto value = line.substr(last + 1);
        if (a == b || value.size() != 8 || value[0] != '0' || value[1] != '.')
            throw std::runtime_error("invalid pair or distance format");
        unsigned distance = 0;
        for (size_t j = 2; j < value.size(); ++j) {
            if (value[j] < '0' || value[j] > '9') throw std::runtime_error("invalid distance");
            distance = distance * 10 + unsigned(value[j] - '0');
        }
        if (distance > 50000) throw std::runtime_error("distance outside default cutoff");
        std::string_view canonical = line;
        if (b < a) {
            reversed.assign(b.data(), b.size());
            reversed += '\t'; reversed.append(a.data(), a.size());
            reversed.append(line.data() + last, line.size() - last);
            canonical = reversed;
        }
        const auto hash = XXH3_128bits(canonical.data(), canonical.size());
        const unsigned __int128 wide = (static_cast<unsigned __int128>(hash.high64) << 64) | hash.low64;
        result.sum += wide; result.xors ^= wide;
        ++result.rows; result.bytes += line.size() + 1;
        lo = end - data + 1;
    }
    return result;
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("usage: output_digest output.tsv threads");
        const int threads = std::stoi(argv[2]);
        if (threads < 1) throw std::runtime_error("invalid thread count");
        const int fd = open(argv[1], O_RDONLY);
        if (fd < 0) throw std::runtime_error("cannot open output");
        struct stat info{};
        if (fstat(fd, &info) || info.st_size < 0) throw std::runtime_error("cannot stat output");
        const size_t size = info.st_size;
        Digest result;
        if (size) {
            const char* data = static_cast<const char*>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
            if (data == MAP_FAILED) throw std::runtime_error("mmap failed");
            std::vector<Digest> partial(threads);
            std::string error;
            #pragma omp parallel for num_threads(threads) schedule(static)
            for (int part = 0; part < threads; ++part) {
                try {
                    partial[part] = scan(data, size, size * part / threads, size * (part + 1) / threads);
                } catch (const std::exception& exception) {
                    #pragma omp critical
                    error = exception.what();
                }
            }
            munmap(const_cast<char*>(data), size);
            if (!error.empty()) throw std::runtime_error(error);
            for (const auto& part : partial) {
                result.rows += part.rows; result.bytes += part.bytes;
                result.sum += part.sum; result.xors ^= part.xors;
            }
        }
        close(fd);
        if (result.bytes != size) throw std::runtime_error("digest missed or repeated input bytes");
        std::cout << "{\"rows\":" << result.rows << ",\"bytes\":" << result.bytes
                  << ",\"canonical_xxh3_128_sum\":\"" << hex128(result.sum)
                  << "\",\"canonical_xxh3_128_xor\":\"" << hex128(result.xors) << "\"}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
