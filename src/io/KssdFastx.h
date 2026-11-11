#ifndef RABBITSKETCH_KSSD_FASTX_H
#define RABBITSKETCH_KSSD_FASTX_H

#include "KssdBuild.h"
#include "kseq.h"
#include <zlib.h>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>

namespace Sketch {
namespace KssdBuild {
namespace FastxDetail {
KSEQ_INIT(gzFile, gzread)
}

// Per-worker input workspace. Plain FASTA spans go directly to the encoder;
// FASTQ and gzip use kseq. Format sniffing reads only two bytes, not the entire
// compressed file. No k-mer ever spans a record boundary.
class FastxReader {
public:
    template<typename Emit>
    bool read(const std::string& path, Stream& stream, Emit&& emit) {
        FILE* input = std::fopen(path.c_str(), "rb");
        if (!input) return false;
        unsigned char magic[2] = {};
        const size_t got = std::fread(magic, 1, 2, input);
        if (got == 0 && !std::ferror(input)) { std::fclose(input); return true; }
        if (got == 0) { std::fclose(input); return false; }
        const bool gzip = got == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
        if (gzip || magic[0] != '>') {
            std::fclose(input);
            gzFile compressed = gzopen(path.c_str(), "rb");
            if (!compressed) return false;
            gzbuffer(compressed, 1 << 18);
            auto* record = FastxDetail::kseq_init(compressed);
            int status;
            while ((status = FastxDetail::kseq_read(record)) >= 0) {
                stream.reset();
                stream.append(record->seq.s, record->seq.l, emit);
            }
            FastxDetail::kseq_destroy(record);
            int gzip_error = Z_OK;
            gzerror(compressed, &gzip_error);
            const int closed = gzclose(compressed);
            return status == -1 && (gzip_error == Z_OK || gzip_error == Z_STREAM_END)
                && closed == Z_OK;
        }
        struct stat info;
        if (fstat(fileno(input), &info) != 0 || info.st_size < 0 ||
            std::fseek(input, 0, SEEK_SET) != 0) {
            std::fclose(input); return false;
        }
        bytes_.resize(static_cast<size_t>(info.st_size));
        const size_t count = std::fread(bytes_.data(), 1, bytes_.size(), input);
        const bool okay = count == bytes_.size() && !std::ferror(input);
        std::fclose(input);
        if (!okay) return false;
        size_t pos = 0;
        while (pos < count) {
            const char* newline = static_cast<const char*>(
                std::memchr(bytes_.data() + pos, '\n', count - pos));
            if (!newline) break; // A final header without sequence is empty.
            const size_t begin = static_cast<size_t>(newline - bytes_.data()) + 1;
            size_t end = begin;
            while (end < count) {
                const char* next = static_cast<const char*>(
                    std::memchr(bytes_.data() + end, '>', count - end));
                if (!next) { end = count; break; }
                end = static_cast<size_t>(next - bytes_.data());
                if (end == 0 || bytes_[end - 1] == '\n' || bytes_[end - 1] == '\r') break;
                ++end;
            }
            stream.reset();
            stream.append<true>(bytes_.data() + begin, end - begin, emit);
            pos = end;
        }
        return true;
    }
private:
    std::vector<char> bytes_;
};

} // namespace KssdBuild
} // namespace Sketch
#endif
