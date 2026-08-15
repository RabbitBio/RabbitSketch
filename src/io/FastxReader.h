#ifndef RABBITSKETCH_IO_FASTX_READER_H
#define RABBITSKETCH_IO_FASTX_READER_H

#include <cstdint>
#include <memory>
#include <string>

namespace Sketch {
namespace IO {

enum class FastxFormat : uint8_t {
    Unknown = 0,
    Fasta = 1,
    Fastq = 2
};

struct FastxRecord {
    uint64_t index = 0;
    std::string name;
    std::string comment;
    std::string sequence;
    std::string quality;

    bool hasQuality() const noexcept { return !quality.empty(); }
};

struct FastxReadStats {
    uint64_t records = 0;
    uint64_t bases = 0;
    uint64_t quality_bases = 0;
};

/** Streaming FASTA/FASTQ reader.  zlib transparent mode accepts plain or
 * gzip-compressed input; path "-" reads standard input. */
class FastxReader {
public:
    explicit FastxReader(std::string path);
    ~FastxReader();
    FastxReader(FastxReader&&) noexcept;
    FastxReader& operator=(FastxReader&&) noexcept;
    FastxReader(const FastxReader&) = delete;
    FastxReader& operator=(const FastxReader&) = delete;

    bool next(FastxRecord& record);
    const std::string& path() const noexcept { return path_; }
    FastxFormat format() const noexcept { return format_; }
    const FastxReadStats& stats() const noexcept { return stats_; }

private:
    struct Impl;
    std::string path_;
    std::unique_ptr<Impl> impl_;
    FastxFormat format_ = FastxFormat::Unknown;
    FastxReadStats stats_;
};

const char* fastxFormatName(FastxFormat format) noexcept;

} // namespace IO
} // namespace Sketch

#endif // RABBITSKETCH_IO_FASTX_READER_H
