#include "io/FastxReader.h"

#include "kseq.h"

#include <stdexcept>
#include <utility>

#include <unistd.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

namespace Sketch {
namespace IO {

struct FastxReader::Impl {
    gzFile file = nullptr;
    kseq_t* sequence = nullptr;

    ~Impl() {
        if (sequence != nullptr) kseq_destroy(sequence);
        if (file != nullptr) gzclose(file);
    }
};

FastxReader::FastxReader(std::string path)
    : path_(std::move(path)), impl_(new Impl()) {
    if (path_.empty())
        throw std::invalid_argument("FASTX input path is empty");
    if (path_ == "-") {
        const int descriptor = ::dup(STDIN_FILENO);
        if (descriptor < 0)
            throw std::runtime_error("failed to duplicate standard input");
        impl_->file = gzdopen(descriptor, "rb");
        if (impl_->file == nullptr) ::close(descriptor);
    } else {
        impl_->file = gzopen(path_.c_str(), "rb");
    }
    if (impl_->file == nullptr)
        throw std::runtime_error("failed to open FASTX input: " + path_);
    impl_->sequence = kseq_init(impl_->file);
    if (impl_->sequence == nullptr)
        throw std::runtime_error("failed to initialize FASTX parser: " + path_);
}

FastxReader::~FastxReader() = default;
FastxReader::FastxReader(FastxReader&&) noexcept = default;
FastxReader& FastxReader::operator=(FastxReader&&) noexcept = default;

bool FastxReader::next(FastxRecord& record) {
    if (!impl_ || impl_->sequence == nullptr)
        throw std::logic_error("FASTX reader is not initialized");
    const int status = kseq_read(impl_->sequence);
    if (status == -1) return false;
    if (status < -1)
        throw std::invalid_argument(
            "truncated or malformed FASTQ record in: " + path_);
    kseq_t* sequence = impl_->sequence;
    record.index = stats_.records;
    record.name.assign(sequence->name.s == nullptr ? "" : sequence->name.s,
                       sequence->name.l);
    record.comment.assign(
        sequence->comment.s == nullptr ? "" : sequence->comment.s,
        sequence->comment.l);
    record.sequence.assign(sequence->seq.s == nullptr ? "" : sequence->seq.s,
                           sequence->seq.l);
    record.quality.assign(sequence->qual.s == nullptr ? "" : sequence->qual.s,
                          sequence->qual.l);
    if (!record.quality.empty() &&
        record.quality.size() != record.sequence.size())
        throw std::invalid_argument(
            "FASTQ sequence/quality lengths differ in record: " + record.name);
    const FastxFormat current = record.quality.empty()
        ? FastxFormat::Fasta : FastxFormat::Fastq;
    if (format_ == FastxFormat::Unknown) format_ = current;
    else if (format_ != current)
        throw std::invalid_argument(
            "FASTX stream mixes FASTA and FASTQ records: " + path_);
    ++stats_.records;
    stats_.bases += record.sequence.size();
    stats_.quality_bases += record.quality.size();
    return true;
}

const char* fastxFormatName(FastxFormat format) noexcept {
    switch (format) {
        case FastxFormat::Unknown: return "unknown";
        case FastxFormat::Fasta: return "fasta";
        case FastxFormat::Fastq: return "fastq";
    }
    return "unknown";
}

} // namespace IO
} // namespace Sketch
