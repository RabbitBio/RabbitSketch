#ifndef RABBITSKETCH_RANK_METADATA_H
#define RABBITSKETCH_RANK_METADATA_H

#include "rank/CanonicalKmer.h"
#include "rank/RankStream.h"

#include <cstdint>
#include <array>
#include <cstddef>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Sketch {
namespace Rank {

enum class Backend : uint8_t {
    FastKMV = 1,
    KSSD = 2,
    ProbMinHash = 3,
    SetSketch = 4,
    HLL = 5,
    BinDash = 6,
    LegacyMinHash = 7,
    FracMinHash = 8,
    OrderMinHash = 9
};

enum class HashProfile : uint8_t {
    CanonicalBits = 0,
    Murmur3Fmix64V1 = 1,
    Murmur3Fmix64Twice = 2,
    LegacyMurmur3X64V1 = 3,
    BinDashDoubleFmixV1 = 4,
    ProbMinHashRaceV1 = 5,
    FracRankStreamV1 = 6,
    OrderMinHashV1 = 7
};

enum class WeightSemantics : uint8_t {
    UnweightedSet = 1,
    Frequency = 2,
    UserSupplied = 3
};

enum class ResolutionKind : uint8_t {
    BottomK = 1,
    RegisterPrecisionBits = 2,
    KssdDrLevel = 3,
    RegisterCount = 4,
    Scaled = 5,
    OrderSampleCount = 6
};

struct RankMetadata {
    // Little-endian wire bytes are the ASCII string "RSK2".
    static constexpr uint32_t MAGIC = UINT32_C(0x324b5352);
    static constexpr uint16_t SCHEMA_VERSION = 1;
    static constexpr size_t WIRE_SIZE = 40;

    uint32_t magic = MAGIC;
    uint16_t schema_version = SCHEMA_VERSION;
    uint16_t rank_stream_version = RankStream::VERSION;
    Canonicalization canonicalization =
        Canonicalization::Lexicographic2BitReverseComplementV1;
    HashProfile hash_profile = HashProfile::Murmur3Fmix64V1;
    Backend backend = Backend::FastKMV;
    WeightSemantics weight_semantics = WeightSemantics::UnweightedSet;
    ResolutionKind resolution_kind = ResolutionKind::BottomK;
    uint8_t fingerprint_bits = 64;
    uint16_t kmer_size = 0;
    uint32_t resolution = 0;
    uint64_t seed = 42;
    uint64_t parameter_fingerprint = 0;

    /** Fixed-width little-endian header; never serialize the C++ struct raw. */
    std::array<uint8_t, WIRE_SIZE> toBytes() const noexcept {
        std::array<uint8_t, WIRE_SIZE> bytes{};
        writeLittleEndian(bytes.data() + 0, magic, 4);
        writeLittleEndian(bytes.data() + 4, schema_version, 2);
        writeLittleEndian(bytes.data() + 6, rank_stream_version, 2);
        writeLittleEndian(bytes.data() + 8,
                          static_cast<uint16_t>(canonicalization), 2);
        bytes[10] = static_cast<uint8_t>(hash_profile);
        bytes[11] = static_cast<uint8_t>(backend);
        bytes[12] = static_cast<uint8_t>(weight_semantics);
        bytes[13] = static_cast<uint8_t>(resolution_kind);
        bytes[14] = fingerprint_bits;
        writeLittleEndian(bytes.data() + 16, kmer_size, 2);
        writeLittleEndian(bytes.data() + 20, resolution, 4);
        writeLittleEndian(bytes.data() + 24, seed, 8);
        writeLittleEndian(bytes.data() + 32, parameter_fingerprint, 8);
        return bytes;
    }

    static RankMetadata fromBytes(const uint8_t* bytes, size_t size) {
        if (bytes == nullptr || size < WIRE_SIZE) {
            throw std::invalid_argument("RankMetadata header is truncated");
        }
        RankMetadata meta;
        meta.magic = static_cast<uint32_t>(readLittleEndian(bytes + 0, 4));
        meta.schema_version = static_cast<uint16_t>(
            readLittleEndian(bytes + 4, 2));
        meta.rank_stream_version = static_cast<uint16_t>(
            readLittleEndian(bytes + 6, 2));
        meta.canonicalization = static_cast<Canonicalization>(
            readLittleEndian(bytes + 8, 2));
        meta.hash_profile = static_cast<HashProfile>(bytes[10]);
        meta.backend = static_cast<Backend>(bytes[11]);
        meta.weight_semantics = static_cast<WeightSemantics>(bytes[12]);
        meta.resolution_kind = static_cast<ResolutionKind>(bytes[13]);
        meta.fingerprint_bits = bytes[14];
        meta.kmer_size = static_cast<uint16_t>(
            readLittleEndian(bytes + 16, 2));
        meta.resolution = static_cast<uint32_t>(
            readLittleEndian(bytes + 20, 4));
        meta.seed = readLittleEndian(bytes + 24, 8);
        meta.parameter_fingerprint = readLittleEndian(bytes + 32, 8);
        if (meta.magic != MAGIC) {
            throw std::invalid_argument("RankMetadata header has invalid magic");
        }
        if (meta.schema_version != SCHEMA_VERSION) {
            throw std::invalid_argument("unsupported RankMetadata schema version");
        }
        if (meta.rank_stream_version != RankStream::VERSION) {
            throw std::invalid_argument("unsupported RankStream version");
        }
        const uint16_t canonical_value =
            static_cast<uint16_t>(meta.canonicalization);
        const uint8_t hash_value = static_cast<uint8_t>(meta.hash_profile);
        const uint8_t backend_value = static_cast<uint8_t>(meta.backend);
        const uint8_t weight_value = static_cast<uint8_t>(meta.weight_semantics);
        const uint8_t resolution_kind_value =
            static_cast<uint8_t>(meta.resolution_kind);
        const bool valid_zero_resolution =
            meta.backend == Backend::KSSD &&
            meta.resolution_kind == ResolutionKind::KssdDrLevel;
        if (canonical_value < 1 || canonical_value > 2 || hash_value > 7 ||
            backend_value < 1 || backend_value > 9 ||
            weight_value < 1 || weight_value > 3 ||
            resolution_kind_value < 1 || resolution_kind_value > 6 ||
            meta.fingerprint_bits == 0 || meta.fingerprint_bits > 64 ||
            meta.kmer_size == 0 || meta.kmer_size > 32 ||
            (meta.resolution == 0 && !valid_zero_resolution)) {
            throw std::invalid_argument(
                "RankMetadata header contains invalid field values");
        }
        return meta;
    }

    std::string coordinationMismatch(const RankMetadata& other) const {
        if (magic != MAGIC || other.magic != MAGIC) return "invalid metadata magic";
        if (schema_version != SCHEMA_VERSION ||
            other.schema_version != SCHEMA_VERSION ||
            schema_version != other.schema_version) {
            return "metadata schema version";
        }
        if (rank_stream_version != RankStream::VERSION ||
            other.rank_stream_version != RankStream::VERSION ||
            rank_stream_version != other.rank_stream_version) {
            return "rank stream version";
        }
        if (canonicalization != other.canonicalization) return "canonicalization";
        if (hash_profile != other.hash_profile) return "hash profile";
        if (kmer_size != other.kmer_size) return "k-mer size";
        if (seed != other.seed) return "seed";
        if (weight_semantics != other.weight_semantics) return "weight semantics";
        return std::string();
    }

    bool coordinatedWith(const RankMetadata& other) const {
        return coordinationMismatch(other).empty();
    }

    std::string sameFamilyMismatch(const RankMetadata& other,
                                   bool allow_resolution_change) const {
        const std::string rank_mismatch = coordinationMismatch(other);
        if (!rank_mismatch.empty()) return rank_mismatch;
        if (backend != other.backend) return "backend";
        if (parameter_fingerprint != other.parameter_fingerprint) {
            return "backend parameters";
        }
        if (resolution_kind != other.resolution_kind) return "resolution kind";
        if (!allow_resolution_change && resolution != other.resolution) {
            return "resolution";
        }
        return std::string();
    }

    void requireCoordinated(const RankMetadata& other,
                            const char* operation) const {
        const std::string mismatch = coordinationMismatch(other);
        if (!mismatch.empty()) throwMismatch(operation, mismatch);
    }

    void requireSameFamily(const RankMetadata& other,
                           const char* operation,
                           bool allow_resolution_change = false) const {
        const std::string mismatch =
            sameFamilyMismatch(other, allow_resolution_change);
        if (!mismatch.empty()) throwMismatch(operation, mismatch);
    }

private:
    static void writeLittleEndian(uint8_t* destination,
                                  uint64_t value,
                                  size_t width) noexcept {
        for (size_t i = 0; i < width; ++i) {
            destination[i] = static_cast<uint8_t>(value >> (8 * i));
        }
    }

    static uint64_t readLittleEndian(const uint8_t* source,
                                     size_t width) noexcept {
        uint64_t value = 0;
        for (size_t i = 0; i < width; ++i) {
            value |= static_cast<uint64_t>(source[i]) << (8 * i);
        }
        return value;
    }

    [[noreturn]] static void throwMismatch(const char* operation,
                                           const std::string& mismatch) {
        std::ostringstream message;
        message << operation << ": incompatible sketches (" << mismatch << ')';
        throw std::invalid_argument(message.str());
    }
};

inline uint64_t doubleBits(double value) noexcept {
    uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "unexpected double width");
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

inline uint64_t parameterFingerprint(double first,
                                     double second,
                                     uint64_t integer_parameter = 0) noexcept {
    const uint64_t a = RankStream::fmix64(doubleBits(first),
                                          UINT64_C(0x6a09e667f3bcc909));
    const uint64_t b = RankStream::fmix64(doubleBits(second),
                                          UINT64_C(0xbb67ae8584caa73b));
    return RankStream::fmix64(a ^ b ^ integer_parameter,
                              UINT64_C(0x3c6ef372fe94f82b));
}

} // namespace Rank
} // namespace Sketch

#endif // RABBITSKETCH_RANK_METADATA_H
