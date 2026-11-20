#ifndef RABBITSKETCH_POSTING_CODEC_H
#define RABBITSKETCH_POSTING_CODEC_H

#include <cstddef>
#include <cstdint>

namespace Sketch { namespace PostingCodec {

// Exact unsigned LEB128 deltas between nondecreasing sample IDs. A zero
// delta preserves repeated occurrences; sample 0 is not a sentinel.
inline uint32_t encodedSize(uint32_t value) {
    return 1 + (value >= (1u << 7)) + (value >= (1u << 14))
             + (value >= (1u << 21)) + (value >= (1u << 28));
}

inline uint32_t encode(uint8_t* output, uint32_t value) {
    uint8_t* start = output;
    while (value >= 128) {
        *output++ = static_cast<uint8_t>(value | 128);
        value >>= 7;
    }
    *output++ = static_cast<uint8_t>(value);
    return static_cast<uint32_t>(output - start);
}

// These readers consume trusted, internally built lists, not external files.
inline uint32_t decode(const uint8_t*& p) {
    uint32_t byte = *p++, value = byte & 127, shift = 7;
    while (byte & 128) {
        byte = *p++;
        value |= (byte & 127) << shift;
        shift += 7;
    }
    return value;
}

template<class Visit>
inline void forEachCount(const uint8_t* p, uint32_t count, Visit visit) {
    uint32_t sample = 0;
    for (uint32_t i = 0; i < count; ++i) {
        sample += decode(p);
        visit(sample);
    }
}

template<class Visit>
inline void forEachBytes(const uint8_t* p, size_t bytes, Visit visit) {
    const uint8_t* end = p + bytes;
    uint32_t sample = 0;
    while (p != end) {
        sample += decode(p);
        visit(sample);
    }
}

}} // namespace Sketch::PostingCodec
#endif
