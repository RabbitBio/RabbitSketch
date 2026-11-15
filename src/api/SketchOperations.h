#ifndef RABBITSKETCH_API_SKETCH_OPERATIONS_H
#define RABBITSKETCH_API_SKETCH_OPERATIONS_H

#include "api/SketchIndex.h"
#include "Sketch.h"
#include "BinDash.h"
#include "SetSketch.h"
#include "fastkmv.h"
#include "fracminhash.h"
#include "probmh.h"

namespace Sketch { namespace API {

inline BuiltSketch toBuiltSketch(const BuiltSketch& value) { return value; }
inline BuiltSketch toBuiltSketch(const FastKMV& value) { return BuiltSketch::fromFastKMV(value); }
inline BuiltSketch toBuiltSketch(const FracMinHash& value) { return BuiltSketch::fromFracMinHash(value); }
inline BuiltSketch toBuiltSketch(const HyperLogLog& value) { return BuiltSketch::fromHyperLogLog(value); }
inline BuiltSketch toBuiltSketch(const SetSketch& value) { return BuiltSketch::fromSetSketch(value); }
inline BuiltSketch toBuiltSketch(const BinDash& value) { return BuiltSketch::fromBinDash(value); }
inline BuiltSketch toBuiltSketch(const ProbMinHash4& value) { return BuiltSketch::fromProbMinHash(value); }
inline BuiltSketch toBuiltSketch(const Kssd& value) {
    return BuiltSketch::fromKssd(std::make_unique<Kssd>(value.project(value.get_drlevel())));
}

/** Native collections use the same CSR implementation as BuiltSketch. */
template<class T>
SketchIndex indexDict(const std::vector<T>& sketches, int num_threads = 1) {
    std::vector<BuiltSketch> snapshots;
    snapshots.reserve(sketches.size());
    for (const auto& sketch : sketches) snapshots.push_back(toBuiltSketch(sketch));
    return indexDict(snapshots, num_threads);
}

}} // namespace Sketch::API
#endif
