# RabbitSketch 2

RabbitSketch is a high-performance C++17/Python library for genomic sketching.
It provides one streaming FASTA/FASTQ build contract and one query result model
across nine algorithms while preserving each algorithm's native estimator.

The high-level API focuses on in-memory workflows: configure one or more
algorithms, read the input once, receive typed `BuiltSketch` objects, and query
compatible sketches. Unsupported comparisons return a reason instead of being
silently interpreted as zero similarity.

## Algorithms

| Algorithm | Sampling model | Primary queries | Exact common resolution |
|---|---|---|---|
| Legacy MinHash | bottom-k | set metrics, ANI | smaller K |
| FastKMV | bottom-k | set metrics, ANI | smaller K |
| FracMinHash | scaled | set metrics, ANI | larger scaled value |
| HyperLogLog | register precision | cardinality-derived set metrics | lower precision |
| SetSketch | register precision | cardinality-derived set metrics | same precision only |
| KSSD | reduction level | Jaccard, distance, ANI | sparser reduction |
| ProbMinHash | native registers | set or weighted metrics | same register count |
| BinDash | b-bit bins | set metrics, ANI | same bin count |
| OrderMinHash | ordered samples | order similarity/distance | shorter prefix |

The FASTX builder accepts DNA and RNA from plain or gzip FASTA/FASTQ, standard
input, or in-memory records. It supports multiline records, canonical k-mers,
ambiguous-base policies, Phred filtering, homopolymer compression,
frequency-aware ProbMinHash, and file/record/collection aggregation. Multiple
algorithms consume the same normalized input pass.

Queries are same-backend by design. The planner validates k-mer size, seed,
canonicalization, hash profile, weight semantics, backend parameters, and
resolution before selecting a native estimator.

## Build and install

Requirements: CMake 3.16+, a C++17 compiler, zlib, and optionally OpenMP.
Linux is the tested deployment target.

```bash
cmake -S . -B build \
  -DCXXAPI=ON \
  -DRABBITSKETCH_BUILD_TESTS=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
cmake --install build --prefix /your/prefix
```

Installed CMake consumers can use:

```cmake
find_package(RabbitSketch 2 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE RabbitSketch::rabbitsketch)
```

The install also provides `rabbitsketch.pc`. Python 3.8+ uses an isolated
PEP 517 build and has no Python runtime dependency on NumPy, `fastx`, or `pymp`:

```bash
python -m pip install .
```

## Python example

```python
import rabbitsketch as rs

fast = rs.SketchConfig()
fast.algorithm = rs.Algorithm.FastKMV
fast.sampling = rs.SamplingMode.BottomK
fast.kmer_size = 21
fast.resolution = 1024

frac = rs.SketchConfig()
frac.algorithm = rs.Algorithm.FracMinHash
frac.sampling = rs.SamplingMode.Scaled
frac.kmer_size = 21
frac.scaled = 1000
frac.aggregation = rs.AggregationMode.OneSketchPerRecord

results = rs.build_fastx("sample.fq.gz", [fast, frac])
query = results[0].sketch.query(results[0].sketch)
assert query.estimate_available and query.jaccard == 1.0
```

`FastxReader` is an iterator, and `MultiSketchBuilder` accepts in-memory
`FastxRecord` objects. Expensive native updates, builds, and queries release the
Python GIL.

## C++ example

```cpp
#include "api/SketchBuilder.h"

Sketch::API::SketchConfig config;
config.algorithm = Sketch::API::Algorithm::FastKMV;
config.sampling = Sketch::API::SamplingMode::BottomK;
config.kmer_size = 21;
config.resolution = 1024;

auto built = Sketch::API::buildFastx("sample.fa.gz", {config});
auto result = built.at(0).sketch.query(built.at(0).sketch);
if (!result.estimate_available)
    throw std::runtime_error(result.plan.reason);
```

Low-level algorithm classes remain available for applications that need direct
updates, projections, merges, or native query methods.

## Correctness and portability

- A rejected query has `estimate_available == false`; inspect `plan.reason`.
- Weighted ProbMinHash and unweighted set sketches have different estimands.
- OrderMinHash exposes order metrics, not set Jaccard.
- ANI and Mash distance are model-derived from k-mer similarity, not alignments.
- SetSketch precision folding is rejected until a valid native rule exists.

Release builds are portable by default and do not use `-march=native`. Hot
comparison loops select scalar, SSE2, AVX2, or AVX-512 implementations at
runtime. Python exposes the selection through `runtime_info()`;
`RABBITSKETCH_SIMD=scalar` forces the scalar path for diagnosis. Local builds
may opt into `-DRABBITSKETCH_NATIVE_ARCH=ON`.

More detail is available in the [practical API guide](docs/rabbitsketch2_practical_api.md),
[unified query contract](docs/rabbitsketch2_unified_query.md), and
[compatibility matrix](docs/rabbitsketch2_compatibility_matrix.md).

## License

RabbitSketch is distributed under the MIT License; see [LICENSE.txt](LICENSE.txt).
