#!/usr/bin/env python3
"""End-to-end smoke test for the public Python FASTX build surface."""

import gzip
import pathlib
import random
import sys
import tempfile
import threading
import time


def check_setsketch(rs, directory):
    """The FASTX builder and native binding must share the corrected estimator."""
    rng = random.Random(20261105)
    small = "".join(rng.choices("ACGT", k=10_000))
    large = small[:130] + "".join(rng.choices("ACGT", k=99_870))
    path = directory / "setsketch.fa"
    path.write_text(f">small\n{small}\n>large\n{large}\n", encoding="ascii")
    config = rs.SketchConfig()
    config.algorithm = rs.Algorithm.SetSketch
    config.sampling = rs.SamplingMode.RegisterPrecision
    config.aggregation = rs.AggregationMode.OneSketchPerRecord
    config.kmer_size = 20
    config.resolution = 16
    built = rs.build_fastx(str(path), [config])
    left, right = (rs.SetSketch(16, 2.0, 5.0, 20, 42) for _ in range(2))
    left.update(small)
    right.update(large)
    native = left.jaccard(right)
    result = built[0].sketch.query(built[1].sketch)
    assert result.estimate_available
    assert result.plan.method == rs.QueryMethod.SetSketchNative
    assert abs(result.jaccard - native) < 1e-12
    assert abs(rs.query(left, right).jaccard - native) < 1e-12
    assert abs(left.distance(right) - (1.0 - native)) < 1e-12
    assert left.jaccard(left) == 1.0

    complement = str.maketrans("ACGT", "TGCA")
    def kmers(sequence):
        return {min(kmer, kmer.translate(complement)[::-1])
                for i in range(len(sequence) - 19)
                for kmer in [sequence[i:i + 20]]}
    a, b = kmers(small), kmers(large)
    truth = len(a & b) / len(a | b)
    assert abs(native - truth) < 0.005, "SetSketch sparse-state bias regression"


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_fastx_python.py BUILD_DIRECTORY")
    sys.path.insert(0, sys.argv[1])
    import rabbitsketch as rs

    assert rs.__version__ == "2.0.0"
    runtime = rs.runtime_info()
    assert runtime.architecture
    assert runtime.portable_baseline
    assert runtime.selected_byte_path in {"scalar", "sse2", "avx2", "avx512bw"}
    assert runtime.selected_rank_path in {"scalar", "avx2", "avx512dq"}

    with tempfile.TemporaryDirectory(prefix="rabbitsketch-fastx-python-") as tmp:
        check_setsketch(rs, pathlib.Path(tmp))
        path = pathlib.Path(tmp) / "reads.fa.gz"
        with gzip.open(path, "wt", encoding="ascii") as output:
            output.write(">alpha comment\nACGTACGTACGT\n>beta\nTTTTACGTACGT\n")

        reader = rs.FastxReader(str(path))
        records = list(reader)
        assert [record.name for record in records] == ["alpha", "beta"]
        assert records[0].comment == "comment"
        assert reader.format == rs.FastxFormat.Fasta
        assert reader.stats.records == 2
        assert reader.stats.bases == 24
        assert reader.read() is None

        fast = rs.SketchConfig()
        fast.kmer_size = 5
        fast.resolution = 64
        frac = rs.SketchConfig()
        frac.algorithm = rs.Algorithm.FracMinHash
        frac.sampling = rs.SamplingMode.Scaled
        frac.kmer_size = 5
        frac.scaled = 2
        frac.aggregation = rs.AggregationMode.OneSketchPerRecord

        built = rs.build_fastx(str(path), [fast, frac])
        assert len(built) == 3
        assert built[0].stats.records == 2
        assert built[0].stats.accepted_kmers == 16
        assert built[0].label
        assert built[0].sketch.algorithm == rs.Algorithm.FastKMV
        assert [item.record_name for item in built[1:]] == ["alpha", "beta"]
        comparison = built[0].sketch.query(built[0].sketch)
        assert comparison.jaccard == 1.0

        native_frac = rs.FracMinHash(2, 5, 42)
        native_frac.update("ACGTACGTACGT")
        native_frac.finalize()
        assert native_frac.cardinality_at(4) >= 0.0

        builder = rs.MultiSketchBuilder([fast])
        builder.update(records[0])
        manual = builder.finish("manual", str(path), records[0].name)
        assert builder.finished and builder.lanes == 1
        assert manual[0].stats.records == 1

        # A second Python thread must make progress during a heavy native
        # update; this catches accidental removal of the GIL release guard.
        ready = threading.Event()
        stop = threading.Event()
        ticks = [0]

        def ticker():
            ready.set()
            while not stop.is_set():
                ticks[0] += 1

        thread = threading.Thread(target=ticker)
        thread.start()
        ready.wait()
        time.sleep(0.01)
        before = ticks[0]
        native = rs.FastKMV(256, 21, 42)
        native.update("ACGT" * 1_000_000)
        after = ticks[0]
        stop.set()
        thread.join()
        assert after > before, "FastKMV.update held the Python GIL"

        order_config = rs.SketchConfig()
        order_config.algorithm = rs.Algorithm.OrderMinHash
        order_config.sampling = rs.SamplingMode.AlgorithmNative
        order_config.aggregation = rs.AggregationMode.OneSketchPerRecord
        order_config.ambiguous_policy = rs.AmbiguousPolicy.RejectRecord
        order_config.kmer_size = 5
        order_config.order_l = 2
        order_config.order_m = 32
        order_built = rs.build_fastx(str(path), [order_config])
        request = rs.QueryRequest()
        request.metric = rs.QueryMetric.OrderSimilarity
        order_result = order_built[0].sketch.query(
            order_built[0].sketch, request
        )
        assert order_result.order_similarity == 1.0

    print("test_fastx_python: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
