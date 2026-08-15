#!/usr/bin/env python3
"""End-to-end smoke test for the public Python FASTX build surface."""

import gzip
import pathlib
import sys
import tempfile
import threading
import time


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

    with tempfile.TemporaryDirectory(prefix="rabbitsketch-fastx-python-") as tmp:
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
