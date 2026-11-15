#!/usr/bin/env python3
"""Exercise the public native and unified operation APIs against direct builds."""
import pathlib
import random
import sys
import tempfile


def rejects(exception, function):
    try:
        function()
    except exception:
        return
    raise AssertionError(f"expected {exception.__name__}")


def main():
    sys.path.insert(0, sys.argv[1])
    import rabbitsketch as rs

    rng = random.Random(20261115)
    a = "".join(rng.choices("ACGT", k=12000))
    b = a[:8000] + "".join(rng.choices("ACGT", k=4000))
    params = rs.kssd_parameter_t(10, 5, 1)
    specifications = [
        (rs.Algorithm.FastKMV, lambda: rs.FastKMV(128, 15), 64),
        (rs.Algorithm.FracMinHash, lambda: rs.FracMinHash(4, 15), 8),
        (rs.Algorithm.HyperLogLog, lambda: rs.HyperLogLog(8, 15, 42), 6),
        (rs.Algorithm.SetSketch, lambda: rs.SetSketch(8, 2.0, 5.0, 15, 42), 6),
        (rs.Algorithm.ProbMinHash, lambda: rs.ProbMinHash(128, 15), 64),
        (rs.Algorithm.BinDash, lambda: rs.BinDash(4, 15), 128),
        (rs.Algorithm.KSSD, lambda: rs.Kssd(params), 2),
    ]
    for algorithm, factory, target in specifications:
        left, right, union = factory(), factory(), factory()
        left.update(a)
        right.update(b)
        union.update(a)
        union.update(b)
        for sketch in (left, right, union):
            sketch.finalize()
            sketch.finalize()
            assert sketch.sealed
            rejects(RuntimeError, lambda: sketch.update(a))
        if algorithm in (rs.Algorithm.HyperLogLog, rs.Algorithm.SetSketch):
            rejects(RuntimeError, lambda: left.add_fingerprint(123))
        for method in ("finalize", "merge", "project", "store_hashes", "index_dict"):
            assert hasattr(left, method), (algorithm, method)
        assert hasattr(left, "cardinality") == (
            algorithm in (rs.Algorithm.HyperLogLog, rs.Algorithm.SetSketch)
        )
        assert left.merge(right).store_hashes() == union.store_hashes()
        projected = left.project(target)
        assert projected.merge(right.project(target)).store_hashes() == union.project(target).store_hashes()
        for index in (rs.index_dict([left, right, left], 3), type(left).index_dict([left, right, left], 2)):
            assert len(index) == 3 and index.key_count > 0
            pairs = index.pairs(min_similarity=0.1, num_threads=3)
            assert any((m.left, m.right) == (0, 2) for m in pairs)
            for match in pairs:
                source = (left, right, left)
                assert match.estimate.value == source[match.left].jaccard(source[match.right])
        assert len(rs.index_dict([projected, projected]).pairs(0.5)) == 1

    with tempfile.TemporaryDirectory(prefix="rabbitsketch-operations-") as temporary:
        path = pathlib.Path(temporary) / "input.fa"
        path.write_text(f">a\n{a}\n>b\n{b}\n", encoding="ascii")
        for algorithm, _, target in specifications:
            config = rs.SketchConfig()
            config.algorithm = algorithm
            config.kmer_size = 15
            config.resolution = 128
            config.aggregation = rs.AggregationMode.OneSketchPerRecord
            if algorithm == rs.Algorithm.FracMinHash:
                config.sampling = rs.SamplingMode.Scaled
                config.scaled = 4
            elif algorithm in (rs.Algorithm.HyperLogLog, rs.Algorithm.SetSketch):
                config.sampling = rs.SamplingMode.RegisterPrecision
                config.resolution = 8
            elif algorithm in (rs.Algorithm.BinDash, rs.Algorithm.ProbMinHash):
                config.sampling = rs.SamplingMode.AlgorithmNative
                if algorithm == rs.Algorithm.BinDash:
                    config.resolution = 256
            elif algorithm == rs.Algorithm.KSSD:
                config.sampling = rs.SamplingMode.KssdReduction
                config.kmer_size = 20
                config.kssd_half_subk = 5
                config.resolution = 1
            built = rs.build_fastx(str(path), [config])
            left, right = [result.sketch for result in built]
            config.aggregation = rs.AggregationMode.OneSketchPerFile
            direct = rs.build_fastx(str(path), [config])[0].sketch
            assert left.merge(right).store_hashes() == direct.store_hashes()
            assert left.project(target).merge(right.project(target)).store_hashes() == direct.project(target).store_hashes()
            assert len(rs.index_dict([left, left]).pairs(0.5)) == 1
            assert len(rs.BuiltSketch.index_dict([left, left]).search(0, 0.5)) == 1
            caps = rs.capabilities_for(algorithm)
            for capability in (rs.Capability.Merge, rs.Capability.Project,
                               rs.Capability.StoreHashes, rs.Capability.IndexDict):
                assert caps.supports(capability)
            cardinality = algorithm in (rs.Algorithm.HyperLogLog, rs.Algorithm.SetSketch)
            assert caps.supports(rs.Capability.Cardinality) == cardinality
            if cardinality:
                assert left.cardinality() > 1000
            else:
                rejects(RuntimeError, left.cardinality)

    # Smaller PMH views must not silently compare with independently generated races.
    parent = rs.ProbMinHash(128, 15, 42, 8)
    parent.update(a)
    projected = parent.project(32)
    independent = rs.ProbMinHash(32, 15, 42, 8)
    independent.update(a)
    assert projected.source_register_count == 128
    rejects(ValueError, lambda: projected.merge(independent))
    rejects(RuntimeError, lambda: projected.update(b))
    assert len(rs.index_dict([])) == 0
    rejects(ValueError, lambda: rs.index_dict([parent], 0))
    rejects(TypeError, lambda: rs.index_dict([parent, rs.FastKMV()]))

    config = rs.SketchConfig()
    record = rs.FastxRecord()
    record.name = "record"
    record.sequence = a
    modern = rs.MultiSketchBuilder([config])
    legacy = rs.MultiSketchBuilder([config])
    modern.update(record)
    legacy.update(record)
    finalized = modern.finalize()
    finished = legacy.finish("sketch")
    assert modern.finished and legacy.finished
    assert finalized[0].stats.sealed
    assert finalized[0].label == finished[0].label
    assert finalized[0].sketch.store_hashes() == finished[0].sketch.store_hashes()
    rejects(RuntimeError, lambda: modern.update(record))
    rejects(RuntimeError, lambda: modern.finish("again"))
    rejects(RuntimeError, legacy.finalize)
    print("Python native and unified operations passed for all seven algorithms")


if __name__ == "__main__":
    main()
