#!/usr/bin/env python3
"""Public cardinality bindings, including unsupported/saturated states."""
import math
import random
import sys


def rejects(exception, function):
    try:
        function()
    except exception:
        return
    raise AssertionError(f"expected {exception.__name__}")


def main():
    sys.path.insert(0, sys.argv[1])
    import rabbitsketch as rs

    rng = random.Random(20261117)
    sequence = "".join(rng.choices("ACGT", k=5000))
    factories = [lambda: rs.FastKMV(256), lambda: rs.FracMinHash(4),
                 lambda: rs.HyperLogLog(10, 21, 42),
                 lambda: rs.SetSketch(10, 2.0, 5.0, 21, 42),
                 lambda: rs.BinDash(128),
                 lambda: rs.Kssd(rs.kssd_parameter_t(10, 5, 2))]
    for factory in factories:
        assert factory().cardinality() == 0
        sketch = factory()
        sketch.update(sequence)
        sketch.finalize()
        estimate = sketch.cardinality()
        assert math.isfinite(estimate) and estimate > 0
        assert sketch.merge(sketch).cardinality() == estimate
        sketch.finalize()
        assert sketch.cardinality() == estimate

    saturated = rs.BinDash(1)
    saturated.update(sequence)
    saturated.finalize()
    assert saturated.cardinality() == math.inf

    assert not hasattr(rs.ProbMinHash(), "cardinality")
    assert not rs.capabilities_for(rs.Algorithm.ProbMinHash).supports(rs.Capability.Cardinality)
    config = rs.SketchConfig()
    config.algorithm = rs.Algorithm.ProbMinHash
    config.sampling = rs.SamplingMode.AlgorithmNative
    config.resolution = 128
    record = rs.FastxRecord()
    record.name = "sample"
    record.sequence = sequence
    builder = rs.MultiSketchBuilder([config])
    builder.update(record)
    rejects(RuntimeError, builder.finalize()[0].sketch.cardinality)
    print("Python cardinality bindings passed")


if __name__ == "__main__":
    main()
