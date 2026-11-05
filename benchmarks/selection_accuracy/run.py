#!/usr/bin/env python3
"""Run the registered accuracy experiment with the existing native adapter."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
SCREENING = HERE.parent / "screening"
ROOT = Path("/home/bigssd/zt/rabbitsketch2_benchmark")
sys.path.insert(0, str(SCREENING))
import run_main  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "runs/selection_accuracy_20261012")
    parser.add_argument("--binary", type=Path, default=ROOT / "build_validated/rabbitsketch_pilot")
    parser.add_argument("--source", type=Path, default=ROOT / "source/RabbitSketch-validated")
    parser.add_argument("--threads", type=int, default=64)
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    lock = (ROOT / "runs/.active_benchmark.lock").open("a")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    files = [HERE / "config.json", HERE / "PROTOCOL.md", Path(__file__),
             SCREENING / "run.py", SCREENING / "run_main.py",
             SCREENING / "synthetic_pilot.cpp", SCREENING / "main_batch.cpp",
             SCREENING / "patches/setsketch_partitioned_joint.patch"]
    snapshot = output / "code_snapshot"
    hashes = {}
    for source in files:
        relative = source.resolve().relative_to(HERE.parent)
        target = snapshot / relative
        content = source.read_bytes()
        if target.exists() and target.read_bytes() != content:
            raise ValueError(f"registered source changed: {relative}; use a new output")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
        hashes[str(relative)] = hashlib.sha256(content).hexdigest()
    (snapshot / "sha256.json").write_text(json.dumps(hashes, indent=2, sort_keys=True) + "\n")
    argv = ["--config", str(HERE / "config.json"), "--output", str(output),
            "--binary", str(args.binary), "--source", str(args.source),
            "--threads", str(args.threads)]
    if args.plan_only:
        argv.append("--plan-only")
    return run_main.main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
