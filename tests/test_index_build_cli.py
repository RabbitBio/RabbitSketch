#!/usr/bin/env python3
"""Check native-key indexed results against direct CLI comparisons.

Usage: python tests/test_index_build_cli.py examples/rabbitsketch [previous_binary]
The optional previous binary additionally checks construction-only equivalence.
"""
import gzip
import os
import pathlib
import random
import subprocess
import sys
import tempfile


def rows(path):
    result = []
    for line in path.read_text().splitlines():
        left, right, distance = line.split("\t")
        result.append((*sorted((left, right)), distance))
    return sorted(result)


def main():
    binary = pathlib.Path(sys.argv[1]).resolve()
    previous = pathlib.Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else None
    rng = random.Random(20261112)
    base = "".join(rng.choices("ACGT", k=20000))
    mutated = "".join(rng.choice("ACGT") if rng.random() < .02 else b for b in base)
    sequences = [base, base, base + base, mutated, base[:200], base[:21],
                 "A" * 100, "".join(rng.choices("ACGT", k=20000)), "ACGT", "", ""]
    with tempfile.TemporaryDirectory(prefix="rabbitsketch_index_build_cli_") as tmp:
        directory = pathlib.Path(tmp)
        files = []
        for i, sequence in enumerate(sequences):
            path = directory / (f"{i}.fa.gz" if i == 3 else f"{i}.fa")
            text = f">seq\n{sequence}\n"
            if i == 3:
                with gzip.open(path, "wt") as stream:
                    stream.write(text)
            else:
                path.write_text(text)
            files.append(str(path))
        listing = directory / "files.list"
        # The existing weighted direct comparator requires populated sketches
        # with compatible weight metadata. Empty key arrays are independently
        # checked by inverted_index_test. Use plain FASTA for its input parser.
        weighted_mutation = directory / "weighted_mutation.fa"
        weighted_mutation.write_text(f">seq\n{mutated}\n")
        weighted_files = [files[0], files[1], files[2], str(weighted_mutation), files[7]]
        output = directory / "result.tsv"

        def run(executable, options, cutoff, threads, indexed, limited=False):
            env = dict(os.environ, OMP_DYNAMIC="FALSE")
            if limited:
                env["OMP_THREAD_LIMIT"] = "1"
            command = [str(executable), *options, "-l", "-i", str(listing),
                       "-d", str(cutoff), "-t", str(threads), "-o", str(output)]
            if indexed:
                command.append("--index")
            result = subprocess.run(command, env=env, capture_output=True, text=True)
            assert result.returncode == 0, (command, result.stderr)
            if indexed:
                assert "build CSR index:" in result.stderr
            return rows(output)

        for options in [
            ["--fastkmv", "-K", "2"], ["--fastkmv", "-K", "1024"],
            ["--fastkmv", "-K", "64", "-k", "17", "--seed", "7"],
            ["--probminhash", "-m", "64"], ["--probminhash", "-m", "1024"],
            ["--fracminhash", "--scaled", "1"], ["--fracminhash", "--scaled", "1000"],
        ]:
            selected = weighted_files if options[0] == "--probminhash" else files
            listing.write_text("\n".join(selected) + "\n")
            for cutoff in [0.0, 0.02, 0.05, 0.25]:
                expected = run(binary, options, cutoff, 1, False)
                if previous:
                    assert run(previous, options, cutoff, 4, True) == expected, (
                        "pre-existing index/direct mismatch", options, cutoff)
                for threads, limited in [(1, False), (16, False), (8, True)]:
                    actual = run(binary, options, cutoff, threads, True, limited)
                    assert actual == expected, (options, cutoff, threads, limited,
                                                set(actual) ^ set(expected))
            print("Index construction CLI:", " ".join(options), "PASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
