#!/usr/bin/env python3
"""Compare the BinDash --index entry point with independent pair-mode runs.

Usage: python tests/test_bindash_index_cli.py examples/rabbitsketch
Exact floating cutoff boundaries and native full-width values are additionally
covered by the C++ bindash_index_test.
"""
import gzip
import itertools
import pathlib
import random
import subprocess
import sys
import tempfile


def main():
    binary = pathlib.Path(sys.argv[1]).resolve()
    rng = random.Random(731)
    base = "".join(rng.choices("ACGT", k=20000))
    mutated = "".join(rng.choice("ACGT") if rng.random() < .03 else c for c in base)
    sequences = [base, base, mutated, "".join(rng.choices("ACGT", k=20000)),
                 "ACGT", "", "A" * 100]
    with tempfile.TemporaryDirectory(prefix="rabbitsketch_bindash_index_cli_") as tmp:
        directory = pathlib.Path(tmp)
        files = []
        for i, sequence in enumerate(sequences):
            path = directory / (f"{i}.fa.gz" if i == 2 else f"{i}.fa")
            if i == 2:
                with gzip.open(path, "wt") as stream:
                    stream.write(f">seq\n{sequence}\n")
            else:
                path.write_text(f">seq\n{sequence}\n")
            files.append(str(path))
        listing = directory / "files.list"
        listing.write_text("\n".join(files) + "\n")
        output = directory / "result.tsv"

        # Exercise 32-bit, 64-bit and composite keys, nondefault k/seed,
        # threads > input count, plus counters exceeding 16 bits.
        for bits, groups, kmer, seed, threads in [
            (16, 16, 21, 42, 4), (1, 1, 17, 7, 12), (32, 1, 21, 91, 2),
            (58, 1, 21, 42, 2), (64, 1, 23, 3, 4), (16, 1024, 21, 42, 2),
        ]:
            options = ["--bindash", "--bbits", str(bits), "--sketch64", str(groups),
                       "-k", str(kmer), "--seed", str(seed), "-t", str(threads)]
            oracle = []
            for left, right in itertools.combinations(files, 2):
                subprocess.run([str(binary), *options, "-i", left, right, "-o", str(output)],
                               check=True, capture_output=True, text=True)
                oracle.append(output.read_text().strip())
            for cutoff in [0.0, 0.05, 0.2, 1.0]:
                result = subprocess.run([str(binary), *options, "--index", "-l", "-i",
                    str(listing), "-d", str(cutoff), "-o", str(output)],
                    check=True, capture_output=True, text=True)
                assert "BinDash index:" in result.stderr and "build CSR index:" in result.stderr
                assert "falling back" not in result.stderr
                actual = sorted(output.read_text().splitlines())
                expected = sorted(row for row in oracle if float(row.split("\t")[-1]) < cutoff)
                assert actual == expected, (bits, groups, cutoff, set(actual) ^ set(expected))
            print(f"BinDash index CLI: bbits={bits}, bins={groups * 64}, threads={threads}: PASS")

        # The CLI must fail explicitly, never substitute another comparison.
        invalid = subprocess.run([str(binary), "--bindash", "--index", "-l", "-i",
            str(listing), "--bbits", "65", "-o", str(output)], capture_output=True, text=True)
        assert invalid.returncode != 0 and "invalid BinDash" in invalid.stderr
        single = directory / "single.list"
        single.write_text(files[0] + "\n")
        subprocess.run([str(binary), "--bindash", "--index", "-l", "-i", str(single),
            "-t", "16", "-o", str(output)], check=True, capture_output=True)
        assert output.read_text() == ""
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
