#!/usr/bin/env python3
"""Check full-frequency defaults and explicit legacy presampling end to end.

Usage: python tests/test_probminhash_cli.py examples/rabbitsketch
The weighted two-element fixture has analytical probability Jaccard 0.2.
Its second canonical k-mer is rejected by the old 1/8 admission rule, which
turns that similarity into 1.0; this distinguishes both input semantics.
"""
import gzip
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def fasta(records):
    return "".join(f">r{i}\n{seq}\n" for i, seq in enumerate(records))


def fastq(records, wrapped=False):
    pieces = []
    for i, seq in enumerate(records):
        quality = "I" * len(seq)
        if wrapped:
            half = len(seq) // 2
            seq = seq[:half] + "\n" + seq[half:]
            quality = quality[:half] + "\n" + quality[half:]
        pieces.append(f"@r{i}\n{seq}\n+\n{quality}\n")
    return "".join(pieces)


def main():
    binary = Path(sys.argv[1]).resolve()
    env = dict(os.environ, OMP_DYNAMIC="FALSE", OMP_THREAD_LIMIT="4")
    with tempfile.TemporaryDirectory(prefix="rabbitsketch_pmh_cli_") as tmp:
        directory = Path(tmp)
        output = directory / "result.tsv"
        listing = directory / "inputs.list"
        k = 21
        # Canonical integer identities 0 and 1. The latter was prefiltered out.
        zero, one = "A" * k, "A" * (k - 1) + "C"
        records_a, records_b = [zero] * 9 + [one], [zero] + [one] * 9
        a, b, scaled = (directory / name for name in ("a.fa", "b.fa", "scaled.fa"))
        a.write_text(fasta(records_a))
        b.write_text(fasta(records_b))
        scaled.write_text(fasta(records_a * 2))

        def run(paths, divisor=None, mode="pair", threads=1, size=k):
            command = [str(binary), "--probminhash", "-m", "4096", "-k", str(size),
                       "-t", str(threads), "-d", "1", "-o", str(output)]
            if divisor is not None:
                command += ["--pmh-sample-divisor", str(divisor)]
            if mode == "pair":
                command += ["-i", *map(str, paths)]
            else:
                listing.write_text("\n".join(map(str, paths)) + "\n")
                command += ["-l", "-i", str(listing)]
                if mode == "index":
                    command.append("--index")
            result = subprocess.run(command, env=env, capture_output=True, text=True)
            assert result.returncode == 0, (command, result.stderr)
            expected_divisor = 1 if divisor is None else divisor
            assert f"sample divisor={expected_divisor}" in result.stderr, result.stderr
            rows = []
            for line in output.read_text().splitlines():
                left, right, distance = line.split("\t")[:3]
                rows.append((*sorted((left, right)), float(distance)))
            return sorted(rows)

        full = run([a, b])
        assert full == run([a, b], 1), "default must be full-frequency input"
        exp_distance = math.exp(-k * full[0][2])
        estimated_jp = exp_distance / (2 - exp_distance)
        assert abs(estimated_jp - 0.2) < 0.06, ("full frequency truth", estimated_jp)
        assert run([a, b], 8)[0][2] == 0, "explicit legacy sample should contain one identity"
        assert run([a, scaled])[0][2] == 0, "uniform abundance scaling must preserve J_P"

        # Plain/gzip and fast-parser/fallback paths must count identical records.
        variants = {
            "a.fa.gz": (fasta(records_a), True),
            "a.fq": (fastq(records_a), False),
            "a.fq.gz": (fastq(records_a), True),
            "wrapped.fq.gz": (fastq(records_a, wrapped=True), True),
            "mixed.fq.gz": (fastq(records_a[:1]) + fastq(records_a[1:], wrapped=True), True),
            "crlf.fq.gz": (fastq(records_a).replace("\n", "\r\n"), True),
            "no_final_newline.fq.gz": (fastq(records_a).rstrip("\n"), True),
            "case_ambiguous.fa": (fasta([s.lower() for s in records_a] + ["N" * k, "A" * (k - 1)]), False),
        }
        for name, (body, compressed) in variants.items():
            path = directory / name
            path.write_bytes(gzip.compress(body.encode()) if compressed else body.encode())
            for divisor in (1, 8):
                assert run([path, b], divisor)[0][2] == run([a, b], divisor)[0][2], name
        multi = directory / "split_member.fq.gz"
        body = fastq(records_a).encode()
        split = len(fastq(records_a[:1]).encode()) + 12
        multi.write_bytes(gzip.compress(body[:split]) + gzip.compress(body[split:]))
        assert run([multi, b])[0][2] == full[0][2], "gzip member boundary lost or duplicated counts"

        # A sole k-mer rejected by the previous admission filter must survive.
        single, reverse = directory / "single.fa", directory / "reverse.fa"
        single.write_text(fasta([one]))
        reverse.write_text(fasta(["G" + "T" * (k - 1)]))
        assert run([single, reverse])[0][2] == 0, "default dropped a valid canonical k-mer"
        single.write_text(fasta(["A" * 31 + "C"]))
        reverse.write_text(fasta(["G" + "T" * 31]))
        assert run([single, reverse], size=32)[0][2] == 0, "k=32 mask or canonicalization failed"

        for divisor in (None, 1, 8):
            paths = [a, b, scaled, directory / "a.fa.gz"]
            expected = run(paths, divisor, "list")
            assert len(expected) == 6, expected
            for mode in ("list", "index"):
                for threads in (1, 4):
                    assert run(paths, divisor, mode, threads) == expected, (divisor, mode, threads)

        for bad in ("0", "-1", "1.5", "eight", "4294967296", "99999999999999999999999"):
            result = subprocess.run([str(binary), "--probminhash", "-i", str(a), str(b),
                                     "--pmh-sample-divisor", bad], env=env, capture_output=True, text=True)
            assert result.returncode != 0 and "--pmh-sample-divisor" in result.stderr, bad
        result = subprocess.run([str(binary), "--bindash", "-i", str(a), str(b),
                                 "--pmh-sample-divisor", "8"], env=env, capture_output=True, text=True)
        assert result.returncode != 0 and "requires --probminhash" in result.stderr
    print("ProbMinHash CLI: full frequencies, explicit sampling, FASTX formats, pair/list/index, k=32: PASS")


if __name__ == "__main__":
    main()
