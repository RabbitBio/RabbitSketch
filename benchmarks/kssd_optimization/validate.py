#!/usr/bin/env python3
"""Compare the actual CLI with the preserved native KSSD estimator."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess

REPO = Path(__file__).resolve().parents[2]
ROOT = Path('/home/bigssd/zt/rabbitsketch2_benchmark')
RUN = ROOT / 'runs/kssd_optimization_20261111'
OLD = ROOT / 'build_kssd_optimization_reference/rabbitsketch_main'
NEW = ROOT / 'build_kssd_optimization_after/rabbitsketch_main'
PREVIOUS = ROOT / 'runs/thread_scaling_index_233_20261109'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def environment(threads, indexed=True):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_THREAD_LIMIT=str(threads),
               OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close', OMP_PLACES='cores',
               OMP_MAX_ACTIVE_LEVELS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', LC_ALL='C')
    for name in ['RABBIT_FORCE_INDEX', 'RABBIT_FORCE_ALLPAIRS', 'RABBITSKETCH_SIMD']:
        env.pop(name, None)
    env['RABBIT_FORCE_INDEX' if indexed else 'RABBIT_FORCE_ALLPAIRS'] = '1'
    return env


def pairs(path):
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith('genome0\tgenome1\t'): continue
        row = line.split('\t')
        if len(row) not in (3, 5): raise ValueError('bad columns')
        key = tuple(sorted(row[:2]))
        if key in result or key[0] == key[1]: raise ValueError('duplicate/self pair')
        result[key] = row[-1]
    return result


def main():
    directory = RUN / 'validation'
    directory.mkdir(exist_ok=True)
    rng = random.Random(110311)
    sequence = ''.join(rng.choices('ACGT', k=120000))
    paths = []
    for i in range(8):
        bases = list(sequence)
        if i >= 2:
            for pos in range(31, len(bases), 101 if i < 5 else 4096):
                bases[pos] = 'SDWsdwN'[pos % 7]
        if i == 1: bases = list(sequence.lower())
        records = [''.join(bases[:60001]), 'AC', ''.join(bases[60001:])]
        fastq = i in [4, 5, 6, 7]
        compressed = i in [2, 3, 6, 7]
        text = ''
        for j, record in enumerate(records):
            if fastq:
                text += f'@r{j}\n{record}\n+\n' + 'I' * len(record) + '\n'
            else:
                text += f'>r{j}\r\n' + '\r\n'.join(record[p:p+79] for p in range(0, len(record), 79)) + '\r\n'
        path = directory / f'input_{i}.fastx'
        if compressed:
            with gzip.open(path, 'wt') as stream: stream.write(text)
        else: path.write_text(text)
        paths.append(path)
    listing = directory / 'formats.list'
    listing.write_text(''.join(str(p) + '\n' for p in paths))
    panels = [
        ('real128', PREVIOUS / 'smoke.list', 3, .05),
        ('ambiguous', PREVIOUS / 'revisions/kssd_encoder_v2/validation/ambiguous.list', 3, .05),
        ('boundaries', PREVIOUS / 'boundaries/panel.list', 3, .05),
        *[(f'formats_dr{dr}', listing, dr, .05) for dr in [0, 1, 2, 3]],
        ('zero_cutoff', listing, 3, 0),
        ('wide_cutoff', listing, 3, .2),
    ]
    rows = []
    for name, panel, dr, cutoff in panels:
        reference = None
        for binary, threads, indexed in [(OLD, 1, False), (NEW, 1, True), (NEW, 8, True), (NEW, 1, False)]:
            tag = f'{name}_{"after" if binary == NEW else "before"}_t{threads}_{"index" if indexed else "native"}'
            output = directory / (tag + '.tsv')
            command = ['taskset', '--cpu-list', ','.join(map(str, range(threads))), str(binary),
                       '--kssd', *(['--index'] if indexed else []), '--drlevel', str(dr),
                       '-d', str(cutoff), '-i', str(panel), '-l', '-t', str(threads), '-o', str(output)]
            with (directory / (tag + '.log')).open('w') as log:
                subprocess.run(command, cwd=REPO / 'examples', env=environment(threads, indexed),
                               stdout=log, stderr=subprocess.STDOUT, check=True)
            observed = pairs(output)
            if reference is None: reference = observed
            if observed != reference:
                raise ValueError(f'{tag}: missing={len(reference.keys()-observed.keys())} '
                                 f'extra={len(observed.keys()-reference.keys())} '
                                 f'changed={sum(reference[k]!=observed[k] for k in reference.keys()&observed.keys())}')
            rows.append(dict(case=tag, pairs=len(observed), match=True, command=command, output_sha256=sha(output)))
        print(name, 'PASS', len(reference), 'pairs', flush=True)
    (RUN / 'validation.json').write_text(json.dumps(dict(all_passed=True, cases=rows,
        before_sha256=sha(OLD), after_sha256=sha(NEW)), indent=2) + '\n')


if __name__ == '__main__': main()
