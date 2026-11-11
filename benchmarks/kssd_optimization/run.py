#!/usr/bin/env python3
"""Sequential paired benchmarks. Timed commands never overlap each other."""
from validate import REPO, ROOT, RUN, OLD, NEW, PREVIOUS, environment
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''): h.update(block)
    return h.hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def setup():
    panels = RUN / 'panels'; panels.mkdir(exist_ok=True)
    original = (PREVIOUS / '233_original.list').read_text().splitlines()
    subset = [original[i * (len(original)-1) // 1023] for i in range(1024)]
    (panels / 'real1024.list').write_text('\n'.join(subset) + '\n')
    real128 = (PREVIOUS / 'smoke.list').read_text().splitlines()
    (panels / 'real128.list').write_text('\n'.join(real128) + '\n')
    (panels / 'one.list').write_text(real128[0] + '\n')
    (panels / 'builder8.list').write_text('\n'.join(real128[:8]) + '\n')
    return panels


def execute(tag, argv, threads, directory, timed=True):
    directory.mkdir(parents=True, exist_ok=True)
    output = directory / (tag + '.stdout')
    error = directory / (tag + '.stderr')
    command = ['/usr/bin/time', '-v', '-o', str(directory / (tag + '.time')),
               'taskset', '--cpu-list', ','.join(map(str, range(threads))), *map(str, argv)]
    write(directory / (tag + '.command.json'), dict(argv=command, cwd=str(REPO/'examples'),
          environment={k:v for k,v in environment(threads).items() if k.startswith(('OMP_', 'RABBIT_'))}))
    write(RUN / 'progress.json', dict(status='running', phase=tag, runner_pid=os.getpid(),
          utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print('START', tag, flush=True)
    start = time.monotonic()
    with output.open('w') as stdout, error.open('w') as stderr:
        subprocess.run(command, cwd=REPO/'examples', env=environment(threads),
                       stdout=stdout, stderr=stderr, check=True)
    elapsed = time.monotonic() - start
    memory = re.search(r'Maximum resident set size \(kbytes\): (\d+)', (directory/(tag+'.time')).read_text())
    print('DONE', tag, f'{elapsed:.3f}s', flush=True)
    return elapsed, int(memory[1])


def micro(panels):
    directory = RUN / 'micro'; directory.mkdir(exist_ok=True)
    configurations = [
        ('real128', 'api', 1, 0), ('real128', 'api', 8, 0),
        ('one', 'api', 1, 1000), ('builder8', 'builder', 1, 0),
    ]
    rows = []; fingerprints = {}
    for panel, mode, threads, fragment in configurations:
        for version, binary in [('before', OLD.parent/'kssd_bench'), ('after', NEW.parent/'kssd_bench')]:
            tag = f'{panel}_{mode}_f{fragment}_t{threads}_{version}'
            dump = directory / (tag + '.bin')
            wall, rss = execute(tag, [binary, panels/(panel+'.list'), mode, threads, 3, fragment, dump, 1], threads, directory)
            values = [json.loads(line) for line in (directory/(tag+'.stdout')).read_text().splitlines()]
            digest = sha(dump); key = (panel, mode, fragment)
            if key in fingerprints and fingerprints[key] != digest: raise ValueError('micro output mismatch: '+tag)
            fingerprints[key] = digest
            rows.append(dict(tag=tag, panel=panel, mode=mode, version=version, threads=threads,
                             fragment=fragment, peak_rss_kib=rss, dump_sha256=digest, trials=values))
            write(RUN/'micro.json', rows)
    reference = fingerprints[('real128', 'api', 0)]
    for mode in ['bitmap', 'map', 'direct', 'bitmap-set', 'file']:
        tag = 'real128_' + mode
        dump = directory / (tag + '.bin')
        wall, rss = execute(tag, [NEW.parent/'kssd_bench', panels/'real128.list', mode, 1, 3, 0, dump, 1], 1, directory)
        values = [json.loads(line) for line in (directory/(tag+'.stdout')).read_text().splitlines()]
        digest = sha(dump)
        if digest != reference: raise ValueError('lookup/dedup/input changed hashes: '+tag)
        rows.append(dict(tag=tag, panel='real128', mode=mode, version='after', threads=1,
                         fragment=0, peak_rss_kib=rss, dump_sha256=digest, trials=values))
        write(RUN/'micro.json', rows)
    return rows


def cli(panels, full_only=False):
    directory = RUN / 'cli'; directory.mkdir(exist_ok=True)
    output_root = REPO.parent / '.rabbitsketch_benchmark_outputs/kssd_optimization_20261111'
    output_root.mkdir(exist_ok=True)
    fingerprint = {}; rows = []
    digest_binary = NEW.parent / 'output_digest'
    jobs = [('real1024', panels/'real1024.list', threads, repeat)
            for threads in [1, 8, 64] for repeat in range(3)]
    if full_only: jobs = []
    jobs += [('full233', PREVIOUS/'233_original.list', 64, repeat) for repeat in range(1)]
    previous = json.loads((PREVIOUS/'jobs/kssd_t64_r1/complete.json').read_text())['output_digest']
    for panel, listing, threads, repeat in jobs:
        order = [('before', OLD), ('after', NEW)]
        if (repeat % 2) == (0 if full_only else 1): order.reverse()
        for version, binary in order:
            tag = f'{panel}_t{threads}_r{repeat}_{version}'
            output = output_root / (tag+'.tsv')
            wall, rss = execute(tag, [binary, '--kssd', '--index', '-i', listing, '-l',
                                     '-t', threads, '-o', output], threads, directory)
            log = (directory/(tag+'.stderr')).read_text()
            def seconds(label):
                match = re.search(re.escape(label)+r': ([0-9.eE+-]+) s', log)
                return float(match[1]) if match else None
            digest = json.loads(subprocess.check_output([digest_binary, output], text=True))
            write(directory/(tag+'.digest.json'), digest)
            if panel in fingerprint and fingerprint[panel] != digest: raise ValueError('CLI output mismatch: '+tag)
            fingerprint[panel] = digest
            if panel == 'full233':
                for key in ['rows','canonical_xxh3_128_sum','canonical_xxh3_128_xor']:
                    if digest[key] != previous[key]: raise ValueError('full output differs from preserved 233 reference')
            # Retain a reviewable sample; remove only this task's verified TSV.
            with output.open() as stream:
                sample = ''.join(next(stream, '') for _ in range(101))
            (directory/(tag+'.sample.tsv')).write_text(sample)
            row = dict(tag=tag, panel=panel, threads=threads, repeat=repeat, version=version,
                       process_seconds=wall, peak_rss_kib=rss, sketch_seconds=seconds('sketch time'),
                       legacy_sketch_local_index_seconds=seconds('sketch + local index'),
                       index_seconds=seconds('build CSR index'), similarity_seconds=seconds('dist time'),
                       internal_total_seconds=seconds('total'), output_digest=digest,
                       binary_sha256=sha(binary), input_list_sha256=sha(listing))
            if row['internal_total_seconds'] is not None and row['similarity_seconds'] is not None:
                row['prepare_seconds'] = row['internal_total_seconds'] - row['similarity_seconds']
            rows.append(row); write(RUN/'cli.json', rows)
            output.unlink()
    return rows


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('phase', choices=['micro','cli','full'])
    args = parser.parse_args()
    lock = (ROOT/'runs/.active_benchmark.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    panels = setup()
    write(RUN/(args.phase+'_provenance.json'), dict(
        before_binary_sha256=sha(OLD), after_binary_sha256=sha(NEW),
        before_bench_sha256=sha(OLD.parent/'kssd_bench'), after_bench_sha256=sha(NEW.parent/'kssd_bench'),
        cpu_affinity='0..threads-1', distance_cutoff=.05, effective_k=20, drlevel=3,
        cache_policy='ordinary filesystem cache; no eviction; paired order alternates',
        timing='CLI wall includes read, sketch, index, comparison and output; validation excluded'))
    if args.phase == 'micro': micro(panels)
    else: cli(panels, full_only=args.phase == 'full')
    write(RUN/'progress.json', dict(status='completed', phase=args.phase))


if __name__ == '__main__': main()
