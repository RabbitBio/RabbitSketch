#!/usr/bin/env python3
"""Run one existing-main process at a time, with OpenMP over a sample manifest."""
import argparse
from collections import OrderedDict
import datetime as dt
import fcntl
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import threading
import time

from run import (DEFAULT_ROOT, HERE, atomic_json, canonical, digest, file_sha,
                 make_plan, source_identity)

PROFILE = 'main_openmp_batch'
FIELDS = ('job_id', 'algorithm', 'parameter', 'k', 'n_a', 'n_b', 'intersection',
          'data_seed', 'hash_seed', 'query_min_ms', 'half_subk', 'bbits')


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def group_jobs(jobs):
    groups = OrderedDict()
    for job in jobs:
        groups.setdefault((job['config_id'], job['k']), []).append(job)
    return list(groups.values())


def write_manifest(path, jobs):
    Path(path).write_text('#' + '\t'.join(FIELDS) + '\n' + ''.join(
        '\t'.join(str(job[field]) for field in FIELDS) + '\n' for job in jobs))


def combine(job, measured, threads, batch_id):
    result = dict(job, **{k: v for k, v in measured.items() if k not in job},
                  measured_ground_truth=measured.get('ground_truth'), batch_id=batch_id)
    if result.get('status') == 'ok':
        mismatches = [k for k in job if k in measured and measured[k] != job[k]]
        if (mismatches or measured.get('measurement_profile') != PROFILE or
                measured.get('batch_threads') != threads):
            result.update(status='error', reason='worker metadata/truth differs from registered job')
        elif not math.isclose(job['ground_truth'], measured.get('ground_truth', -1), abs_tol=1e-15):
            result.update(status='error', reason='worker/plan truth mismatch')
    return result


def export_rows(output, jobs):
    temporary = output / 'raw.jsonl.tmp'
    with temporary.open('w') as stream:
        for job in jobs:
            path = output / 'rows' / (job['job_id'] + '.json')
            if path.exists():
                stream.write(canonical(json.loads(path.read_text())) + '\n')
    temporary.replace(output / 'raw.jsonl')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=HERE / 'configs/pilot_six.json')
    p.add_argument('--binary', type=Path, default=DEFAULT_ROOT / 'build_validated/rabbitsketch_pilot')
    p.add_argument('--source', type=Path, default=DEFAULT_ROOT / 'source/RabbitSketch-validated')
    p.add_argument('--output', type=Path, default=DEFAULT_ROOT / 'runs/pilot_six')
    p.add_argument('--threads', type=int, default=len(os.sched_getaffinity(0)))
    p.add_argument('--limit', type=int, help='maximum unfinished pairs to run, then resume later')
    p.add_argument('--retry-failed', action='store_true')
    p.add_argument('--plan-only', action='store_true')
    args = p.parse_args(argv)
    affinity = sorted(os.sched_getaffinity(0))
    if not 1 <= args.threads <= len(affinity):
        p.error('--threads must be between 1 and the number of available CPUs')
    if args.limit is not None and args.limit < 1:
        p.error('--limit must be positive')
    plan = make_plan(json.loads(args.config.read_text()))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    lock = (output / '.runner.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    plan_path = output / 'plan.json'
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        raise ValueError('output contains a different plan; use a new directory')
    atomic_json(plan_path, plan)
    groups = group_jobs(plan['jobs'])
    print(json.dumps(dict(jobs=len(plan['jobs']), batches=len(groups), threads=args.threads,
                          feasible_regimes=len(plan['regimes']),
                          skipped_regimes=len(plan['skipped_regimes']),
                          execution=PROFILE)), flush=True)
    if args.plan_only:
        return 0

    binary, source = args.binary.resolve(), args.source.resolve()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC='FALSE',
               OMP_THREAD_LIMIT=str(args.threads), OMP_PROC_BIND='spread',
               OMP_PLACES='threads', OPENBLAS_NUM_THREADS='1')
    cli_files = ('examples/main.cpp', 'examples/sub_command.cpp', 'examples/sub_command.h')
    identity = dict(binary_sha256=file_sha(binary), source=source_identity(source),
                    cli_sha256=digest([(name, file_sha(source / name)) for name in cli_files]),
                    runner_sha256=file_sha(__file__), helpers_sha256=file_sha(HERE / 'run.py'),
                    worker_source_sha256=file_sha(HERE / 'synthetic_pilot.cpp'),
                    adapter_source_sha256=file_sha(HERE / 'main_batch.cpp'),
                    machine=platform.platform(), cpu_affinity=affinity,
                    cpu_model=next((s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines()
                                    if s.startswith('model name')), 'unknown'),
                    simd_override=env.get('RABBITSKETCH_SIMD', ''),
                    measurement_profile=PROFILE, threads=args.threads,
                    openmp={k: env[k] for k in ('OMP_NUM_THREADS', 'OMP_DYNAMIC', 'OMP_THREAD_LIMIT',
                                               'OMP_PROC_BIND', 'OMP_PLACES')})
    manifest_path = output / 'manifest.json'
    if manifest_path.exists():
        if json.loads(manifest_path.read_text())['identity'] != identity:
            raise ValueError('source, runner, binary, execution profile or host changed; use a new directory')
    else:
        atomic_json(manifest_path, dict(identity=identity, config_sha256=plan['config_sha256'],
                    created_utc=utc(), rss_scope='shared_batch_process',
                    timeout_policy='single-pair limit * ceil(remaining batch pairs / threads)',
                    timing_scope='pair elapsed time under concurrent sample-list load; setup shared per batch'))

    pending, completed, failed = [], 0, 0
    for job in plan['jobs']:
        path = output / 'rows' / (job['job_id'] + '.json')
        if path.exists():
            previous = json.loads(path.read_text())
            if any(previous.get(k) != v for k, v in job.items()):
                raise ValueError(f'corrupt registered row: {path}')
            if previous['status'] == 'ok' or not args.retry_failed:
                completed += 1
                failed += previous['status'] != 'ok'
                continue
        pending.append(job)
    if args.limit is not None:
        pending = pending[:args.limit]
    batch_dir = output / 'batches'
    batch_dir.mkdir(exist_ok=True)
    executed = 0
    progress = dict(runner_pid=os.getpid(), planned_jobs=len(plan['jobs']),
                    measurement_profile=PROFILE, threads=args.threads,
                    completed=completed, failed=failed, phase='running')

    def save_progress(**changes):
        progress.update(changes, updated_utc=utc())
        atomic_json(output / 'progress.json', progress)

    save_progress()
    for batch_number, jobs in enumerate(group_jobs(pending), 1):
        # Repeated launches get distinct attempt files; completed row IDs remain stable.
        batch_id = digest([identity, [j['job_id'] for j in jobs], time.time_ns()])[:24]
        tsv = batch_dir / (batch_id + '.tsv')
        write_manifest(tsv, jobs)
        expected = {j['job_id']: j for j in jobs}
        received = set()
        cmd = [str(binary), '--screening-manifest', str(tsv), '-t', str(args.threads)]
        record = dict(batch_id=batch_id, command=cmd, started_utc=utc(),
                      planned_pairs=len(jobs), algorithm=jobs[0]['algorithm'],
                      parameter=jobs[0]['parameter'], k=jobs[0]['k'], phase='running')
        timeout = plan['config']['timeout_seconds'] * math.ceil(len(jobs) / args.threads)
        record['timeout_seconds'] = timeout
        start = time.monotonic()
        timed_out = threading.Event()
        print(f"Batch {batch_number}: {jobs[0]['algorithm']} p={jobs[0]['parameter']} "
              f"k={jobs[0]['k']}, {len(jobs)} pairs, -t {args.threads}", flush=True)
        with (batch_dir / (batch_id + '.stderr')).open('w', buffering=1) as stderr:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr,
                                    text=True, bufsize=1, env=env)
            record['pid'] = proc.pid
            atomic_json(batch_dir / (batch_id + '.json'), record)
            save_progress(batch_id=batch_id, batch_pid=proc.pid, batch_algorithm=jobs[0]['algorithm'],
                          batch_parameter=jobs[0]['parameter'], batch_k=jobs[0]['k'])

            def expire():
                if proc.poll() is None:
                    timed_out.set()
                    proc.kill()

            timer = threading.Timer(timeout, expire)
            timer.daemon = True
            timer.start()
            try:
                for line in proc.stdout:
                    if not line.startswith('{'):
                        stderr.write(line)
                        continue
                    measured = json.loads(line)
                    if 'batch_event' in measured:
                        if measured.get('threads') != args.threads or measured.get('pairs') != len(jobs):
                            raise ValueError('main did not start the requested OpenMP batch')
                        record.update(actual_threads=measured['threads'],
                                      shared_setup_seconds=measured['shared_setup_seconds'])
                        atomic_json(batch_dir / (batch_id + '.json'), record)
                        continue
                    job_id = measured.get('job_id')
                    if job_id not in expected or job_id in received:
                        raise ValueError('main emitted an unknown or duplicate job ID')
                    result = combine(expected[job_id], measured, args.threads, batch_id)
                    atomic_json(output / 'rows' / (job_id + '.json'), result)
                    received.add(job_id)
                    completed += 1
                    executed += 1
                    failed += result['status'] != 'ok'
                    save_progress(completed=completed, failed=failed)
                    print(f"[{completed}/{len(plan['jobs'])}] {result['algorithm']} p={result['parameter']} "
                          f"k={result['k']} rep={result['replicate']} {result['status']}", flush=True)
                exit_code = proc.wait()
            finally:
                timer.cancel()
                if proc.poll() is None:
                    proc.kill()
                    proc.wait()
                proc.stdout.close()
        for job_id in expected.keys() - received:
            result = dict(expected[job_id], status='timeout' if timed_out.is_set() else 'error',
                          reason=f'batch ended without a result (exit {exit_code})',
                          measurement_profile=PROFILE, batch_threads=args.threads,
                          rss_scope='shared_batch_process', batch_id=batch_id)
            atomic_json(output / 'rows' / (job_id + '.json'), result)
            completed += 1
            executed += 1
            failed += 1
        wall_seconds = time.monotonic() - start
        record.update(phase='finished', finished_utc=utc(), exit_code=exit_code,
                      timed_out=timed_out.is_set(), received_pairs=len(received),
                      missing_results=len(jobs)-len(received), wall_seconds=wall_seconds,
                      pairs_per_wall_second=len(received)/wall_seconds)
        atomic_json(batch_dir / (batch_id + '.json'), record)
        save_progress(completed=completed, failed=failed, batch_pid=None)
    export_rows(output, plan['jobs'])
    save_progress(phase='measured' if completed == len(plan['jobs']) else 'partial', batch_pid=None)
    print(f'Executed {executed}; observed {completed}/{len(plan["jobs"])}; failed {failed}.', flush=True)
    return 1 if failed else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
