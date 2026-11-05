#!/usr/bin/env python3
"""Sequential, resumable native execution with enforced freeze-before-validation."""
import argparse
import datetime
import fcntl
import importlib.util
import shutil

from scenario_common import *

spec = importlib.util.spec_from_file_location('oct30_execution', HERE.parent / 'comprehensive_accuracy/run.py')
old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
old.RUN = RUN; old.BINARY = BINARY


def jobs_for(split, pilot):
    jobs = []
    grid = [('minhash', CONFIG['pilot_kmv_size'], 6)] if pilot else [('exact', 1, 6)] + native_grid()
    for k in CONFIG['k']:
        for algorithm, parameter, half in grid:
            seeds = ([CONFIG['pilot_seed']] if pilot else [42] if algorithm in ['exact', 'kssd']
                     else CONFIG[split + '_seeds'])
            for seed in seeds:
                job = dict(split=split, purpose='pilot' if pilot else 'outcomes', algorithm=algorithm,
                    parameter=parameter, half_subk=half, k=k, seed=seed, mode='accuracy',
                    objects=str(DATA / (split + '_objects.tsv')), pairs=str(DATA / (split + '_pairs.tsv')))
                job['job_id'] = key(json.dumps(job, sort_keys=True))[:24]
                jobs.append(job)
    return jobs


def check_freeze(stage):
    if stage == 'development': return
    path = RUN / 'model_frozen.json'
    if not path.exists(): raise ValueError('development model must be frozen before validation pilot')
    freeze = json.loads(path.read_text())
    for filename, digest in freeze['files'].items():
        if sha(filename) != digest: raise ValueError('frozen model or fitting code changed: ' + filename)
    if stage == 'validation-outcomes':
        path = RUN / 'decisions_frozen.json'
        if not path.exists(): raise ValueError('freeze validation decisions before measuring outcomes')
        freeze = json.loads(path.read_text())
        for filename, digest in freeze['files'].items():
            if sha(filename) != digest: raise ValueError('frozen decisions changed: ' + filename)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage', choices=['development', 'validation-pilot', 'validation-outcomes'], required=True)
    ap.add_argument('--threads', type=int, default=64)
    args = ap.parse_args()
    RUN.mkdir(parents=True, exist_ok=True)
    lock = (ROOT / 'runs/.active_benchmark.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    check_freeze(args.stage)
    split = 'development' if args.stage == 'development' else 'validation'
    prep = json.loads((DATA / (split + '_preparation.json')).read_text())
    for filename, digest in prep['files'].items():
        if sha(filename) != digest: raise ValueError('prepared manifest changed')
    if args.stage == 'development': jobs = jobs_for(split, True) + jobs_for(split, False)
    else: jobs = jobs_for(split, args.stage == 'validation-pilot')
    plan_path = RUN / (args.stage + '_plan.json')
    identity = dict(jobs=jobs, binary_sha256=sha(BINARY), run_sha256=sha(__file__),
        helper_sha256=sha(HERE.parent / 'comprehensive_accuracy/run.py'),
        config_sha256=sha(HERE / 'config.json'), preparation_sha256=sha(DATA / (split + '_preparation.json')))
    if plan_path.exists() and json.loads(plan_path.read_text()) != identity:
        raise ValueError('native execution plan identity changed')
    atomic(plan_path, identity)
    snapshot = RUN / 'code_snapshot'; snapshot.mkdir(exist_ok=True)
    for filename in ['run.py', 'prepare.py', 'scenario_common.py', 'config.json', 'PROTOCOL.md']:
        dst = snapshot / filename
        if not dst.exists(): shutil.copyfile(HERE / filename, dst)
    pairs = {r['pair_id'] for r in read_tsv(DATA / (split + '_pairs.tsv'))}
    objects = {r['object_id'] for r in read_tsv(DATA / (split + '_objects.tsv'))}
    for i, job in enumerate(jobs, 1):
        destination = RUN / 'native' / job['job_id']; destination.mkdir(parents=True, exist_ok=True)
        manifest = destination / 'manifest.tsv'
        content = '#regimes-v1\n' + ''.join(f'{k}\t{v}\n' for k, v in job.items()
            if k not in ['split', 'purpose', 'job_id'])
        if not manifest.exists(): manifest.write_text(content)
        elif manifest.read_text() != content: raise ValueError('native manifest changed')
        print(f'{args.stage} {i}/{len(jobs)} {job["purpose"]} {job["algorithm"]} '
              f'p={job["parameter"]} k={job["k"]} seed={job["seed"]}', flush=True)
        old.execute(manifest, destination, args.threads, expected_pairs=pairs, expected_objects=objects)
        atomic(RUN / 'progress.json', dict(stage=args.stage, completed_jobs=i, planned_jobs=len(jobs), last_job=job))
    atomic(RUN / (args.stage + '_complete.json'), dict(plan_sha256=sha(plan_path), jobs=len(jobs),
        completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat()))


if __name__ == '__main__': main()
