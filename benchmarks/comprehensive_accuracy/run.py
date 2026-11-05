#!/usr/bin/env python3
"""Registered, resumable accuracy suite; a single native process at a time."""
import argparse
import collections
import csv
import fcntl
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = Path('/home/bigssd/zt/rabbitsketch2_benchmark')
RUN = ROOT / 'runs/comprehensive_accuracy_20261030'
BUILD = ROOT / 'build_comprehensive_20261030'
BINARY = BUILD / 'cmake/rabbitsketch_accuracy'
CONFIG = json.loads((HERE / 'config.json').read_text())


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def atomic(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temp.replace(path)


def synthetic_units():
    rows = []
    for k, nb, ratio, containment in itertools.product(CONFIG['k'], CONFIG['larger_cardinalities'],
            CONFIG['size_ratios'], CONFIG['small_set_containments']):
        na = nb // ratio
        common = int(containment * na + .5)
        base = dict(k=k, n_a=na, n_b=nb, size_ratio=ratio, nominal_containment=containment,
                    intersection=common, truth_j=common / (na + nb - common), truth_c=common / na)
        regime_id = digest(base)[:16]
        for rep in range(CONFIG['replicates']):
            seed = int(digest([CONFIG['data_seed'], regime_id, rep])[:16], 16)
            row = dict(base, regime_id=regime_id, replicate=rep, data_seed=seed,
                       hash_seed=CONFIG['hash_seeds'][rep % len(CONFIG['hash_seeds'])],
                       split='calibration' if rep < CONFIG['replicates'] // 2 else 'validation')
            row['unit_id'] = digest(row)[:24]; rows.append(row)
    return pd.DataFrame(rows)


def real_inventory():
    data = RUN / 'data'; data.mkdir(exist_ok=True)
    old = ROOT / 'data/robustness_20261001'
    taxonomy = pd.read_csv(ROOT / 'analysis/robustness_20261001/natural_taxonomy.csv')
    genera = sorted(taxonomy.genus.unique(), key=lambda x: digest([20261030, int(x)]))
    split = {g: 'calibration' if i < len(genera) // 2 else 'validation' for i, g in enumerate(genera)}
    taxonomy['split'] = taxonomy.genus.map(split)
    taxonomy.to_csv(data / 'natural_taxonomy.csv', index=False)
    old_pairs = pd.read_csv(old / 'natural_pairs.tsv', sep='\t')
    old_lookup = {frozenset([r.a, r.b]): r.pair_id for r in old_pairs.itertuples()}
    indexed = taxonomy.set_index('object_id')
    pairs = []
    within = [pair for _, group in taxonomy.groupby('genus')
              for pair in itertools.combinations(sorted(group.object_id), 2)]
    controls = [(r.a, r.b) for r in old_pairs.itertuples() if r.genus_a != r.genus_b]
    for a, b in within + controls:
        ta, tb = indexed.loc[a], indexed.loc[b]
        pair_id = old_lookup.get(frozenset([a, b]), 'expanded_' + digest([a, b])[:18])
        pairs.append(dict(pair_id=pair_id, a=a, b=b, genus_a=int(ta.genus), genus_b=int(tb.genus),
                          species_a=int(ta.species), species_b=int(tb.species),
                          cluster=str(int(ta.genus)) if ta.genus == tb.genus else 'cross_genus',
                          split=ta['split'] if ta['split'] == tb['split'] else 'cross_split'))
    pair_frame = pd.DataFrame(pairs)
    if pair_frame.duplicated(['a', 'b']).any() or pair_frame.pair_id.duplicated().any():
        raise ValueError('duplicate real pair')
    pair_frame.to_csv(data / 'natural_pairs.tsv', sep='\t', index=False)
    shutil.copyfile(old / 'natural_objects.tsv', data / 'natural_objects.tsv')
    inventories = {'natural': (data / 'natural_objects.tsv', data / 'natural_pairs.tsv')}
    for phase, source, objects, pairs in [
        ('mixed', old, 'mixed_objects', 'mixed_references'),
        ('containment', ROOT / 'data/containment_holdout_20261001', 'contain_objects', 'contain_pairs'),
        ('weighted', ROOT / 'data/regimes_20261001', 'weighted_objects', 'weighted_pairs')]:
        op, pp = data / (phase + '_objects.tsv'), data / (phase + '_pairs.tsv')
        shutil.copyfile(source / (objects + '.tsv'), op)
        shutil.copyfile(source / (pairs + '.tsv'), pp)
        inventories[phase] = (op, pp)
    return inventories


def real_jobs(inventories):
    jobs = []
    configs = [('minhash', b // 16, 6) for b in CONFIG['pair_budgets']]
    configs += [('probminhash_unweighted', b // 16, 6) for b in CONFIG['pair_budgets']]
    configs += [('bindash', b // 4, 6) for b in CONFIG['pair_budgets']]
    configs += [(a, (b // 2).bit_length() - 1, 6) for a in ['hll', 'setsketch'] for b in CONFIG['pair_budgets']]
    configs += [('fracminhash', p, 6) for p in CONFIG['real_frac_scaled']]
    configs += [('kssd', p, half) for p, half in CONFIG['kssd_resolutions']]
    for phase, (objects, pairs) in inventories.items():
        ks = CONFIG['real_k'] if phase == 'natural' else [24]
        if phase == 'weighted':
            grid = [('probminhash', p, 6) for p in [32, 128, 512, 2048, 8192]]
        else:
            grid = [('exact', 1, 6)] + configs
        for k in ks:
            for alg, param, half in grid:
                seeds = [42] if alg in ['kssd', 'exact'] else (CONFIG['hash_seeds'] if phase == 'weighted' else CONFIG['real_seeds'])
                for seed in seeds:
                    row = dict(phase=phase, algorithm=alg, parameter=param, k=k, seed=seed,
                               half_subk=half, objects=str(objects), pairs=str(pairs),
                               mode='weighted' if phase == 'weighted' else 'accuracy')
                    row['job_id'] = digest(row)[:24]; jobs.append(row)
    return jobs


def write_units(path, units, verify=False):
    fields = ['unit_id', 'k', 'n_a', 'n_b', 'intersection', 'replicate', 'data_seed', 'hash_seed']
    with path.open('w') as stream:
        stream.write('#comprehensive-verify-v1\n' if verify else '#comprehensive-v1\n')
        units[fields].to_csv(stream, sep='\t', index=False, header=False)


def execute(manifest, destination, threads, expected_units=None, expected_pairs=None, expected_objects=None):
    destination.mkdir(parents=True, exist_ok=True)
    completed = destination / 'completion.json'
    if completed.exists():
        previous = json.loads(completed.read_text())
        if previous['manifest_sha256'] != sha(manifest) or previous['binary_sha256'] != sha(BINARY):
            raise ValueError('completed job identity changed')
        if previous['raw_sha256'] != sha(destination / 'raw.jsonl.gz'):
            raise ValueError('completed raw archive digest changed')
        return previous
    env = os.environ.copy(); env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE',
        OMP_PROC_BIND='spread', OMP_PLACES='threads', OPENBLAS_NUM_THREADS='1')
    command = [str(BINARY), '--screening-manifest', str(manifest), '-t', str(threads)]
    counts, units, pair_ids, object_ids = collections.Counter(), collections.Counter(), set(), set()
    configs = collections.defaultdict(set)
    last = 0.; invalid = 0
    with (destination / 'stderr.log').open('w') as err, gzip.open(destination / 'raw.jsonl.gz.tmp', 'wt', compresslevel=3) as raw:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err, text=True, env=env)
        try:
            for line in proc.stdout:
                if not line.startswith('{'):
                    err.write(line); continue
                row = json.loads(line)
                if row.get('status') == 'error':
                    raise ValueError(f'native exception: {row}')
                if expected_units is not None:
                    unit = row['unit_id']
                    if unit not in expected_units: raise ValueError('unregistered synthetic input')
                    key = (row['algorithm'], row['parameter'], row['requested_budget'], row['half_subk'])
                    if key in configs[unit]: raise ValueError('duplicate synthetic configuration')
                    configs[unit].add(key); units[unit] += 1
                    invalid += row['status'] != 'ok'
                else:
                    typ = row['type']; counts[typ] += 1
                    if typ == 'pair':
                        if row['pair_id'] in pair_ids: raise ValueError('duplicate real pair')
                        pair_ids.add(row['pair_id'])
                    if typ == 'object':
                        if row['object_id'] in object_ids: raise ValueError('duplicate real object')
                        object_ids.add(row['object_id'])
                raw.write(line)
                if time.monotonic() - last > 3:
                    atomic(RUN / 'active_job.json', dict(job=destination.name, native_pid=proc.pid,
                        synthetic_units_complete=sum(n == 29 for n in units.values()),
                        synthetic_rows=sum(units.values()), real_counts=dict(counts)))
                    last = time.monotonic()
            code = proc.wait()
        finally:
            if proc.poll() is None: proc.terminate(); proc.wait()
        if code != 0: raise ValueError(f'native batch failed ({code}): {destination}')
    if expected_units is not None:
        if set(units) != expected_units or any(n != 29 for n in units.values()):
            raise ValueError('incomplete synthetic batch')
    elif pair_ids != expected_pairs or object_ids != expected_objects:
        raise ValueError('incomplete real batch')
    (destination / 'raw.jsonl.gz.tmp').replace(destination / 'raw.jsonl.gz')
    result = dict(command=command, rows=dict(counts), synthetic_units=len(units),
                  synthetic_rows=sum(units.values()), invalid_estimates=invalid,
                  manifest_sha256=sha(manifest), binary_sha256=sha(BINARY),
                  raw_sha256=sha(destination / 'raw.jsonl.gz'))
    atomic(completed, result)
    return result


def preflight(threads):
    out = RUN / 'preflight'; out.mkdir(exist_ok=True)
    units = []
    for k in CONFIG['k']:
        group = []
        for ratio, c in itertools.product([1, 10], [0., .5, 1.]):
            row = dict(unit_id=digest(['preflight', k, ratio, c])[:24], k=k,
                       n_a=10000 // ratio, n_b=10000, intersection=int(10000 // ratio * c),
                       replicate=0, data_seed=20261030, hash_seed=101)
            group.append(row); units.append(row)
        directory = out / str(k); directory.mkdir(exist_ok=True)
        manifest = directory / 'inputs.tsv'; write_units(manifest, pd.DataFrame(group), verify=True)
        execute(manifest, directory, threads, expected_units={r['unit_id'] for r in group})
    reference = next(r for r in units if r['k'] == 24 and r['n_a'] == 1000 and r['intersection'] == 500)
    with gzip.open(out / '24/raw.jsonl.gz', 'rt') as stream:
        rows = [json.loads(line) for line in stream]
    comparisons = []
    for alg in ['minhash', 'fracminhash', 'hll', 'setsketch', 'bindash', 'probminhash', 'kssd']:
        selected = [r for r in rows if r['unit_id'] == reference['unit_id'] and r['algorithm'] == alg
                    and (r['parameter'] == 2 if alg == 'kssd' else r['requested_budget'] == 8192)]
        assert len(selected) == 1
        row = selected[0]
        command = [str(ROOT / 'build_validated/synthetic_pilot'), '--algorithm', alg,
            '--parameter', str(row['parameter']), '--k', '24', '--n-a', '1000', '--n-b', '10000',
            '--intersection', '500', '--data-seed', '20261030', '--hash-seed', str(row['hash_seed']),
            '--half-subk', str(row['half_subk']), '--query-min-ms', '1', '--verify-data', '1']
        proc = subprocess.run(command, capture_output=True, text=True, check=True)
        old = json.loads([s for s in proc.stdout.splitlines() if s.startswith('{')][-1])
        if abs(old['estimate'] - row['jaccard']) > 1e-14:
            raise ValueError(f'frozen-adapter estimate mismatch: {alg}')
        for key in ['payload_bytes_a', 'payload_bytes_b', 'input_digest_a', 'input_digest_b']:
            if row[key] != old[key]: raise ValueError(f'frozen adapter {key} mismatch: {alg}')
        comparisons.append(dict(algorithm=alg, estimate=row['jaccard'], matches_frozen=True))
    atomic(out / 'verification.json', dict(exact_verified_inputs=len(units), native_comparisons=comparisons))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage', choices=['plan', 'preflight', 'synthetic', 'real', 'all'], default='all')
    ap.add_argument('--threads', type=int, default=64)
    args = ap.parse_args()
    RUN.mkdir(exist_ok=True)
    lock = (ROOT / 'runs/.active_benchmark.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    units = synthetic_units()
    plan_path = RUN / 'synthetic_units.tsv'
    if plan_path.exists():
        old = pd.read_csv(plan_path, sep='\t', dtype={'data_seed': 'uint64'})
        if old.unit_id.tolist() != units.unit_id.tolist(): raise ValueError('registered plan changed')
    else: units.to_csv(plan_path, sep='\t', index=False)
    inventory = real_inventory(); jobs = real_jobs(inventory)
    plan = dict(config=CONFIG, synthetic_units=len(units), synthetic_measurements=len(units) * 29,
                real_jobs=jobs, input_manifests={str(p): sha(p) for paths in inventory.values() for p in paths})
    old_plan = RUN / 'plan.json'
    if old_plan.exists() and json.loads(old_plan.read_text()) != plan: raise ValueError('plan changed')
    atomic(old_plan, plan)
    print(json.dumps(dict(synthetic_units=len(units), native_synthetic_rows=len(units) * 29,
                          real_jobs=len(jobs), natural_pairs=len(pd.read_csv(inventory['natural'][1], sep='\t')))), flush=True)
    if args.stage == 'plan': return
    snapshot = RUN / 'code_snapshot'; snapshot.mkdir(exist_ok=True)
    files = [HERE / n for n in ['config.json', 'PROTOCOL.md', 'run.py', 'build.py', 'accuracy_adapter.cpp']]
    files += [BUILD / n for n in ['build_manifest.json', 'real_adapter.cpp']]
    for path in files:
        target = snapshot / path.name
        if target.exists() and sha(target) != sha(path): raise ValueError('registered execution code changed')
        shutil.copyfile(path, target)
    native_source = ROOT / 'source/RabbitSketch-validated'
    source_hashes = {str(p.relative_to(native_source)): sha(p) for p in sorted((native_source / 'src').rglob('*')) if p.is_file()}
    identity = dict(binary_sha256=sha(BINARY), native_source_hashes=source_hashes,
        native_commit=subprocess.check_output(['git', '-C', str(native_source), 'rev-parse', 'HEAD'], text=True).strip(),
        config_sha256=sha(HERE / 'config.json'), threads=args.threads)
    identity_path = RUN / 'identity.json'
    if identity_path.exists() and json.loads(identity_path.read_text()) != identity: raise ValueError('native identity changed')
    atomic(identity_path, identity)
    if args.stage in ['all', 'preflight']:
        preflight(args.threads)
        print('Native preflight passed', flush=True)
    if args.stage == 'preflight': return
    if not (RUN / 'preflight/verification.json').exists(): raise ValueError('preflight required')
    if args.stage in ['all', 'synthetic']:
        for (k, nb), group in units.groupby(['k', 'n_b'], sort=True):
            directory = RUN / 'synthetic' / f'k{k}_n{nb}'; directory.mkdir(parents=True, exist_ok=True)
            manifest = directory / 'inputs.tsv'; write_units(manifest, group)
            print(f'Synthetic k={k}, larger cardinality={nb}, inputs={len(group)}', flush=True)
            execute(manifest, directory, args.threads, expected_units=set(group.unit_id))
    if args.stage in ['all', 'real']:
        # Full file inventory is hashed once, before any real measurements.
        paths = sorted({r['path'] for objects, _ in inventory.values()
                        for r in csv.DictReader(objects.open(), delimiter='\t')})
        hashes = {p: sha(p) for p in paths}
        original = RUN / 'real_input_sha256.json'
        if original.exists() and json.loads(original.read_text()) != hashes: raise ValueError('real inputs changed')
        atomic(original, hashes)
        for number, job in enumerate(jobs, 1):
            directory = RUN / 'real' / job['job_id']; directory.mkdir(parents=True, exist_ok=True)
            manifest = directory / 'manifest.tsv'
            manifest.write_text('#regimes-v1\n' + ''.join(f'{k}\t{v}\n' for k, v in job.items() if k not in ['phase', 'job_id']))
            pairs = set(pd.read_csv(job['pairs'], sep='\t').pair_id)
            objects = set(pd.read_csv(job['objects'], sep='\t').object_id)
            print(f'Real {number}/{len(jobs)}: {job["phase"]} {job["algorithm"]} p={job["parameter"]} k={job["k"]} seed={job["seed"]}', flush=True)
            execute(manifest, directory, args.threads, expected_pairs=pairs, expected_objects=objects)
            atomic(RUN / 'real_progress.json', dict(completed_jobs=number, planned_jobs=len(jobs), last_job=job))
    atomic(RUN / (args.stage + '_finished.json'), dict(status='complete', stage=args.stage))


if __name__ == '__main__':
    main()
