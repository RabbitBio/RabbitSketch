#!/usr/bin/env python3
"""Plan, run and resume isolated native sketch measurements (stdlib only)."""
import argparse
import datetime as dt
from fractions import Fraction
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import random
import subprocess
import sys

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = Path('/home/bigssd/zt/rabbitsketch2_benchmark')
ALGORITHMS = {'minhash', 'fracminhash', 'hll', 'setsketch', 'kssd', 'probminhash', 'bindash'}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    tmp.replace(path)


def positive_int(x):
    return isinstance(x, int) and not isinstance(x, bool) and x > 0


def validate(config):
    if config.get('schema_version') != 1:
        raise ValueError('unsupported schema_version')
    for key in ('k', 'larger_cardinalities', 'size_ratios', 'jaccards'):
        if not config.get(key) or len(set(config[key])) != len(config[key]):
            raise ValueError(f'{key} must be nonempty and unique')
    if any(not positive_int(k) or k > 32 for k in config['k']):
        raise ValueError('k must be an integer in [1,32]')
    for key in ('larger_cardinalities', 'size_ratios'):
        if any(not positive_int(x) for x in config[key]):
            raise ValueError(f'{key} must contain positive integers')
    if any(not isinstance(j, (int, float)) or not math.isfinite(j) or not 0 <= j <= 1
           for j in config['jaccards']):
        raise ValueError('Jaccard must be finite and in [0,1]')
    for key in ('replicates', 'query_min_ms', 'timeout_seconds'):
        if not positive_int(config.get(key)):
            raise ValueError(f'{key} must be a positive integer')
    if config['query_min_ms'] > 1000:
        raise ValueError('query_min_ms exceeds worker limit')
    for key in ('data_seed', 'hash_seed'):
        x = config.get(key)
        if not isinstance(x, int) or isinstance(x, bool) or not 0 <= x < 2**64:
            raise ValueError(f'{key} must be a uint64')
    if not config.get('algorithms') or set(config['algorithms']) - ALGORITHMS:
        raise ValueError('unknown/empty algorithm list; MinHash/FastKMV must be one entry')
    for name, spec in config['algorithms'].items():
        params = spec.get('parameters', [])
        if not params or len(params) != len(set(params)):
            raise ValueError(f'{name}: parameters must be nonempty and unique')
        for p in params:
            if not isinstance(p, int) or isinstance(p, bool) or not 0 <= p < 2**32:
                raise ValueError(f'{name}: native parameter must be uint32')
            valid = p > 0
            if name == 'probminhash':
                valid = p >= 2
            elif name == 'hll':
                valid = 4 <= p <= 20
            elif name == 'setsketch':
                valid = 4 <= p <= 16
            elif name == 'bindash':
                valid = p > 0 and p % 64 == 0 and spec.get('bbits', 16) in range(1, 17)
            elif name == 'kssd':
                h = spec.get('half_subk', 6)
                valid = isinstance(h, int) and 3 <= h <= 7 and p + 3 <= h
                valid = valid and all(k % 2 == 0 and k // 2 >= h for k in config['k'])
            if not valid:
                raise ValueError(f'{name}: illegal native parameter {p}')
    analysis = config['analysis']
    if not positive_int(analysis['minimum_replicates']) or not positive_int(analysis['bootstrap_samples']):
        raise ValueError('positive replicate/bootstrap counts required')
    for key in ('dominated_regime_fraction', 'resource_advantage_fraction'):
        if not 0 < analysis[key] <= 1:
            raise ValueError(f'invalid {key}')
    if not math.isfinite(analysis['error_equivalence_margin']) or analysis['error_equivalence_margin'] < 0:
        raise ValueError('invalid error margin')
    if any(not math.isfinite(x) or x <= 0 for x in analysis['error_constraints']):
        raise ValueError('invalid error constraint')
    if any(not positive_int(q) for q in analysis['query_reuse']):
        raise ValueError('invalid query reuse')


def make_plan(config):
    validate(config)
    jobs, regimes, skipped = [], [], []
    for k, nb, ratio, j in itertools.product(config['k'], config['larger_cardinalities'],
                                            config['size_ratios'], config['jaccards']):
        # "Cardinality" is that of the LARGER set: no unintended billion-element jobs.
        if nb % ratio:
            raise ValueError(f'larger cardinality {nb} not divisible by ratio {ratio}')
        na = nb // ratio
        if not na:
            raise ValueError('size ratio creates an empty set')
        item = dict(k=k, n_a=na, n_b=nb, size_ratio=ratio, target_jaccard=j)
        target = Fraction(str(j))
        if target > Fraction(na, nb):
            skipped.append(dict(item, reason='J > min(|A|,|B|)/max(|A|,|B|)'))
            continue
        exact_i = target * (na + nb) / (1 + target)
        intersection = (2 * exact_i.numerator + exact_i.denominator) // (2 * exact_i.denominator)
        if na + nb - intersection > (2**(2 * k) - 1) // 4:
            raise ValueError('set too large for k-mer universe')
        item.update(intersection=intersection,
                    ground_truth=intersection / (na + nb - intersection))
        item['regime_id'] = digest(item)[:16]
        if any(r['k'] == k and r['n_a'] == na and r['n_b'] == nb and
               r['intersection'] == intersection for r in regimes):
            raise ValueError('two J targets round to the same realized regime')
        regimes.append(item)
        for rep in range(config['replicates']):
            seed = int(digest([config['data_seed'], item['regime_id'], rep])[:16], 16)
            for name, spec in config['algorithms'].items():
                for p in spec['parameters']:
                    native = dict(algorithm=name, parameter=p, bbits=spec.get('bbits', 16),
                                  half_subk=spec.get('half_subk', 6))
                    job = dict(item, **native, replicate=rep, data_seed=seed,
                               hash_seed=config['hash_seed'], query_min_ms=config['query_min_ms'])
                    job['config_id'] = digest(native)[:16]
                    job['job_id'] = digest(job)[:24]
                    jobs.append(job)
    if not jobs:
        raise ValueError('no feasible regimes')
    random.Random(config['data_seed']).shuffle(jobs)  # avoid fixed algorithm/time confounding
    return dict(schema_version=1, config=config, config_sha256=digest(config),
                regimes=regimes, skipped_regimes=skipped, jobs=jobs)


def command(binary, job):
    fields = ('algorithm', 'parameter', 'k', 'n_a', 'n_b', 'intersection', 'data_seed',
              'hash_seed', 'query_min_ms', 'half_subk', 'bbits')
    return [str(binary)] + [v for key in fields for v in ('--' + key.replace('_', '-'), str(job[key]))]


def source_identity(source):
    def git(*args):
        return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()
    files = sorted(p for p in (source / 'src').rglob('*') if p.is_file())
    return dict(path=str(source), commit=git('rev-parse', 'HEAD'), status=git('status', '--porcelain'),
                src_sha256=digest([(str(p.relative_to(source)), file_sha(p)) for p in files]))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=HERE / 'configs/smoke.json')
    parser.add_argument('--binary', type=Path, default=DEFAULT_ROOT / 'build/synthetic_pilot')
    parser.add_argument('--source', type=Path, default=DEFAULT_ROOT / 'source/RabbitSketch')
    parser.add_argument('--output', type=Path, default=DEFAULT_ROOT / 'runs/smoke')
    parser.add_argument('--plan-only', action='store_true')
    parser.add_argument('--limit', type=int, help='run at most this many unfinished jobs (resume later)')
    parser.add_argument('--retry-failed', action='store_true')
    parser.add_argument('--cpu', type=int, help='pin native workers to one allowed CPU')
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit < 1:
        parser.error('--limit must be positive')
    plan = make_plan(json.loads(args.config.read_text()))
    summary = dict(feasible_regimes=len(plan['regimes']), skipped_regimes=len(plan['skipped_regimes']),
                   configurations=sum(len(s['parameters']) for s in plan['config']['algorithms'].values()),
                   jobs=len(plan['jobs']), input_element_visits=sum(j['n_a'] + j['n_b'] for j in plan['jobs']),
                   maximum_input_carrier_bytes=max((j['n_a'] + j['n_b']) * (j['k'] + 1) for j in plan['jobs']))
    print(json.dumps(summary, indent=2), flush=True)
    output = args.output.resolve()
    plan_path = output / 'plan.json'
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        raise ValueError('output contains a different plan; choose a new output directory')
    atomic_json(plan_path, plan)
    if args.plan_only:
        return 0
    binary, source = args.binary.resolve(), args.source.resolve()
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '1'
    env['OPENBLAS_NUM_THREADS'] = '1'
    if args.cpu is not None:
        if args.cpu not in os.sched_getaffinity(0):
            raise ValueError('--cpu is outside this process affinity mask')
        os.sched_setaffinity(0, {args.cpu})
    identity = dict(binary_sha256=file_sha(binary), source=source_identity(source),
                    runner_sha256=file_sha(__file__), worker_source_sha256=file_sha(HERE / 'synthetic_pilot.cpp'),
                    machine=platform.platform(), cpu_affinity=sorted(os.sched_getaffinity(0)),
                    simd_override=env.get('RABBITSKETCH_SIMD', ''), omp_num_threads=1,
                    cpu_model=next((s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines()
                                    if s.startswith('model name')), 'unknown'))
    manifest_path = output / 'manifest.json'
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        if old['identity'] != identity:
            raise ValueError('source, runner, binary or host changed; use a new output directory')
    else:
        atomic_json(manifest_path, dict(identity=identity, config_sha256=plan['config_sha256'],
                                       created_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
    completed, failed, executed = 0, 0, 0
    for job in plan['jobs']:
        path = output / 'rows' / (job['job_id'] + '.json')
        if path.exists():
            previous = json.loads(path.read_text())
            if previous.get('job_id') != job['job_id']:
                raise ValueError(f'corrupt row: {path}')
            if previous['status'] == 'ok' or not args.retry_failed:
                completed += 1
                failed += previous['status'] != 'ok'
                continue
        if args.limit is not None and executed >= args.limit:
            break
        cmd = command(binary, job)
        stderr = ''
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, env=env,
                                  timeout=plan['config']['timeout_seconds'])
            stderr = proc.stderr
            # Some native backends emit diagnostic lines on stdout.
            lines = [s for s in proc.stdout.splitlines() if s.startswith('{')]
            measured = json.loads(lines[-1]) if lines else dict(status='error', reason='worker emitted no JSON')
            if proc.returncode and measured.get('status') == 'ok':
                measured.update(status='error', reason=f'worker exit {proc.returncode}')
            measured['exit_code'] = proc.returncode
            atomic_json(path, dict(job, **{key: value for key, value in measured.items() if key not in job},
                                   measured_ground_truth=measured.get('ground_truth')))
        except subprocess.TimeoutExpired:
            atomic_json(path, dict(job, status='timeout', reason='worker time limit exceeded'))
        if stderr:
            (output / 'rows' / (job['job_id'] + '.stderr')).write_text(stderr)
        result = json.loads(path.read_text())
        if result['status'] == 'ok' and not math.isclose(result['ground_truth'], result['measured_ground_truth'], abs_tol=1e-15):
            result.update(status='error', reason='worker/plan truth mismatch')
            atomic_json(path, result)
        completed += 1; executed += 1
        failed += result['status'] != 'ok'
        print(f"[{completed}/{len(plan['jobs'])}] {job['algorithm']} p={job['parameter']} "
              f"k={job['k']} rep={job['replicate']} {result['status']}", flush=True)
    with (output / 'raw.jsonl').open('w') as out:
        for job in plan['jobs']:
            path = output / 'rows' / (job['job_id'] + '.json')
            if path.exists():
                out.write(canonical(json.loads(path.read_text())) + '\n')
    print(f'Executed {executed}; observed {completed}/{len(plan["jobs"])}; failed {failed}.', flush=True)
    return 1 if failed else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as exc:
        sys.exit(str(exc))
