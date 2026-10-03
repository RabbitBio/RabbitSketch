"""Check that main's OpenMP batch path preserves the statistical experiment."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from analyze import aggregate
from run import DEFAULT_ROOT, HERE, digest, make_plan
from run_main import group_jobs, write_manifest
from test_screening import artificial_rows


@pytest.fixture
def main_binary():
    binary = Path(os.environ.get('SCREENING_MAIN_BINARY', DEFAULT_ROOT / 'build_main/rabbitsketch_pilot'))
    if not binary.exists():
        pytest.skip('main adapter not built')
    return binary


def batch(binary, path, jobs, threads):
    write_manifest(path, jobs)
    proc = subprocess.run([str(binary), '--screening-manifest', str(path), '-t', str(threads)],
                          text=True, capture_output=True, timeout=60)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    records = [json.loads(line) for line in proc.stdout.splitlines() if line.startswith('{')]
    event = records.pop(0)
    assert event['threads'] == threads and event['pairs'] == len(jobs)
    assert {r['job_id'] for r in records} == {j['job_id'] for j in jobs}
    assert all(r['status'] == 'ok' and r['batch_threads'] == threads for r in records)
    return {r['job_id']: r for r in records}


@pytest.mark.parametrize('algorithm,parameter', [('minhash', 64), ('fracminhash', 50), ('hll', 8),
                        ('setsketch', 8), ('kssd', 1), ('probminhash', 64), ('bindash', 64)])
def test_thread_count_preserves_sets_and_estimators(main_binary, tmp_path, algorithm, parameter):
    jobs = []
    for seed in range(8):
        job = dict(algorithm=algorithm, parameter=parameter, k=20, n_a=1000, n_b=10000,
                   intersection=500, data_seed=seed, hash_seed=42, query_min_ms=1,
                   half_subk=4, bbits=16)
        job['job_id'] = digest(job)[:24]
        jobs.append(job)
    serial = batch(main_binary, tmp_path/'serial.tsv', jobs, 1)
    parallel = batch(main_binary, tmp_path/'parallel.tsv', jobs, 4)
    stable_fields = ('input_digest_a', 'input_digest_b', 'ground_truth', 'estimate',
                     'payload_bytes_a', 'payload_bytes_b', 'empty_sample_a', 'empty_sample_b',
                     'cardinality_estimate_a', 'cardinality_estimate_b')
    for job_id, row in parallel.items():
        assert {k: row[k] for k in stable_fields} == {k: serial[job_id][k] for k in stable_fields}
        assert row['measurement_profile'] == 'main_openmp_batch'
        assert row['rss_scope'] == 'shared_batch_process'
    assert len({r['worker_thread'] for r in parallel.values()}) == 4


@pytest.mark.parametrize('case', ['duplicate', 'mixed', 'malformed'])
def test_manifest_rejects_ambiguous_batches(main_binary, tmp_path, case):
    config = json.loads((HERE/'configs/smoke.json').read_text())
    jobs = group_jobs(make_plan(config)['jobs'])[0][:2]
    if case == 'duplicate':
        jobs[1] = jobs[0]
    if case == 'mixed':
        jobs[1] = dict(jobs[1], k=28)
    path = tmp_path/'bad.tsv'
    write_manifest(path, jobs)
    if case == 'malformed':
        path.write_text('invalid\tmanifest\n')
    proc = subprocess.run([str(main_binary), '--screening-manifest', str(path), '-t', '4'],
                          capture_output=True, text=True, timeout=15)
    assert proc.returncode != 0
    assert 'screening batch:' in proc.stderr
    assert 'batch_event' not in proc.stdout


def test_runner_resumes_and_waits_for_each_main_process(main_binary, tmp_path):
    config = json.loads((HERE/'configs/smoke.json').read_text())
    config.update(k=[20], larger_cardinalities=[1000], size_ratios=[1], jaccards=[.5],
                  replicates=4, query_min_ms=1)
    config['algorithms'] = {'minhash': {'parameters': [64]}, 'fracminhash': {'parameters': [50]}}
    cfg = tmp_path/'config.json'
    cfg.write_text(json.dumps(config))
    out = tmp_path/'output'
    cmd = [sys.executable, str(HERE/'run_main.py'), '--config', str(cfg), '--output', str(out),
           '--binary', str(main_binary), '--threads', '4']
    first = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    assert first.returncode == 0, first.stdout + first.stderr
    rows = {p.name: p.read_bytes() for p in (out/'rows').glob('*.json')}
    assert len(rows) == 8
    batches = sorted((json.loads(p.read_text()) for p in (out/'batches').glob('*.json')),
                     key=lambda b: b['started_utc'])
    assert len(batches) == 2
    assert batches[0]['finished_utc'] <= batches[1]['started_utc']
    assert all(b['actual_threads'] == 4 and b['exit_code'] == 0 for b in batches)
    second = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    assert second.returncode == 0, second.stdout + second.stderr
    assert 'Executed 0;' in second.stdout
    assert rows == {p.name: p.read_bytes() for p in (out/'rows').glob('*.json')}
    changed = subprocess.run(cmd[:-1]+['2'], capture_output=True, text=True, timeout=30)
    assert changed.returncode != 0 and 'execution profile' in changed.stderr


def test_analysis_rejects_mixing_serial_and_batch_timings():
    config = json.loads((HERE/'configs/smoke.json').read_text())
    plan = make_plan(config)
    rows = artificial_rows(plan)
    rows[0].update(measurement_profile='main_openmp_batch', batch_threads=4)
    with pytest.raises(ValueError, match='mix measurement profiles'):
        aggregate(plan, rows)
