"""Scientific invariants: true sets, feasible regimes, and conservative decisions."""
import json
import os
from pathlib import Path
import subprocess

import numpy as np
import pytest

from analyze import aggregate, crossover, decisions, dominates, paired_error_upper
from run import HERE, DEFAULT_ROOT, command, make_plan


@pytest.fixture
def config():
    return json.loads((HERE / 'configs/smoke.json').read_text())


def test_feasible_jaccard_and_realized_truth(config):
    plan = make_plan(config)
    assert plan['skipped_regimes']
    for r in plan['regimes']:
        assert 0 <= r['intersection'] <= min(r['n_a'], r['n_b'])
        assert r['target_jaccard'] <= r['n_a'] / r['n_b']
        assert abs(r['ground_truth'] - r['target_jaccard']) <= 2 / (r['n_a'] + r['n_b'])
    assert all(j['n_b'] == 10000 for j in plan['jobs'])


def test_replicates_are_paired_across_every_backend(config):
    plan = make_plan(config)
    for r in plan['regimes']:
        seeds = []
        for rep in range(config['replicates']):
            one = {j['data_seed'] for j in plan['jobs'] if j['regime_id'] == r['regime_id'] and j['replicate'] == rep}
            assert len(one) == 1
            seeds.append(next(iter(one)))
        assert len(set(seeds)) == config['replicates']
    assert make_plan(config) == plan


@pytest.mark.parametrize('name,parameter', [('kssd', 4), ('bindash', 65), ('hll', 21),
                                           ('setsketch', 17), ('probminhash', 1)])
def test_native_resolution_validation(config, name, parameter):
    config['algorithms'][name]['parameters'] = [parameter]
    with pytest.raises(ValueError):
        make_plan(config)


def test_alias_cannot_be_a_second_algorithm(config):
    config['algorithms']['fastkmv'] = {'parameters': [64]}
    with pytest.raises(ValueError):
        make_plan(config)


def test_dominance_requires_all_costs_and_one_strict_gain():
    a = dict(mae=.01, payload_bytes=1024, build_seconds=1, query_seconds=.001)
    b = dict(a, mae=.02)
    assert dominates(a, b)
    assert not dominates(a, a)
    assert not dominates(dict(a, query_seconds=.002), b)


def artificial_rows(plan):
    rows = []
    for job in plan['jobs']:
        factor = 1 if job['algorithm'] == 'minhash' else 2
        rows.append(dict(job, status='ok', estimate=job['ground_truth'] + .001,
                         input_digest_a=job['data_seed'], input_digest_b=job['data_seed'] + 1,
                         sketch_payload_bytes=64*factor, build_seconds=.01*factor,
                         query_seconds=.00001*factor, setup_seconds=0,
                         query_cold_seconds=.00002*factor, process_peak_rss_bytes=100000))
    return rows


def test_pruning_needs_complete_independent_replicates(config):
    config.update(purpose='test_screening', k=[20], size_ratios=[1], jaccards=[.5], replicates=20)
    config['algorithms'] = {'minhash': {'parameters': [64]}, 'fracminhash': {'parameters': [50]}}
    plan = make_plan(config)
    rows = artificial_rows(plan)
    report = decisions(plan, aggregate(plan, rows))
    assert report[1]['status'] == 'review_pruning_candidate'
    report = decisions(plan, aggregate(plan, rows[:-1]))
    assert all(r['status'] != 'review_pruning_candidate' for r in report)
    config['purpose'] = 'smoke_only'
    plan = make_plan(config)
    assert all(r['status'] == 'insufficient_replicates_or_smoke'
               for r in decisions(plan, aggregate(plan, artificial_rows(plan))))


def test_missing_failed_empty_and_mismatched_inputs_are_not_silent(config):
    config.update(k=[20], size_ratios=[1], jaccards=[.5])
    config['algorithms'] = {'minhash': {'parameters': [64]}, 'fracminhash': {'parameters': [50]}}
    plan = make_plan(config)
    rows = artificial_rows(plan)
    rows[0]['status'] = 'timeout'
    summaries = aggregate(plan, rows)
    assert sum(s['failed'] for s in summaries) == 1
    assert not next(s for s in summaries if s['failed'])['pareto_frontier']
    rows = artificial_rows(plan)
    rows[0]['empty_sample_a'] = 1
    summaries = aggregate(plan, rows)
    assert not next(s for s in summaries if s['empty_samples'])['eligible']
    rows = artificial_rows(plan)
    rows[0]['input_digest_a'] = -1
    with pytest.raises(ValueError, match='same input sets'):
        aggregate(plan, rows)


def test_bootstrap_preserves_pairing():
    a = np.array([0., 1., 0., 1.])
    indices = np.random.default_rng(42).integers(0, 4, (2000, 4))
    assert paired_error_upper(a, a, indices) == 0
    assert paired_error_upper(a, a[::-1], indices) > 0


def test_positive_query_reuse_crossovers_only():
    assert crossover(dict(build_seconds=1, query_seconds=2), dict(build_seconds=11, query_seconds=1)) == 10
    assert crossover(dict(build_seconds=1, query_seconds=1), dict(build_seconds=11, query_seconds=1)) is None
    assert crossover(dict(build_seconds=1, query_seconds=1), dict(build_seconds=11, query_seconds=2)) is None


def native(algorithm, parameter, common=64, na=128, nb=256):
    binary = Path(os.environ.get('SCREENING_BINARY', DEFAULT_ROOT / 'build/synthetic_pilot'))
    if not binary.is_file():
        pytest.skip('native worker not built; set SCREENING_BINARY')
    job = dict(algorithm=algorithm, parameter=parameter, k=20, n_a=na, n_b=nb,
               intersection=common, data_seed=12345, hash_seed=42, query_min_ms=1,
               half_subk=3, bbits=16)
    proc = subprocess.run(command(binary, job) + ['--verify-data', '1'], capture_output=True, text=True, timeout=30)
    assert proc.returncode == 0, proc.stderr + proc.stdout
    return json.loads([s for s in proc.stdout.splitlines() if s.startswith('{')][-1])


@pytest.mark.parametrize('algorithm,parameter', [('minhash', 512), ('fracminhash', 1), ('kssd', 0)])
@pytest.mark.parametrize('common', [0, 64, 128])
def test_unsampled_backends_recover_exact_sets(algorithm, parameter, common):
    row = native(algorithm, parameter, common)
    assert row['verified_distinct_a'] == 128 and row['verified_distinct_b'] == 256
    assert row['verified_intersection'] == common
    assert row['estimate'] == pytest.approx(common / (384-common), abs=1e-12)


@pytest.mark.parametrize('algorithm,parameter', [('minhash', 64), ('fracminhash', 10), ('hll', 8),
                         ('setsketch', 8), ('kssd', 0), ('probminhash', 64), ('bindash', 64)])
def test_native_identical_sets_ignore_element_order(algorithm, parameter):
    row = native(algorithm, parameter, common=128, na=128, nb=128)
    assert row['input_digest_a'] != row['input_digest_b']  # independent permutations
    assert row['estimate'] == pytest.approx(1., abs=1e-12)


def test_sparse_empty_state_is_explicit():
    row = native('fracminhash', 4294967295)
    assert row['empty_sample_a'] and row['empty_sample_b']
    assert row['payload_bytes_a'] == row['payload_bytes_b'] == 0
