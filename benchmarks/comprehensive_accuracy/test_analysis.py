"""Guard independent cohorts, missing outputs and calibration leakage."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import pytest
from common import add_errors, summaries, paired_family, candidates, bounded_containment
from analyze_synthetic import choose_kssd
from analyze_real import fixed_parameters
from downstream import adjusted_rand, single_linkage, linkage_curve, nearest_scores, matrix, calibrate_threshold
import analyze_real


def frame():
    return pd.DataFrame(dict(algorithm=['fracminhash'] * 4, jaccard=[0., .5, .2, 0.],
        truth_j=[0., .4, .2, 0.], truth_c=[0., .8, .4, 0.],
        containment_ab=[0., 1., .4, 0.], native_containment_ab=[0., 1., .4, 0.],
        cardinality_a=[10., 10., 10., 0.], cardinality_b=[10., 10., 10., 0.],
        n_a=[10] * 4, n_b=[10] * 4, payload_bytes_a=[8, 8, 8, 0],
        payload_bytes_b=[8, 8, 8, 0], empty_a=[0, 0, 0, 1], empty_b=[0, 0, 0, 1],
        budget=[16] * 4, group=['x'] * 4))


def test_empty_zero_is_not_a_perfect_estimate():
    d = add_errors(frame())
    assert d.j_error.iloc[0] == 0
    assert np.isnan(d.j_error.iloc[3])
    s = summaries(d, ['group'])
    row = s[s.metric == 'j_error'].iloc[0]
    assert row.availability == .75 and row.valid == 3
    assert row.mae == pytest.approx(.1 / 3)


def test_unsupported_metrics_stay_unavailable():
    f = frame(); f['algorithm'] = 'kssd'
    d = add_errors(f)
    assert d.c_error.isna().all() and d.card_relative_error.isna().all()
    assert d.j_error.notna().sum() == 3


def test_infinite_cardinality_does_not_create_containment():
    c = bounded_containment([.2, .9, .2], [np.inf, 10., 0.], [10., 100., 10.])
    assert np.isnan(c[0]) and np.isnan(c[2])
    assert c[1] == 1


def test_relative_error_undefined_at_zero_truth():
    f = frame(); f['truth_j'] = .5
    s = summaries(add_errors(f), ['group'])
    j = s[s.metric == 'j_error'].iloc[0]
    assert j.relative_mae == pytest.approx(j.mae / .5)
    f['truth_j'] = 0.
    s = summaries(add_errors(f), ['group'])
    assert np.isnan(s[s.metric == 'j_error'].relative_mae.iloc[0])


def test_kssd_parameter_never_uses_validation_size_or_error():
    rows = []
    for p, size in [(1, 600), (2, 100), (3, 8)]:
        for split in ['calibration', 'validation']:
            rows.append(dict(regime_id='r', parameter=p, half_subk=6, split=split,
                algorithm='kssd', pair_bytes=size if split == 'calibration' else size * 10,
                requested_budget=0))
    d = pd.DataFrame(rows)
    selected, choices = choose_kssd(d)
    chosen = choices[choices.budget == 512].iloc[0]
    assert chosen.parameter == 2
    # Validation overflow is retained, not repaired with a lower-resolution sketch.
    assert selected[(selected.budget == 512) & (selected.split == 'validation')].pair_bytes.iloc[0] == 1000


def test_duplicate_budget_preserves_bootstrap_dependence():
    a, b = np.array([.2, .1, .3, .15]), np.array([.3, .4, .2, .25])
    one, c1 = paired_family({8: {'a': a, 'b': b}}, 42, samples=1000)
    two, c2 = paired_family({8: {'a': a, 'b': b}, 32: {'a': a, 'b': b}}, 42, samples=1000)
    assert c1 == c2
    assert one[0]['simultaneous_low'] == two[1]['simultaneous_low']


def test_missing_paired_value_is_rejected():
    with pytest.raises(ValueError):
        paired_family({8: {'a': np.array([.1, np.nan]), 'b': np.array([.2, .3])}}, 1, 100)
    with pytest.raises(ValueError):
        paired_family({8: {'a': np.array([.1]), 'b': np.array([.2, .3])}}, 1, 100)


def test_candidate_set_requires_win_over_all_other_candidates():
    def contrast(a, b):
        return dict(left=a, right=b, mae_left=.01, mae_right=.1,
                    simultaneous_low=-.12, simultaneous_high=-.06, se=.01)
    remaining, winner = candidates(['a', 'b', 'c'], [contrast('a', 'b')])
    assert remaining == ['a', 'c'] and winner == ''
    remaining, winner = candidates(['a', 'b', 'c'], [contrast('a', 'b'), contrast('a', 'c')])
    assert remaining == ['a'] and winner == 'a'


def test_tiny_or_degenerate_advantages_do_not_create_winners():
    base = dict(left='a', right='b', mae_left=.00001, mae_right=.00002,
                simultaneous_low=-.00002, simultaneous_high=-.000001, se=.00001)
    assert candidates(['a', 'b'], [base]) == (['a', 'b'], '')
    base.update(mae_left=.01, mae_right=.1, simultaneous_low=-.1, simultaneous_high=-.08, se=0)
    assert candidates(['a', 'b'], [base]) == (['a', 'b'], '')


def test_real_parameter_selection_keeps_heldout_failure():
    rows = []
    for p, size in [(1000, 1000), (5000, 200), (10000, 100)]:
        for split in ['calibration', 'validation']:
            rows.append(dict(algorithm='fracminhash', parameter=p, half_subk=6, split=split,
                             pair_bytes=size if split == 'calibration' else size * 3))
    d, choices = fixed_parameters(pd.DataFrame(rows))
    assert choices[choices.budget == 512].parameter.tolist() == [5000]
    assert d[(d.budget == 512) & (d.split == 'validation')].pair_bytes.tolist() == [600]


@pytest.mark.parametrize('labels,clusters,expected', [
    ([0, 0, 1, 1], [3, 3, 7, 7], 1.),
    ([0, 0, 1, 1], [0, 1, 0, 1], -.5),
    ([0, 0, 1, 1], [0, 0, 0, 0], 0.),
    ([0, 0, 0], [1, 1, 1], 1.),
    ([0, 1, 2], [4, 5, 6], 1.),
])
def test_ari_known_partitions(labels, clusters, expected):
    assert adjusted_rand(labels, clusters) == pytest.approx(expected)


def test_linkage_curve_matches_independent_threshold_components():
    rng = np.random.default_rng(42)
    s = rng.random((15, 15)); s = (s + s.T) / 2; np.fill_diagonal(s, 1)
    labels = np.repeat([0, 1, 2], 5); thresholds = np.arange(1, 100) / 100
    result = linkage_curve(s, labels, thresholds)
    expected = [adjusted_rand(labels, single_linkage(s, t)) for t in thresholds]
    assert np.allclose(result, expected)


def test_neighbor_ties_use_fractional_credit_and_exclude_self():
    scores = np.ones((4, 4)); labels = np.array([0, 0, 1, 1])
    exact = np.array([[1, .9, .1, .1], [.9, 1, .1, .1], [.1, .1, 1, .9], [.1, .1, .9, 1]])
    species, recovery, ties = nearest_scores(scores, labels, exact)
    assert species == pytest.approx(1 / 3) and recovery == pytest.approx(1 / 3) and ties == 3


def test_threshold_never_uses_validation_labels():
    curves = np.array([[.1, .8, .2], [.9, .0, 1.]])
    assert calibrate_threshold(curves, ['calibration', 'validation'], [.1, .2, .3]) == 1
    curves[1] = [100, -100, -100]
    assert calibrate_threshold(curves, ['calibration', 'validation'], [.1, .2, .3]) == 1


def test_missing_or_duplicate_edges_rejected():
    pairs = pd.DataFrame([dict(a='a', b='b', j=.5), dict(a='a', b='c', j=.3), dict(a='b', b='c', j=.7)])
    assert matrix(pairs, ['a', 'b', 'c'], 'j')[0, 1] == .5
    with pytest.raises(ValueError): matrix(pairs.iloc[:2], ['a', 'b', 'c'], 'j')
    pairs.loc[2, ['a', 'b']] = ['a', 'b']
    with pytest.raises(ValueError): matrix(pairs, ['a', 'b', 'c'], 'j')


def real_fixture():
    exact_base = dict(algorithm='exact', parameter=1, half_subk=6, seed=42, k=24)
    native_base = dict(algorithm='minhash', parameter=32, half_subk=6, seed=42, k=24)
    exact = dict(pair=pd.DataFrame([dict(exact_base, pair_id='p', n_a=100, n_b=200, intersection=50)]),
        object=pd.DataFrame([dict(exact_base, object_id='a', n=100, input_digest=2**64 - 1),
                             dict(exact_base, object_id='b', n=200, input_digest=2**64 - 2)]))
    native = dict(pair=pd.DataFrame([dict(native_base, pair_id='p', a='a', b='b', jaccard=.2,
        containment_native_ab=.5, containment_derived_ab=.5, containment_derived_ba=.25,
        payload_bytes_a=256, payload_bytes_b=256, empty_a=0, empty_b=0)]),
        object=pd.DataFrame([dict(native_base, object_id='a', cardinality=100., input_digest=2**64 - 1),
                             dict(native_base, object_id='b', cardinality=200., input_digest=2**64 - 2)]))
    meta = pd.DataFrame([dict(pair_id='p', a='a', b='b', genus_a=1, genus_b=1,
        species_a=2, species_b=3, split='validation', cluster='1')])
    return exact, native, meta


def test_real_loader_joins_exact_counts_and_native_cards_without_oracle(monkeypatch):
    exact, native, meta = real_fixture()
    jobs = [dict(phase='natural', k=24, algorithm='exact'), dict(phase='natural', k=24, algorithm='minhash')]
    monkeypatch.setattr(analyze_real, 'raw_job', lambda j: exact if j['algorithm'] == 'exact' else native)
    monkeypatch.setattr(analyze_real, 'metadata', lambda phase: meta)
    d, truth, _ = analyze_real.load_panel(jobs, 'natural', 24)
    assert d.j_error.iloc[0] == pytest.approx(0.)
    assert d.c_error.iloc[0] == pytest.approx(0.)
    assert d.card_relative_error.iloc[0] == 0.
    # Alter native cardinality: the estimate must not silently use exact cardinality.
    native['object'].loc[0, 'cardinality'] = 200.
    d, _, _ = analyze_real.load_panel(jobs, 'natural', 24)
    assert d.c_error.iloc[0] > .1
    assert d.card_relative_error.iloc[0] == .5


def test_real_loader_rejects_one_bit_fingerprint_mismatch(monkeypatch):
    exact, native, meta = real_fixture()
    native['object'].loc[0, 'input_digest'] = np.uint64(2**64 - 3)
    jobs = [dict(phase='natural', k=24, algorithm='exact'), dict(phase='natural', k=24, algorithm='minhash')]
    monkeypatch.setattr(analyze_real, 'raw_job', lambda j: exact if j['algorithm'] == 'exact' else native)
    monkeypatch.setattr(analyze_real, 'metadata', lambda phase: meta)
    with pytest.raises(ValueError, match='inputs differ'):
        analyze_real.load_panel(jobs, 'natural', 24)
