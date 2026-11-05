"""Independent checks for leakage barriers, failure accounting, and task semantics."""
import numpy as np
import pandas as pd
import pytest

from advisor import capability, nearest_distinct, payload_guard
from measure import observable_features, score_group, square, statistics
from prepare import derive, sample_segments
from scenario_common import CONFIG, bounded_containment, fingerprint


def test_fingerprint_ignores_record_order_and_strand():
    records = [b'AAATGC', b'ACGGT', b'NNAC']
    complement = bytes.maketrans(b'ACGTN', b'TGCAN')
    reverse = [s.translate(complement)[::-1] for s in records[::-1]]
    assert fingerprint(records) == fingerprint(reverse)
    assert fingerprint(records) != fingerprint(records + [b'ACGGT'])


def test_no_sampling_across_record_boundaries():
    records = [b'A' * 200, b'C' * 210]
    reads = sample_segments(records, 150, 100, np.random.default_rng(32))
    assert set(reads) == {b'A' * 150, b'C' * 150}


def test_read_error_variants_share_positions_and_nested_substitutions():
    records = [b'ACGT' * 2000, b'TGCA' * 1700]
    conditions = [dict(name='test', method='reads', coverage=1, length=150, substitution=p)
                  for p in [0, .01, .03]]
    results = [derive(records, records, c, 'source') for c in conditions]
    a, b, c = [np.frombuffer(b''.join(x[0]), dtype=np.uint8) for x in results]
    assert len(a) == len(b) == len(c)
    assert np.all((a != b) <= (a != c))
    assert np.array_equal(b[a != b], c[a != b])
    assert int(np.sum(a != b)) == results[1][1]['substitutions']
    assert int(np.sum(a != c)) == results[2][1]['substitutions']


def test_fragment_truncation_is_explicit():
    result, info = derive([b'ACGT' * 50], [], dict(method='fragment', length=1000), 'short')
    assert result == [b'ACGT' * 50]
    assert info['truncated']


def test_contamination_replaces_tiles_and_preserves_total_bases():
    source, donor = [b'A' * 1000000], [b'C' * 1000000]
    condition = dict(method='tiles', fraction=.1, tile_length=10000, contamination=.2)
    result, info = derive(source, donor, condition, 'contaminated')
    assert sum(map(len, result)) == 100000
    assert sum(s == b'C' * 10000 for s in result) == 2
    assert info['contaminant_records'] == 2


def test_containment_preserves_query_direction_and_common_intersection():
    a = bounded_containment([.1], [100], [1000])[0]
    b = bounded_containment([.1], [1000], [100])[0]
    assert a == pytest.approx(1) and b == pytest.approx(.1)
    assert a * 100 == pytest.approx(b * 1000)
    assert np.isnan(bounded_containment([.1], [0], [100])[0])


def test_missing_and_over_budget_are_failures_even_with_zero_error():
    result = statistics([0., 0., .03], [True, False, True], [True, False, False], [600, 400, 400])
    assert result['loss'] == pytest.approx(2.03 / 3)
    assert result['availability'] == pytest.approx(1 / 3)
    assert result['mae'] == .03
    assert result['budget_violation'] == pytest.approx(1 / 3)


def test_cardinality_clipping_only_changes_selection_loss():
    result = statistics([4., .2], [True, True], [False, False], [20, 20], clip=True)
    assert result['loss'] == pytest.approx(.6)
    assert result['mae'] == pytest.approx(2.1)
    assert result['p95'] > 3


def test_square_edge_order_matches_sorted_query_pairs():
    expected = np.array([[1, .1, .2], [.1, 1, .3], [.2, .3, 1]])
    assert np.array_equal(square([.1, .2, .3], 3), expected)
    with pytest.raises(ValueError): square([.1, .2], 3)


def test_nearest_evidence_uses_distinct_genera_and_matching_domain():
    x = np.array([[0.], [.01], [.02], [.1], [.2], [.3]])
    groups = np.array(['a', 'a', 'a', 'b', 'c', 'd'])
    indices, _, _, _ = nearest_distinct({'x': 0}, x, ['x'], np.array([1.]), 24, 'reads',
        groups, np.array([24, 24, 24, 24, 24, 32]), np.array(['reads'] * 6), count=3)
    assert indices.tolist() == [0, 3, 4]
    assert len(set(groups[indices])) == 3


def test_out_of_support_is_flagged():
    indices, supported, _, outside = nearest_distinct({'x': 100}, np.array([[0.], [1.]]), ['x'],
        np.array([1.]), 24, 'reads', np.array(['a', 'b']), np.array([24, 24]),
        np.array(['reads', 'reads']), count=2)
    assert len(indices) == 2 and not supported and outside == ['x']


@pytest.mark.parametrize('metric', ['cardinality', 'containment'])
@pytest.mark.parametrize('algorithm', ['kssd', 'probminhash_unweighted'])
def test_unsupported_estimands_cannot_be_selected(algorithm, metric):
    assert not capability(algorithm, metric)
    assert capability(algorithm, 'jaccard')


def test_variable_payload_guard_uses_the_task_objects():
    f = dict(log_query_bases_max=3, log_reference_bases_max=6,
             log_pilot_query_n_max=3, log_pilot_reference_n_max=6)
    pair = payload_guard(f, 24, 'fracminhash', 1000, 'jaccard', 'metadata')
    cluster = payload_guard(f, 24, 'fracminhash', 1000, 'clustering', 'metadata')
    assert pair == 8008 and cluster == 16
    assert payload_guard(f, 24, 'fracminhash', 1000, 'jaccard', 'pilot') == pair * 1.5


def test_observable_features_ignore_truth_and_taxonomy_columns():
    q, r = ['q1', 'q2'], ['r1', 'r2']
    work = dict(queries=q, references=r, retrieval_pairs=['a', 'b', 'c', 'd'])
    meta = pd.DataFrame(dict(object_id=q+r, bases=[100, 100, 1000, 1000],
        records=[1, 1, 1, 1], median_record_length=[100, 100, 1000, 1000])).set_index('object_id')
    pairs = pd.DataFrame(dict(pair_id=['a','b','c','d'], a=['q1','q1','q2','q2'],
        b=['r1','r2','r1','r2'], jaccard=[.01,.02,.03,.04])).set_index('pair_id')
    objects = pd.DataFrame(dict(object_id=q+r, cardinality=[90, 90, 900, 900], payload_bytes=[256]*4)).set_index('object_id')
    original = observable_features(work, meta, pairs, objects)
    pairs['truth_j'] = [1, 0, 0, 1]; objects['n'] = [1, 1, 100000000, 1]
    meta['species'] = ['hidden'] * 4; meta['true_error_rate'] = [1, 1, 1, 1]
    assert observable_features(work, meta, pairs, objects) == original


def task_fixture():
    labels = np.repeat([1, 2, 3], 3); ref = np.repeat([1, 2, 3], 2)
    truth = np.where(labels[:, None] == ref[None, :], .8, .1)
    cluster = np.where(labels[:, None] == labels[None, :], .8, .1); np.fill_diagonal(cluster, 1)
    w = dict(workload_id='w', kind='assembly')
    record = dict(workload=w, seed=42, j=truth.copy(), c=truth.copy(), available=np.ones((9,6),bool),
        pair_bytes=np.ones((9,6))*512, qq=cluster, qqvalid=True, qqbytes=np.ones(36)*512,
        cards=np.ones(9)*100, cardgood=np.ones(9,bool), cardbytes=np.ones(9)*512)
    target = dict(j=truth, c=truth, n=np.ones(9)*100, qq=cluster, query_species=labels, reference_species=ref)
    thresholds = dict(assembly=dict(clustering_threshold=.5, acceptance_threshold=.5))
    return record, {'w': target}, thresholds


def test_perfect_tasks_and_open_set_exclusion():
    r,t,th = task_fixture(); scores, extra = score_group([r], t, th, 512)
    assert all(scores[m]['loss'] == 0 for m in scores)
    assert extra['known_acceptance'] == 1
    assert extra['unknown_false_acceptance'] == 0
    assert extra['ari'] == 1


def test_uniform_retrieval_ties_receive_fractional_credit():
    r,t,th = task_fixture(); r['j'][:] = 0
    scores, extra = score_group([r], t, th, 512)
    assert scores['species']['loss'] == pytest.approx(2/3)
    assert scores['neighbor']['loss'] == pytest.approx(2/3)


def test_invalid_edge_fails_whole_query_and_clustering():
    r,t,th = task_fixture(); r['available'][0, 0] = False; r['qqvalid'] = False
    scores, extra = score_group([r], t, th, 512)
    assert scores['neighbor']['loss'] == pytest.approx(1/9)
    assert scores['species']['availability'] == pytest.approx(8/9)
    assert scores['clustering']['loss'] == 1
    assert scores['clustering']['availability'] == 0


def test_budget_failure_is_not_rescued_by_perfect_ranking():
    r,t,th = task_fixture(); scores, extra = score_group([r], t, th, 256)
    assert all(scores[m]['loss'] == 1 for m in scores)
    assert extra['known_acceptance'] == 0
    assert extra['open_availability'] == 0


def test_unknown_gallery_does_not_include_removed_species_failure():
    r,t,th = task_fixture()
    r['available'][0, 0] = False
    th['assembly']['acceptance_threshold'] = .05
    _, extra = score_group([r], t, th, 512)
    assert extra['open_availability'] == pytest.approx(8/9)
    assert extra['unknown_availability'] == 1
    assert extra['unknown_false_acceptance'] == 1


def test_zero_estimator_diagnostic_is_preserved():
    r,t,th = task_fixture(); r['j'][:] = 0
    scores, _ = score_group([r], t, th, 512)
    assert scores['jaccard']['mae'] == pytest.approx(scores['jaccard']['constant0_mae'])
    assert scores['jaccard']['positive_pair_relative_mae'] == 1


def test_condition_aggregates_do_not_inflate_biological_sample_count():
    from evaluate import expanded_conditions
    frame = pd.DataFrame(dict(genus=[1,1,2,2],condition=['full','unseen_fragment_5k']*2))
    result = expanded_conditions(frame)
    assert len(result[result.condition=='all']) == 4
    assert result[result.condition=='all'].genus.nunique() == 2
    assert len(result[result.condition=='original_conditions']) == 2
    assert len(result[result.condition=='unseen_conditions']) == 2


def test_primary_intervals_use_paired_genus_units(tmp_path,monkeypatch):
    import evaluate
    monkeypatch.setattr(evaluate,'OUT',tmp_path)
    rows=[]
    for genus,change in enumerate([-.01,0,.01],1):
        for condition in ['full','unseen_fragment_5k']:
            for policy,loss in [('fixed',.2),('pilot',.2+change),('metadata',.2+2*change)]:
                rows.append(dict(genus=genus,condition=condition,policy=policy,loss=loss,
                                 metric='jaccard',k=24,budget=8192))
    result=evaluate.comparisons(pd.DataFrame(rows))
    primary=result[(result.condition=='all') & (result.policy=='pilot')].iloc[0]
    assert primary.genera==3 and primary.family_size==1
    assert primary.difference==pytest.approx(0)
    assert primary.simultaneous_low<0<primary.simultaneous_high
    assert not primary.material_improvement


def test_degenerate_candidate_contrast_is_not_a_supported_win(tmp_path,monkeypatch):
    import evaluate
    monkeypatch.setattr(evaluate,'OUT',tmp_path)
    rows=[]
    for genus in [1,2,3]:
        for policy,loss in [('minhash',.1),('hll',.2)]:
            rows.append(dict(genus=genus,condition='full',policy=policy,loss=loss,
                mae=loss,availability=1,budget_violation=0,metric='jaccard',k=24,budget=8192))
    result=evaluate.candidate_map(pd.DataFrame(rows))
    assert result.iloc[0].candidate_count==2
    assert result.iloc[0].best_observed=='minhash'


def test_target_failure_distinguishes_selection_from_budget_infeasibility(tmp_path,monkeypatch):
    import evaluate
    monkeypatch.setattr(evaluate,'OUT',tmp_path)
    rows=[];alternatives=[]
    for wid,possible in [('good_alternative',True),('no_alternative',False)]:
        base=dict(workload_id=wid,k=24,budget=8192,metric='jaccard')
        rows.append(dict(base,genus=wid,condition='full',policy='pilot',supported=True,
            predicted_availability=1,predicted_budget_violation=0,predicted_mae=.01,predicted_p95=.02,
            availability=1,budget_violation=0,mae=.1,p95=.2))
        alternatives.extend([dict(base,availability=1,budget_violation=0,mae=.1,p95=.2),
            dict(base,availability=1,budget_violation=0,mae=.001 if possible else .1,p95=.002 if possible else .2)])
    summary=evaluate.precision_profiles(pd.DataFrame(rows),pd.DataFrame(alternatives))
    standard=summary[(summary.condition=='full')&(summary.profile=='standard')].iloc[0]
    assert standard.false_feasible_rate==1
    assert standard.no_measured_feasible_rate==.5
    assert standard.selection_miss_rate==1


def test_candidate_map_marks_no_target_feasible_without_deleting_rankings(tmp_path,monkeypatch):
    import evaluate
    monkeypatch.setattr(evaluate,'OUT',tmp_path)
    rows=[]
    for genus in [1,2,3]:
        rows.append(dict(genus=genus,condition='full',policy='minhash',loss=.2,mae=.2,
            availability=1,budget_violation=0,metric='containment',k=24,budget=8192))
    profiles=pd.DataFrame([dict(condition='full',policy='pilot',profile='standard',metric='containment',
        k=24,budget=8192,feasible_alternative_rate=0.,no_measured_feasible_rate=1.)])
    result=evaluate.candidate_map(pd.DataFrame(rows),profiles).iloc[0]
    assert result.candidates=='minhash'
    assert result.standard_target_status=='no_measured_configuration_meets_target'
