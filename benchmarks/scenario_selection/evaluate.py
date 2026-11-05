#!/usr/bin/env python3
"""Evaluate frozen decisions once; retain abstentions, failures and negative results."""
import itertools

from advisor import read_scores
from scenario_common import *

KEYS = ['workload_id', 'k', 'budget', 'metric']


def verify_freezes():
    for name in ['model_frozen.json', 'decisions_frozen.json']:
        f = json.loads((RUN / name).read_text())
        for path, digest in f['files'].items():
            if sha(path) != digest: raise ValueError('frozen artifact changed: ' + path)
    decisions_time = (RUN / 'decisions_frozen.json').stat().st_mtime
    plan = json.loads((RUN / 'validation-outcomes_plan.json').read_text())
    first_manifest = min((RUN / 'native' / j['job_id'] / 'manifest.tsv').stat().st_mtime for j in plan['jobs'])
    if first_manifest < decisions_time: raise ValueError('validation outcomes preceded frozen choices')
    return dict(model_frozen_sha256=sha(RUN / 'model_frozen.json'),
        decisions_frozen_sha256=sha(RUN / 'decisions_frozen.json'),
        validation_outcomes_after_decisions=True, first_outcome_manifest_timestamp=first_manifest,
        decisions_timestamp=decisions_time)


def match_decisions():
    decisions = pd.read_csv(OUT / 'validation_decisions.csv.gz', keep_default_na=True)
    decisions['config'] = decisions.config.fillna('')
    native = read_scores('validation')
    truth = native[native.algorithm == 'exact'].copy()
    native = native[native.algorithm != 'exact'].copy()
    if native.duplicated(KEYS + ['config']).any(): raise ValueError('duplicated outcome key')
    take = [c for c in native.columns if c not in decisions.columns or c in KEYS + ['config']]
    joined = decisions.merge(native[take], on=KEYS + ['config'], how='left', validate='many_to_one', indicator=True)
    missing = joined['_merge'] != 'both'
    if (missing & joined.config.ne('')).any(): raise ValueError('selected configuration has no measured outcome')
    joined.loc[missing, 'loss'] = 1.; joined.loc[missing, 'availability'] = 0.
    joined.loc[missing, 'budget_violation'] = 0.; joined.loc[missing, 'max_payload'] = 0.
    joined['algorithm'] = joined.algorithm.fillna('unsupported'); joined = joined.drop(columns='_merge')
    oracle = native.groupby(KEYS).loss.min().rename('oracle_loss').reset_index()
    feasible = native[(native.max_payload <= native.budget) & (native.availability >= .99)]
    feasible = feasible.groupby(KEYS).loss.min().rename('feasible_oracle_loss').reset_index()
    joined = joined.merge(oracle, on=KEYS, validate='many_to_one').merge(feasible, on=KEYS, how='left', validate='many_to_one')
    joined['regret'] = joined.loss - joined.oracle_loss
    if joined.regret.min() < -1e-8: raise ValueError('decision outperformed its exhaustive penalized oracle')
    joined['near_oracle'] = joined.regret <= np.maximum(.0001, .1 * joined.oracle_loss)
    joined['beats_constant0'] = joined.loss < joined.constant0_mae - .0001
    joined['beats_constant1'] = joined.loss < joined.constant1_mae - .0001
    # This is recommendation coverage, not an accuracy guarantee.
    joined['supported'] = joined.supported.fillna(False).astype(bool)
    save_csv(OUT / 'decision_outcomes.csv.gz', joined)
    save_csv(OUT / 'exact_task_reference.csv.gz', truth)
    return joined, native, truth


def expanded_conditions(frame):
    unseen = frame.condition.str.startswith('unseen_')
    return pd.concat([frame, frame.assign(condition='all'),
        frame[~unseen].assign(condition='original_conditions'),
        frame[unseen].assign(condition='unseen_conditions')], ignore_index=True)


def summarize(decisions):
    d = expanded_conditions(decisions)
    keys = ['metric', 'k', 'budget', 'condition', 'policy']
    values = dict(loss=('loss', 'mean'), conditional_mae=('mae', 'mean'),
        mean_within_workload_p95=('p95', 'mean'), availability=('availability', 'mean'),
        budget_violation=('budget_violation', 'mean'), max_payload=('max_payload', 'max'),
        mean_regret=('regret', 'mean'), near_oracle_fraction=('near_oracle', 'mean'),
        support_fraction=('supported', 'mean'), mean_predicted_loss=('predicted_loss', 'mean'),
        zero_fraction=('zero_fraction', 'mean'), constant0_mae=('constant0_mae', 'mean'),
        constant1_mae=('constant1_mae', 'mean'), positive_pair_relative_mae=('positive_pair_relative_mae', 'mean'),
        beats_constant0_fraction=('beats_constant0', 'mean'), beats_constant1_fraction=('beats_constant1', 'mean'),
        workload_count=('workload_id', 'size'), genera=('genus', 'nunique'))
    summary = d.groupby(keys).agg(**values).reset_index()
    save_csv(OUT / 'policy_summary.csv', summary)
    freq = decisions.groupby(['metric', 'k', 'budget', 'condition', 'policy', 'algorithm']).size().rename('workloads').reset_index()
    freq['fraction'] = freq.workloads / freq.groupby(['metric', 'k', 'budget', 'condition', 'policy']).workloads.transform('sum')
    save_csv(OUT / 'choice_frequencies.csv', freq)
    return summary


def precision_profiles(decisions, native):
    frames = []
    for metric, profiles in CONFIG['accuracy_profiles'].items():
        d = decisions[decisions.metric == metric].copy()
        for i, (mae, p95) in enumerate(profiles):
            frame = d.copy()
            alternatives=native[native.metric==metric].copy()
            alternatives['any_measured_feasible']=((alternatives.availability>=.99)&
                (alternatives.budget_violation<=.01)&(alternatives.mae<=mae)&(alternatives.p95<=p95))
            possible=alternatives.groupby(KEYS).any_measured_feasible.any().reset_index()
            frame=frame.merge(possible,on=KEYS,how='left',validate='many_to_one')
            if frame.any_measured_feasible.isna().any():raise ValueError('missing feasible-alternative audit')
            frame['profile'] = 'standard' if i == 0 else 'strict'
            frame['mae_target'] = mae; frame['p95_target'] = p95
            frame['predicted_feasible'] = (frame.supported & (frame.predicted_availability >= .99) &
                (frame.predicted_budget_violation <= .01) & (frame.predicted_mae <= mae) & (frame.predicted_p95 <= p95))
            frame['achieved'] = ((frame.availability >= .99) & (frame.budget_violation <= .01) &
                (frame.mae <= mae) & (frame.p95 <= p95))
            frame['false_feasible'] = frame.predicted_feasible & ~frame.achieved
            frame['missed_feasible'] = ~frame.predicted_feasible & frame.achieved
            frame['selection_missed_target'] = frame.any_measured_feasible & ~frame.achieved
            frame['no_measured_feasible'] = ~frame.any_measured_feasible
            frames.append(frame)
    full = pd.concat(frames, ignore_index=True)
    save_csv(OUT / 'precision_profiles.csv.gz', full[KEYS + ['genus', 'condition', 'policy', 'profile',
        'mae_target', 'p95_target', 'predicted_feasible', 'achieved', 'false_feasible', 'missed_feasible',
        'any_measured_feasible','selection_missed_target','no_measured_feasible']])
    full = expanded_conditions(full)
    keys = ['metric', 'k', 'budget', 'condition', 'policy', 'profile', 'mae_target', 'p95_target']
    rows = full.groupby(keys).agg(workloads=('achieved', 'size'),
        recommendations=('predicted_feasible', 'sum'), achieved=('achieved', 'sum'),
        false_feasible=('false_feasible', 'sum'), missed_feasible=('missed_feasible', 'sum'),
        any_measured_feasible=('any_measured_feasible','sum'),
        selection_missed_target=('selection_missed_target','sum'),no_measured_feasible=('no_measured_feasible','sum')).reset_index()
    rows['recommendation_coverage'] = rows.recommendations / rows.workloads
    rows['attainment_rate'] = rows.achieved / rows.workloads
    rows['false_feasible_rate'] = rows.false_feasible / rows.recommendations.replace(0, np.nan)
    rows['feasible_alternative_rate'] = rows.any_measured_feasible / rows.workloads
    rows['selection_miss_rate'] = rows.selection_missed_target / rows.any_measured_feasible.replace(0,np.nan)
    rows['no_measured_feasible_rate'] = rows.no_measured_feasible / rows.workloads
    save_csv(OUT / 'precision_profile_summary.csv', rows)
    return rows


def bootstrap_weights(n):
    return np.random.default_rng(CONFIG['bootstrap_seed']).multinomial(n, np.repeat(1/n, n),
        size=CONFIG['bootstrap_replicates']) / n


def comparisons(decisions):
    rows = []
    d = expanded_conditions(decisions[decisions.policy.isin(['fixed', 'pilot', 'metadata'])])
    grouped = d.groupby(['metric', 'k', 'budget', 'condition', 'genus', 'policy']).loss.mean().unstack('policy')
    for (metric, k, budget, condition), frame in grouped.groupby(level=[0,1,2,3]):
        frame = frame.droplevel([0,1,2,3]); weights = bootstrap_weights(len(frame))
        for policy in ['pilot', 'metadata']:
            diff = (frame[policy] - frame.fixed).to_numpy()
            lo, hi = np.quantile(weights @ diff, [.025, .975])
            rows.append(dict(metric=metric, k=k, budget=budget, condition=condition, policy=policy,
                genera=len(frame), fixed_loss=float(frame.fixed.mean()), policy_loss=float(frame[policy].mean()),
                difference=float(diff.mean()), pointwise_low=float(lo), pointwise_high=float(hi),
                se=float(diff.std(ddof=1) / np.sqrt(len(diff)))))
    result = pd.DataFrame(rows)
    # Primary family: all 72 objective/k/budget contrasts at condition='all'.
    for policy in ['pilot', 'metadata']:
        sub = grouped.xs('all', level='condition')
        matrix = (sub[policy] - sub.fixed).unstack(['metric', 'k', 'budget']).sort_index()
        values = matrix.to_numpy(); means = values.mean(axis=0)
        ses = values.std(axis=0, ddof=1) / np.sqrt(len(matrix)); active = ses > 1e-12
        bootstrap = bootstrap_weights(len(matrix)) @ values
        critical = float(np.quantile(np.max(abs((bootstrap[:, active] - means[active]) / ses[active]), axis=1), .95)) if active.any() else 0.
        for j, (metric, k, budget) in enumerate(matrix.columns):
            mask = ((result.policy == policy) & (result.condition == 'all') & (result.metric == metric) &
                    (result.k == k) & (result.budget == budget))
            result.loc[mask, 'simultaneous_low'] = means[j] - critical * ses[j]
            result.loc[mask, 'simultaneous_high'] = means[j] + critical * ses[j]
            result.loc[mask, 'family_size'] = len(matrix.columns)
    result['material_improvement'] = ((result.simultaneous_high < 0) & (result.se > 1e-12) &
        (-result.difference >= np.maximum(.0001, .1 * result.fixed_loss)))
    save_csv(OUT / 'advisor_vs_fixed_intervals.csv', result)
    return result


def candidate_map(decisions, profile_summary=None):
    """Descriptive held-out candidates; this table never feeds the frozen advisor."""
    d = decisions[decisions.policy.isin(ALGORITHMS)]
    rows, contrasts = [], []
    keys = ['metric', 'k', 'budget', 'condition']
    for key_values, group in d.groupby(keys):
        base = dict(zip(keys, key_values))
        performance = group.groupby('policy').agg(loss=('loss','mean'), availability=('availability','mean'),
            budget_violation=('budget_violation','mean'), mae=('mae','mean'))
        eligible = sorted(performance[(performance.availability >= .99) & (performance.budget_violation <= .01)].index)
        if not eligible:
            rows.append(dict(base, candidates='', best_observed='', best_loss=np.nan, eligible=0,
                             interpretation='no backend meets availability/budget criterion')); continue
        pivot = group[group.policy.isin(eligible)].pivot(index='genus', columns='policy', values='loss')[eligible]
        if pivot.isna().any().any(): raise ValueError('unpaired biological validation cohorts')
        pairs = list(itertools.combinations(eligible, 2))
        beaten = set()
        if pairs:
            differences = np.array([(pivot[a] - pivot[b]).to_numpy() for a,b in pairs]).T
            means = differences.mean(axis=0); ses = differences.std(axis=0, ddof=1)/np.sqrt(len(pivot))
            active = ses > 1e-12; boot = bootstrap_weights(len(pivot)) @ differences
            critical = float(np.quantile(np.max(abs((boot[:,active]-means[active])/ses[active]),axis=1),.95)) if active.any() else 0.
            for i,(a,b) in enumerate(pairs):
                lo, hi = means[i]-critical*ses[i], means[i]+critical*ses[i]
                la, lb = float(performance.loc[a,'loss']), float(performance.loc[b,'loss'])
                a_wins = bool(ses[i] > 1e-12 and hi < 0 and lb-la >= max(.0001,.1*lb))
                b_wins = bool(ses[i] > 1e-12 and lo > 0 and la-lb >= max(.0001,.1*la))
                if a_wins: beaten.add(b)
                if b_wins: beaten.add(a)
                contrasts.append(dict(base,left=a,right=b,difference=means[i],simultaneous_low=lo,
                    simultaneous_high=hi,eligible_contrasts=len(pairs),left_supported=a_wins,right_supported=b_wins))
        remaining = [a for a in eligible if a not in beaten]
        best = min(eligible, key=lambda a:(performance.loc[a,'loss'],a))
        rows.append(dict(base,candidates=';'.join(remaining),best_observed=best,
            best_loss=float(performance.loc[best,'loss']),eligible=len(eligible),candidate_count=len(remaining),
            interpretation='held-out descriptive candidates; within-cell simultaneous contrasts, not global guarantees'))
    result=pd.DataFrame(rows)
    if profile_summary is not None:
        profile=profile_summary[(profile_summary.policy=='pilot')&(profile_summary.profile=='standard')]
        result=result.merge(profile[keys+['feasible_alternative_rate','no_measured_feasible_rate']],
                            on=keys,how='left',validate='one_to_one')
        if result.feasible_alternative_rate.isna().any():raise ValueError('candidate map lacks target-feasibility context')
        result['standard_target_status']=np.where(result.feasible_alternative_rate==0,'no_measured_configuration_meets_target',
            np.where(result.feasible_alternative_rate==1,'some_configuration_meets_target_in_all_workloads','mixed_feasibility'))
    save_csv(OUT / 'scenario_candidates.csv', result)
    save_csv(OUT / 'algorithm_contrasts.csv', pd.DataFrame(contrasts))
    return result


def task_summary(decisions):
    frames = [pd.read_csv(p) for p in sorted((OUT / 'validation/tasks').glob('*.csv.gz'))]
    detail = pd.concat(frames, ignore_index=True)
    selected = decisions[decisions.metric.isin(['neighbor','species','clustering'])].copy()
    take = [c for c in detail if c not in selected.columns or c in ['workload_id','k','budget','config']]
    d = selected.merge(detail[take],on=['workload_id','k','budget','config'],how='left',validate='many_to_one')
    for c in ['known_acceptance','unknown_false_acceptance','open_availability','unknown_availability']:
        d[c] = d[c].fillna(0.)
    save_csv(OUT / 'selected_task_outcomes.csv.gz', d)
    d = expanded_conditions(d)
    summary = d.groupby(['metric','k','budget','condition','policy']).agg(
        loss=('loss','mean'),known_acceptance=('known_acceptance','mean'),
        unknown_false_acceptance=('unknown_false_acceptance','mean'),open_availability=('open_availability','mean'),
        unknown_availability=('unknown_availability','mean'),conditional_ari=('ari','mean'),
        exact_species_accuracy=('exact_species_accuracy','mean')).reset_index()
    save_csv(OUT / 'downstream_summary.csv', summary)
    exact = expanded_conditions(detail[detail.algorithm == 'exact'])
    exact = exact.groupby(['k','budget','condition']).agg(known_acceptance=('known_acceptance','mean'),
        unknown_false_acceptance=('unknown_false_acceptance','mean'),ari=('ari','mean'),
        species_accuracy=('exact_species_accuracy','mean')).reset_index()
    save_csv(OUT / 'exact_downstream_summary.csv', exact)
    return summary


def main():
    identity = verify_freezes()
    decisions, native, truth = match_decisions()
    summarize(decisions); profiles=precision_profiles(decisions,native); intervals = comparisons(decisions)
    candidate_map(decisions,profiles); task_summary(decisions)
    identity.update(evaluate_sha256=sha(__file__), validation_native_score_rows=len(native),
        frozen_decision_rows=len(decisions), validation_genera=int(decisions.genus.nunique()),
        validation_workloads=int(decisions.workload_id.nunique()),
        primary_material_improvements=int(intervals[(intervals.policy=='pilot') & (intervals.condition=='all')].material_improvement.sum()))
    atomic(OUT / 'evaluation_complete.json', identity)
    print(json.dumps(identity,indent=2),flush=True)


if __name__ == '__main__': main()
