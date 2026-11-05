#!/usr/bin/env python3
"""Post-registration constant-baseline audit and illustrative error requirements."""
from common import *

PROFILES = {
    'j_error': [('absolute', x) for x in [.001, .005, .01, .05]] + [('relative', x) for x in [.1, .2, .5]],
    'c_error': [('absolute', x) for x in [.01, .05, .1]] + [('relative', x) for x in [.2, .5]],
    'card_relative_error': [('absolute', x) for x in [.01, .02, .05]],
}


def qualifies(row, limit):
    return (np.isfinite(limit) and limit > 0 and row['availability'] >= .99 and
        row['budget_violation_fraction'] <= .01 and row['mae'] <= limit and row['p95'] <= 2 * limit)


def main():
    s = pd.read_csv(OUT / 'synthetic_accuracy.csv')
    decisions = pd.read_csv(OUT / 'synthetic_decisions.csv')
    validation = s[s.split == 'validation'].copy()
    diagnostic = validation[validation.metric.isin(['j_error', 'c_error'])].copy()
    diagnostic['truth'] = np.where(diagnostic.metric == 'j_error', diagnostic.truth_j, diagnostic.truth_c)
    diagnostic['constant_zero_mae'] = diagnostic.truth
    diagnostic['constant_one_mae'] = 1 - diagnostic.truth
    diagnostic['relative_to_zero_mae'] = diagnostic.mae / diagnostic.truth.replace(0, np.nan)
    diagnostic['beats_constant_zero'] = diagnostic.mae < diagnostic.constant_zero_mae - 1e-14
    diagnostic['beats_constant_one'] = diagnostic.mae < diagnostic.constant_one_mae - 1e-14
    save('synthetic_constant_baselines', diagnostic)
    winner = decisions[decisions.supported_winner.notna()].merge(
        diagnostic[['regime_id', 'budget', 'metric', 'algorithm', 'zero_fraction', 'relative_to_zero_mae',
                    'beats_constant_zero', 'beats_constant_one']],
        left_on=['regime_id', 'budget', 'metric', 'supported_winner'],
        right_on=['regime_id', 'budget', 'metric', 'algorithm'], how='left', validate='one_to_one')
    save('synthetic_winner_diagnostics', winner)
    rows = []
    for (regime, budget, metric), g in s.groupby(['regime_id', 'budget', 'metric'], sort=True):
        meta = g.iloc[0]
        cal = g[g.split == 'calibration'].set_index('algorithm')
        val = g[g.split == 'validation'].set_index('algorithm')
        truth = meta.truth_j if metric == 'j_error' else meta.truth_c
        for scale, target in PROFILES[metric]:
            limit = target * truth if scale == 'relative' else target
            inapplicable = scale == 'relative' and truth == 0
            c = [] if inapplicable else [a for a, r in cal.iterrows() if qualifies(r, limit)]
            v = [] if inapplicable else [a for a, r in val.iterrows() if qualifies(r, limit)]
            nominee = min(c, key=lambda a: (cal.loc[a, 'mae'], a)) if c else ''
            rows.append(dict(regime_id=regime, k=int(meta.k), n_b=int(meta.n_b), size_ratio=int(meta.size_ratio),
                nominal_containment=meta.nominal_containment, truth_j=meta.truth_j, truth_c=meta.truth_c,
                budget=budget, metric=metric, scale=scale, target=target,
                absolute_mae_limit=limit, status='inapplicable_zero_truth' if inapplicable else
                    'no_calibration_candidate' if not nominee else 'pass' if nominee in v else 'validation_failure',
                calibration_candidates=';'.join(sorted(c)), validation_candidates=';'.join(sorted(v)),
                calibration_nominee=nominee, nominee_validates=int(bool(nominee) and nominee in v),
                validation_mae=val.loc[nominee, 'mae'] if nominee else np.nan,
                validation_p95=val.loc[nominee, 'p95'] if nominee else np.nan))
    profiles = pd.DataFrame(rows)
    save('synthetic_constraint_profiles', profiles)
    save('synthetic_constraint_coverage', profiles.groupby(['metric', 'scale', 'target', 'budget', 'status']).size().rename('regimes').reset_index())
    minimum = profiles[profiles.calibration_nominee != ''].sort_values('budget').drop_duplicates(
        ['regime_id', 'metric', 'scale', 'target'])
    # Retain the first calibration-feasible budget even when validation fails.
    save('synthetic_minimum_calibrated_budget', minimum)
    configurations = pd.read_csv(OUT / 'synthetic_all_configurations.csv')
    frontier = configurations[(configurations.split == 'validation') &
        (configurations.availability >= .99) & (configurations.max_pair_bytes <= max(CONFIG['pair_budgets']))].copy()
    frontier['descriptive_frontier'] = False
    for _, group in frontier.groupby(['regime_id', 'metric']):
        costs = group.max_pair_bytes.to_numpy(); errors = group.mae.to_numpy()
        dominates = (costs[:, None] <= costs[None, :]) & (errors[:, None] <= errors[None, :]) & (
            (costs[:, None] < costs[None, :]) | (errors[:, None] < errors[None, :]))
        frontier.loc[group.index, 'descriptive_frontier'] = ~dominates.any(axis=0)
    save('synthetic_error_payload_frontier', frontier)
    atomic(OUT / 'constraint_profiles_complete.json', dict(status='complete', post_registration_diagnostic=True,
        profiles=len(profiles), source_sha256=sha(__file__),
        addendum_sha256=sha(HERE / 'ANALYSIS_ADDENDUM.md')))


if __name__ == '__main__': main()
