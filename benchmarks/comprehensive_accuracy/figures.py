#!/usr/bin/env python3
"""Show effect size, availability and uncertainty separately."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm
from matplotlib.patches import Patch

from common import *

plt.rcParams.update({'pdf.fonttype': 42, 'axes.spines.top': False,
    'axes.spines.right': False, 'font.size': 10})


def finish(fig, name, rect=None):
    fig.tight_layout(rect=rect)
    for extension in ['png', 'pdf']:
        fig.savefig(OUT / (name + '.' + extension), dpi=180)
    plt.close(fig)


def error_surfaces(s):
    ratios, cs = CONFIG['size_ratios'], CONFIG['small_set_containments']
    for metric, title in [('j_error', 'Absolute Jaccard error'), ('c_error', 'Directed containment error')]:
        d = s[(s.k == 24) & (s.n_b == 1000000) & (s.budget == 8192) & (s.metric == metric) & (s.split == 'validation')]
        algorithms = ALGORITHMS if metric == 'j_error' else [a for a in ALGORITHMS if a not in ['kssd', 'probminhash']]
        fig, axes = plt.subplots(2, 4 if metric == 'j_error' else 3, figsize=(15, 8), squeeze=False)
        image = None
        for ax, algorithm in zip(axes.flat, algorithms):
            g = d[d.algorithm == algorithm]
            z = g.pivot(index='nominal_containment', columns='size_ratio', values='mae').reindex(index=cs, columns=ratios)
            image = ax.imshow(np.maximum(z.to_numpy(), 1e-6), aspect='auto', origin='lower',
                              cmap='viridis', norm=LogNorm(vmin=1e-5, vmax=.5))
            ax.set_title(NAMES[algorithm]); ax.set_xticks(range(len(ratios)), ratios)
            ax.set_yticks(range(len(cs)), [f'{c:g}' for c in cs]); ax.set_xlabel('Larger / smaller cardinality')
            ax.set_ylabel('Nominal small-set containment')
            for row in g.itertuples():
                if row.availability < .99 or row.budget_violation_fraction > .01:
                    ax.text(ratios.index(row.size_ratio), cs.index(row.nominal_containment), '×',
                            ha='center', va='center', color='red', fontsize=13, fontweight='bold')
        for ax in list(axes.flat)[len(algorithms):]: ax.axis('off')
        fig.suptitle(f'{title}: validation MAE among available outputs\nk=24, larger set=1,000,000, 8 KiB/pair; × fails availability or budget', fontsize=12)
        fig.subplots_adjust(right=.9, top=.86, hspace=.55, wspace=.55)
        color_ax = fig.add_axes([.92, .2, .015, .55]); fig.colorbar(image, cax=color_ax, label='MAE (log scale; zero clipped for display)')
        for ext in ['png', 'pdf']: fig.savefig(OUT / (metric + '_surfaces.' + ext), dpi=180)
        plt.close(fig)


def candidate_maps(decisions):
    ratios, cs = CONFIG['size_ratios'], CONFIG['small_set_containments']
    colors = [COLORS[a] for a in ALGORITHMS] + ['#D5D7DB', '#222222']
    cmap = ListedColormap(colors); mapping = {a: i for i, a in enumerate(ALGORITHMS)}
    for metric, title in [('j_error', 'Jaccard'), ('c_error', 'Directed containment')]:
        fig, axes = plt.subplots(3, 4, figsize=(15, 10), sharex=True, sharey=True)
        for i, nb in enumerate(CONFIG['larger_cardinalities']):
            for j, budget in enumerate(CONFIG['pair_budgets']):
                ax = axes[i, j]
                d = decisions[(decisions.k == 24) & (decisions.n_b == nb) & (decisions.budget == budget) & (decisions.metric == metric)].copy()
                d['state'] = d.supported_winner.map(mapping).fillna(len(ALGORITHMS))
                d.loc[d.status == 'no_candidate', 'state'] = len(ALGORITHMS) + 1
                matrix = d.pivot(index='nominal_containment', columns='size_ratio', values='state').reindex(index=cs, columns=ratios)
                ax.imshow(matrix.to_numpy(), aspect='auto', origin='lower', cmap=cmap, vmin=0, vmax=len(colors)-1)
                for row in d.itertuples():
                    if row.status != 'supported_winner':
                        count = len(str(row.candidate_set).split(';')) if pd.notna(row.candidate_set) else 0
                        ax.text(ratios.index(row.size_ratio), cs.index(row.nominal_containment), str(count),
                                ha='center', va='center', fontsize=7, color='white' if row.status == 'no_candidate' else '#555555')
                ax.set_title(f'n_large={nb:,}; {budget/1024:g} KiB/pair', fontsize=10)
                ax.set_xticks(range(len(ratios)), ratios, rotation=45)
                ax.set_yticks(range(len(cs)), [f'{c:g}' for c in cs])
                if i == 2: ax.set_xlabel('Larger / smaller cardinality')
                if j == 0: ax.set_ylabel('Nominal small-set containment')
        handles = [Patch(facecolor=COLORS[a], label=NAMES[a]) for a in ALGORITHMS]
        handles += [Patch(facecolor='#D5D7DB', label='Unresolved: number of surviving candidates'), Patch(facecolor='#222222', label='No eligible candidate')]
        fig.legend(handles=handles, loc='lower center', ncol=3, fontsize=9)
        fig.suptitle(f'{title}: supported winner among eligible candidates, k=24\nPaired simultaneous 95% bands within regime; >=10% and >=0.0001 MAE improvement', fontsize=12)
        finish(fig, metric + '_candidate_map', rect=[0, .1, 1, .94])


def relative_precision(s):
    d = s[(s.k == 24) & (s.n_b == 1000000) & (s.split == 'validation') &
          (s.metric == 'j_error') & (s.nominal_containment == .5)]
    fig, axes = plt.subplots(1, 4, figsize=(15, 4), sharex=True, sharey=True)
    for ax, budget in zip(axes, CONFIG['pair_budgets']):
        for algorithm, g in d[d.budget == budget].groupby('algorithm'):
            g = g.sort_values('size_ratio')
            good = (g.availability >= .99) & (g.budget_violation_fraction <= .01)
            ax.plot(g.size_ratio, g.relative_mae.where(good), 'o-', ms=4, color=COLORS[algorithm], label=NAMES[algorithm])
            ax.scatter(g.loc[~good, 'size_ratio'], g.loc[~good, 'relative_mae'], marker='x', s=35, color=COLORS[algorithm])
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_title(f'{budget/1024:g} KiB/pair')
        ax.set_xlabel('Larger / smaller cardinality'); ax.grid(alpha=.2)
    axes[0].set_ylabel('MAE / exact Jaccard'); axes[-1].legend(fontsize=7)
    fig.suptitle('Small absolute J error can still mean poor relative accuracy\nk=24, n_large=1,000,000, small-set containment=0.5; × fails budget/availability')
    finish(fig, 'relative_precision')


def real_figures():
    s = pd.read_csv(OUT / 'real_budget_accuracy.csv')
    d = s[(s.phase == 'natural') & (s.k == 24) & (s.split == 'validation') &
          (s.metric == 'j_error') & s.stratum.isin(['J_0_001', 'J_001_01', 'J_01_05', 'J_05_1'])]
    strata = ['J_0_001', 'J_001_01', 'J_01_05', 'J_05_1']
    fig, axes = plt.subplots(1, 4, figsize=(15, 4), sharey=True)
    for ax, budget in zip(axes, CONFIG['pair_budgets']):
        for algorithm, g in d[d.budget == budget].groupby('algorithm'):
            g = g.set_index('stratum').reindex(strata)
            good = (g.availability >= .99) & (g.budget_violation_fraction <= .01)
            ax.plot(range(4), g.mae.where(good), 'o-', color=COLORS[algorithm], label=NAMES[algorithm])
        ax.set_title(f'{budget/1024:g} KiB/pair'); ax.set_xticks(range(4), ['[0,.01)', '[.01,.1)', '[.1,.5)', '[.5,1]'], rotation=45)
        ax.set_yscale('log'); ax.set_xlabel('Exact Jaccard stratum'); ax.grid(alpha=.2)
    axes[0].set_ylabel('Validation MAE'); axes[-1].legend(fontsize=7)
    fig.suptitle('Real genomes: held-out genera, k=24; parameters fixed on calibration genera\nGaps: no eligible configuration/output in that stratum; same 900-genome inventory')
    finish(fig, 'real_genome_accuracy')

    downstream = pd.read_csv(OUT / 'downstream_summary.csv')
    d = downstream[(downstream.k == 24) & (downstream.split == 'validation')]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, metric, title in zip(axes, ['species_accuracy', 'exact_neighbor_recovery', 'ari'],
                                 ['Species nearest neighbor', 'Exact-J neighbor recovery', 'Single-linkage ARI']):
        for algorithm, g in d[d.metric == metric].groupby('algorithm'):
            g = g.sort_values('budget'); good = (g.valid_genera == 30) & (g.calibration_complete == 1)
            values = g['mean'].where(good)
            ax.plot(g.budget / 1024, values, 'o-' if algorithm != 'exact' else 'k--',
                    **({'color': COLORS[algorithm]} if algorithm != 'exact' else {}), label=NAMES[algorithm])
        ax.set_xscale('log', base=2); ax.set_title(title); ax.set_xlabel('Pair payload budget (KiB)'); ax.grid(alpha=.2)
    axes[0].set_ylabel('Mean across 30 validation genera'); axes[-1].legend(fontsize=7)
    fig.suptitle('Downstream validation at k=24: no threshold refitting on validation genera')
    finish(fig, 'downstream_accuracy')

    choices = pd.read_csv(OUT / 'downstream_calibrated_choices.csv')
    choices = choices[(choices.k == 24) & (choices.budget == 8192)].set_index('metric')
    pairs = pd.read_csv(OUT / 'downstream_paired_intervals.csv')
    pairs = pairs[(pairs.k == 24) & (pairs.budget == 8192) & (pairs.common_genera == 30)]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric, title in zip(axes, ['species_accuracy', 'exact_neighbor_recovery', 'ari'],
                                ['Species nearest neighbor', 'Exact-J neighbor recovery', 'Single-linkage ARI']):
        nominee = choices.loc[metric, 'algorithm']; alternatives = [a for a in ALGORITHMS if a != nominee]
        factor = 1 if metric == 'ari' else 100
        for y, other in enumerate(alternatives):
            row = pairs[(pairs.metric == metric) & (((pairs.left == nominee) & (pairs.right == other)) |
                        ((pairs.right == nominee) & (pairs.left == other)))].iloc[0]
            mean, lo, hi = row.difference, row.ci_low, row.ci_high
            if row.right == nominee: mean, lo, hi = -mean, -hi, -lo
            ax.errorbar(mean * factor, y, xerr=np.array([[max(0, mean-lo)], [max(0, hi-mean)]]) * factor,
                        fmt='o', color=COLORS[other], capsize=3)
        ax.axvline(0, color='#555555', linestyle='--'); ax.grid(axis='x', alpha=.2)
        ax.set_yticks(range(len(alternatives)), [NAMES[a] for a in alternatives])
        ax.set_title(f'{title}\nCalibration choice: {NAMES[nominee]}', fontsize=10)
        ax.set_xlabel('Choice minus alternative (ARI)' if metric == 'ari' else 'Choice minus alternative (percentage points)')
    fig.suptitle('Do fixed calibration choices improve held-out tasks? k=24, 8 KiB/pair\nPaired 95% genus-bootstrap intervals; 30 validation genera; pointwise, no multiplicity correction')
    finish(fig, 'downstream_paired_effects')


def payload_frontier():
    path = OUT / 'synthetic_error_payload_frontier.csv'
    if not path.exists(): return
    d = pd.read_csv(path)
    d = d[(d.k == 24) & (d.n_b == 1000000) & (d.metric == 'j_error')]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=False)
    for ax, (ratio, c) in zip(axes, [(1, .1), (1, .9), (100, .5)]):
        g = d[(d.size_ratio == ratio) & (d.nominal_containment == c)]
        for alg, a in g.groupby('algorithm'):
            ax.scatter(a.max_pair_bytes / 1024, a.mae, color=COLORS[alg], label=NAMES[alg], s=30)
        front = g[g.descriptive_frontier].sort_values('max_pair_bytes')
        ax.plot(front.max_pair_bytes / 1024, front.mae, 'k--', alpha=.6, label='Observed frontier')
        if len(g): ax.axhline(g.truth_j.iloc[0], color='#888888', linestyle=':', label='Constant-zero error')
        ax.set_xscale('log', base=2); ax.set_yscale('log'); ax.set_xlabel('Maximum observed pair payload (KiB)')
        ax.set_title(f'Cardinality ratio={ratio}; containment={c:g}'); ax.grid(alpha=.2)
    axes[0].set_ylabel('Validation Jaccard MAE'); axes[-1].legend(fontsize=7)
    fig.suptitle('Descriptive error/payload frontier: measured configurations at k=24, n_large=1,000,000\nComplete outputs, at most 32 KiB; not a significance test or online selector')
    finish(fig, 'error_payload_frontier')


def controlled_figures():
    s = pd.read_csv(OUT / 'real_budget_accuracy.csv')
    d = s[(s.phase == 'containment') & (s.split == 'validation') & (s.metric == 'c_error')].copy()
    extracted = d.stratum.str.extract(r'^ratio([\d.]+)_containment([\d.]+)$').astype(float)
    d['ratio'], d['nominal_c'] = extracted[0], extracted[1]
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=True, sharey='row')
    for i, containment in enumerate([.5, 1.]):
        for ax, budget in zip(axes[i], CONFIG['pair_budgets']):
            selected = d[(d.budget == budget) & (d.nominal_c == containment)]
            for algorithm, g in selected.groupby('algorithm'):
                if algorithm in ['kssd', 'probminhash']: continue
                g = g.sort_values('ratio')
                good = (g.availability >= .99) & (g.budget_violation_fraction <= .01)
                displayed = g.mae.clip(lower=1e-5)
                ax.plot(g.ratio, displayed.where(good), 'o-', color=COLORS[algorithm], label=NAMES[algorithm])
                ax.scatter(g.loc[~good, 'ratio'], displayed.loc[~good], marker='x', color=COLORS[algorithm])
            ax.set_title(f'Nominal C={containment:g}; {budget/1024:g} KiB/pair')
            ax.set_xscale('log'); ax.set_yscale('log'); ax.grid(alpha=.2)
            if i == 1: ax.set_xlabel('Reference / query sequence length')
        axes[i, 0].set_ylabel('Validation containment MAE')
    axes[0, -1].legend(fontsize=7)
    fig.suptitle('Controlled containment on real sequences, k=24: fixed calibration parameters\nExact k-mer truth; errors below 1e-5 shown at 1e-5; gaps / crosses fail output or budget constraints')
    finish(fig, 'real_containment_accuracy')

    d = s[(s.phase == 'mixed') & (s.split == 'validation') & (s.budget == 8192)].copy()
    extracted = d.stratum.str.extract(r'^identity([\d.]+)_indel([\d.]+)$').astype(float)
    d['identity'], d['indel'] = extracted[0], extracted[1]
    truth = pd.read_csv(OUT / 'real_panels/mixed_k24/exact_truth.csv')
    truth = truth[truth.split == 'validation'].copy()
    truth['model_error'] = abs((2 * truth.truth_j / (1 + truth.truth_j)) ** (1 / 24) - truth.column_identity)
    model = truth.groupby(['target_identity', 'indel_share']).model_error.mean()
    fig, axes = plt.subplots(3, 4, figsize=(15, 10), sharex=True, sharey='row')
    metrics = [('j_error', 'Jaccard sampling MAE'), ('ani_sampling_error', 'ANI sampling MAE'),
               ('column_ani_error', 'Column-identity total MAE')]
    for i, (metric, label) in enumerate(metrics):
        for ax, identity in zip(axes[i], sorted(d.identity.dropna().unique())):
            selected = d[(d.metric == metric) & (d.identity == identity)]
            for algorithm, g in selected.groupby('algorithm'):
                g = g.sort_values('indel')
                good = (g.availability >= .99) & (g.budget_violation_fraction <= .01)
                ax.plot(g.indel, g.mae.where(good), 'o-', color=COLORS[algorithm], label=NAMES[algorithm])
            if i == 2:
                base = model.loc[identity]
                ax.plot(base.index, base.to_numpy(), 'k--', label='Exact-J model residual')
            ax.set_title(f'Target identity={identity:g}'); ax.set_yscale('log'); ax.grid(alpha=.2)
            if i == 2: ax.set_xlabel('Indel share of edits')
        axes[i, 0].set_ylabel(label)
    axes[-1, -1].legend(fontsize=7)
    fig.suptitle('Substitution / indel validation, k=24, 8 KiB/pair\nSampling error and identity-model error are separate outcomes; validation references only')
    finish(fig, 'mixed_mutation_accuracy')

    weighted = pd.read_csv(OUT / 'weighted_accuracy.csv')
    weighted = weighted[weighted.split == 'validation'].copy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    cami = weighted[weighted.family.str.contains('cami', case=False)].sort_values('parameter')
    axes[0].plot(cami.pair_payload / 1024, cami.mae, 'o-', label='MAE')
    axes[0].plot(cami.pair_payload / 1024, cami.p95, 's--', label='P95')
    axes[0].set_title('CAMI marine profiles'); axes[0].legend()
    synthetic = weighted[weighted.family == 'synthetic_weights'].copy()
    synthetic['error_sum'] = synthetic.mae * synthetic.observations
    pooled = synthetic.groupby(['model', 'pair_payload'])[['error_sum', 'observations']].sum().reset_index()
    pooled['mae'] = pooled.error_sum / pooled.observations
    for model, g in pooled.groupby('model'):
        axes[1].plot(g.pair_payload / 1024, g.mae, 'o-', label=model)
    axes[1].set_title('Synthetic abundance models'); axes[1].legend(fontsize=7)
    for ax in axes:
        ax.set_xscale('log', base=2); ax.set_yscale('log'); ax.grid(alpha=.2)
        ax.set_xlabel('Pair payload (KiB)'); ax.set_ylabel('Probability-Jaccard error')
    fig.suptitle('Weighted ProbMinHash: a separate target, eight hash seeds\nNo competing weighted backend was tested; not min/max weighted Jaccard')
    finish(fig, 'weighted_accuracy')


def requirement_map():
    path = OUT / 'synthetic_constraint_profiles.csv'
    if not path.exists(): return
    d = pd.read_csv(path); d = d[(d.k == 24) & (d.n_b == 1000000)]
    mapping = {a: i for i, a in enumerate(ALGORITHMS)}
    colors = [COLORS[a] for a in ALGORITHMS] + ['#222222', '#E58C8C', '#FFFFFF']
    ratios, cs = CONFIG['size_ratios'], CONFIG['small_set_containments']
    fig, axes = plt.subplots(2, 4, figsize=(15, 8), sharex=True, sharey=True)
    requirements = [('j_error', 'relative', .2, 'J: relative MAE <=20%'),
                    ('c_error', 'absolute', .05, 'Containment: MAE <=0.05')]
    for i, (metric, scale, target, title) in enumerate(requirements):
        for j, budget in enumerate(CONFIG['pair_budgets']):
            g = d[(d.metric == metric) & (d.scale == scale) & (d.target == target) & (d.budget == budget)].copy()
            g['state'] = g.calibration_nominee.map(mapping)
            for status, code in [('no_calibration_candidate', 7), ('validation_failure', 8), ('inapplicable_zero_truth', 9)]:
                g.loc[g.status == status, 'state'] = code
            values = g.pivot(index='nominal_containment', columns='size_ratio', values='state').reindex(index=cs, columns=ratios)
            ax = axes[i, j]; ax.imshow(values, aspect='auto', origin='lower', cmap=ListedColormap(colors), vmin=0, vmax=9)
            ax.set_title(f'{title}\n{budget/1024:g} KiB/pair', fontsize=10)
            ax.set_xticks(range(len(ratios)), ratios, rotation=45); ax.set_yticks(range(len(cs)), [f'{c:g}' for c in cs])
            if i == 1: ax.set_xlabel('Larger / smaller cardinality')
            if j == 0: ax.set_ylabel('Nominal small-set containment')
    labels = [NAMES[a] for a in ALGORITHMS] + ['No calibration-feasible candidate', 'Fixed nominee fails validation', 'Relative error undefined at J=0']
    fig.legend(handles=[Patch(facecolor=c, edgecolor='#777777', label=l) for c, l in zip(colors, labels)],
               loc='lower center', ncol=3, fontsize=9)
    fig.suptitle('Illustrative error requirements: calibration nominee with held-out empirical validation\nP95 <=2x MAE limit; output >=99%; budget violation <=1%; post-registration diagnostic', fontsize=12)
    finish(fig, 'error_requirement_map', rect=[0, .13, 1, .92])


def main():
    s = pd.read_csv(OUT / 'synthetic_accuracy.csv')
    decisions = pd.read_csv(OUT / 'synthetic_decisions.csv')
    error_surfaces(s); candidate_maps(decisions); relative_precision(s); real_figures(); payload_frontier(); requirement_map(); controlled_figures()
    print('Accuracy, candidate and downstream figures written', flush=True)


if __name__ == '__main__': main()
