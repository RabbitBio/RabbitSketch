#!/usr/bin/env python3
"""Read a completed pilot and write exploratory, budget-matched diagnostics.

The registered pruning analysis is preserved. This follow-up selects the
largest measured payload fitting each budget, without selecting on accuracy.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from run import DEFAULT_ROOT, HERE, atomic_json, digest, file_sha

ALGORITHMS = ['minhash', 'fracminhash', 'hll', 'setsketch', 'kssd', 'probminhash', 'bindash']
LABELS = dict(zip(ALGORITHMS, ['MinHash', 'FracMinHash', 'HLL', 'SetSketch', 'KSSD',
                             'ProbMinHash (unweighted)', 'BinDash']))
BUDGETS = [1024 * n for n in (1, 2, 4, 8, 16, 32, 64, 128, 256)]
SHARED_BUDGETS = [2048, 8192, 32768]
RELATIVE_LIMIT = .25


def paired_interval(a, b, indices):
    differences = np.asarray(a) - np.asarray(b)
    return np.quantile(differences[indices].mean(axis=1), [.025, .975])


def budget_selections(configurations, budgets):
    """One configuration per algorithm/regime: use capacity, never observed error."""
    selections = []
    good = configurations[configurations.eligible]
    for budget in budgets:
        possible = good[good.max_observed_payload_bytes <= budget]
        selected = (possible.sort_values(['payload_bytes', 'config_id'], ascending=[False, True])
                    .drop_duplicates(['regime_id', 'algorithm']).copy())
        selected['budget_bytes'] = budget
        selected['budget_kib'] = budget / 1024
        selected['relative_mae'] = selected.mae / selected.ground_truth
        selected['relative_mae_ci_upper'] = selected.mae_ci_upper / selected.ground_truth
        selected['budget_used_fraction'] = selected.payload_bytes / budget
        selections.append(selected)
    return pd.concat(selections, ignore_index=True)


def overview(configurations, decisions):
    records = []
    for name in ALGORITHMS:
        g = configurations[configurations.algorithm == name]
        regimes = g.groupby('regime_id')[['pareto_frontier', 'error_space_frontier']].any()
        decision = decisions[decisions.algorithm == name].iloc[0]
        records.append(dict(algorithm=name, native_configurations=g.config_id.nunique(),
                            regimes=len(regimes), four_objective_frontier_regimes=int(regimes.pareto_frontier.sum()),
                            error_space_frontier_regimes=int(regimes.error_space_frontier.sum()),
                            robustly_dominated_regimes=int(decision.robustly_dominated_regimes),
                            empty_sample_pairs=int(g.empty_samples.sum()),
                            ineligible_configuration_regimes=int((~g.eligible).sum()),
                            minimum_mean_pair_payload_bytes=float(g.payload_bytes.min()),
                            maximum_mean_pair_payload_bytes=float(g.payload_bytes.max()),
                            maximum_shared_process_rss_bytes=int(g.process_peak_rss_bytes.max()),
                            registered_decision=decision.status))
    return pd.DataFrame(records)


def compare_budgets(selected, errors, bootstrap):
    comparisons = []
    for (budget, regime), g in selected.groupby(['budget_bytes', 'regime_id']):
        by_algorithm = g.set_index('algorithm')
        if 'minhash' not in by_algorithm.index:
            continue
        reference = by_algorithm.loc['minhash']
        for name in ('bindash', 'fracminhash', 'probminhash'):
            if name not in by_algorithm.index:
                continue
            challenger = by_algorithm.loc[name]
            lo, hi = paired_interval(errors[(regime, challenger.config_id)],
                                     errors[(regime, reference.config_id)], bootstrap[regime])
            comparisons.append(dict(budget_bytes=budget, budget_kib=budget/1024, regime_id=regime,
                                    k=int(reference.k), n_a=int(reference.n_a), n_b=int(reference.n_b),
                                    size_ratio=int(reference.size_ratio), ground_truth=reference.ground_truth,
                                    challenger=name, reference='minhash',
                                    challenger_parameter=int(challenger.parameter),
                                    reference_parameter=int(reference.parameter),
                                    challenger_max_payload=challenger.max_observed_payload_bytes,
                                    reference_max_payload=reference.max_observed_payload_bytes,
                                    challenger_mae=challenger.mae, reference_mae=reference.mae,
                                    mae_ratio=challenger.mae/reference.mae,
                                    paired_mae_difference_ci_lower=float(lo),
                                    paired_mae_difference_ci_upper=float(hi),
                                    challenger_build_ratio=challenger.build_seconds/reference.build_seconds,
                                    challenger_query_ratio=challenger.query_seconds/reference.query_seconds))
    return pd.DataFrame(comparisons)


def hll_setsketch(configurations, errors, card_errors, bootstrap):
    h = configurations[configurations.algorithm == 'hll']
    s = configurations[configurations.algorithm == 'setsketch']
    matched = h.merge(s, on=['regime_id', 'parameter'], suffixes=('_hll', '_ss'), validate='one_to_one')
    records = []
    for r in matched.itertuples():
        assert r.payload_bytes_hll == r.payload_bytes_ss
        a, b = (r.regime_id, r.config_id_ss), (r.regime_id, r.config_id_hll)
        jlo, jhi = paired_interval(errors[a], errors[b], bootstrap[r.regime_id])
        clo, chi = paired_interval(card_errors[a], card_errors[b], bootstrap[r.regime_id])
        records.append(dict(regime_id=r.regime_id, parameter=r.parameter,
                            k=r.k_hll, n_a=r.n_a_hll, n_b=r.n_b_hll, size_ratio=r.size_ratio_hll,
                            ground_truth=r.ground_truth_hll, pair_payload_bytes=r.payload_bytes_hll,
                            hll_jaccard_mae=r.mae_hll, setsketch_jaccard_mae=r.mae_ss,
                            jaccard_mae_ratio_ss_to_hll=r.mae_ss/r.mae_hll,
                            jaccard_difference_ci_lower=float(jlo), jaccard_difference_ci_upper=float(jhi),
                            hll_cardinality_relative_mae=float(card_errors[b].mean()),
                            setsketch_cardinality_relative_mae=float(card_errors[a].mean()),
                            cardinality_mae_ratio_ss_to_hll=float(card_errors[a].mean()/card_errors[b].mean()),
                            cardinality_difference_ci_lower=float(clo), cardinality_difference_ci_upper=float(chi)))
    return pd.DataFrame(records)


def relative_guarded_reuse(configurations, regimes):
    records = []
    for regime in regimes:
        good = configurations[(configurations.regime_id == regime['regime_id']) & configurations.eligible]
        limit = min(.005, RELATIVE_LIMIT * regime['ground_truth'])
        for budget in BUDGETS:
            possible = good[(good.max_observed_payload_bytes <= budget) & (good.mae_ci_upper <= limit)]
            for q in (1, 100, 10000):
                row = dict(regime_id=regime['regime_id'], k=regime['k'], n_b=regime['n_b'],
                           size_ratio=regime['size_ratio'], ground_truth=regime['ground_truth'],
                           budget_bytes=budget, budget_kib=budget/1024, query_reuse=q,
                           mae_upper_bound_limit=limit)
                if possible.empty:
                    row['algorithm'] = 'no_feasible_configuration'
                else:
                    ranked = possible.assign(predicted_seconds=possible.build_seconds + q*possible.query_seconds)
                    ranked = ranked.sort_values(['predicted_seconds', 'config_id'])
                    first = ranked.iloc[0]
                    row.update(algorithm=first.algorithm, parameter=int(first.parameter),
                               predicted_seconds=first.predicted_seconds,
                               within_10_percent_algorithms=';'.join(sorted(set(ranked.loc[
                                   ranked.predicted_seconds <= 1.1*first.predicted_seconds, 'algorithm']))))
                records.append(row)
    return pd.DataFrame(records)


def figures(output, config, selected, frontiers, hll_ss):
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/rabbitsketch_screening_mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    colors = dict(zip(ALGORITHMS, plt.get_cmap('tab10').colors))

    def save(fig, name):
        fig.savefig(output/(name+'.png'), dpi=170, bbox_inches='tight')
        fig.savefig(output/(name+'.pdf'), bbox_inches='tight')
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    x = np.arange(len(ALGORITHMS))
    axes[0].barh(x-.18, frontiers.four_objective_frontier_regimes, height=.35, label='Error, bytes, build, query')
    axes[0].barh(x+.18, frontiers.error_space_frontier_regimes, height=.35, label='Error and bytes')
    axes[0].set_yticks(x, [LABELS[n] for n in ALGORITHMS])
    axes[0].set_xlim(0, 89); axes[0].set_xlabel('Regimes with any frontier configuration (out of 84)')
    axes[0].legend(loc='lower right', fontsize=8); axes[0].grid(axis='x', alpha=.2)
    for name in ALGORITHMS:
        counts = (selected[(selected.algorithm == name) &
                           (selected.relative_mae_ci_upper <= RELATIVE_LIMIT)]
                  .groupby('budget_kib').size().reindex([b/1024 for b in BUDGETS], fill_value=0))
        axes[1].plot(counts.index, counts.values, '-o', label=LABELS[name], color=colors[name], ms=4)
    axes[1].set_xscale('log', base=2); axes[1].set_ylim(0, 87)
    axes[1].set_xlabel('Observed maximum pair payload budget (KiB)')
    axes[1].set_ylabel('Regimes with upper MAE bound / J <= 25% (out of 84)')
    axes[1].legend(fontsize=7); axes[1].grid(alpha=.2)
    fig.suptitle('Exploratory pilot diagnostics: frontier membership and useful resolution')
    fig.tight_layout(); save(fig, 'overview')

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    cmap = ListedColormap([colors[n] for n in ALGORITHMS])
    norm = BoundaryNorm(np.arange(-.5, len(ALGORITHMS)+.5), len(ALGORITHMS))
    code = {name:i for i, name in enumerate(ALGORITHMS)}
    winners = (selected.sort_values(['mae', 'payload_bytes', 'algorithm'])
               .drop_duplicates(['regime_id', 'budget_bytes']))
    for row, k in enumerate((20, 24)):
        for col, ratio in enumerate((1, 10, 100)):
            local = winners[(winners.k == k) & (winners.n_b == 1000000) & (winners.size_ratio == ratio)]
            js = sorted(local.target_jaccard.unique())
            matrix = np.full((len(BUDGETS), len(js)), np.nan)
            for r in local.itertuples():
                matrix[BUDGETS.index(r.budget_bytes), js.index(r.target_jaccard)] = code[r.algorithm]
            ax = axes[row, col]
            ax.imshow(matrix, origin='lower', aspect='auto', cmap=cmap, norm=norm)
            ax.set_xticks(np.arange(len(js)), [f'{j:g}' for j in js], rotation=35)
            ax.set_yticks(np.arange(len(BUDGETS)), [str(b//1024) for b in BUDGETS])
            ax.set_title(f'k={k}, size ratio={ratio}, larger set=1e6')
            ax.set_xlabel('Target Jaccard'); ax.set_ylabel('Pair budget (KiB)')
    fig.legend(handles=[Patch(facecolor=colors[n], label=LABELS[n]) for n in ALGORITHMS],
               loc='outside lower center', ncol=4, fontsize=9)
    fig.suptitle('Lowest observed MAE; parameter chosen only by payload capacity\n'
                 'Exploratory grid, 20 replicates; HLL/SetSketch have even p only; BinDash caps at 64 KiB', fontsize=12)
    save(fig, 'accuracy_budget_map')

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    cases = [(1, .001), (1, .5), (10, .05), (100, .01)]
    for ax, (ratio, j) in zip(axes.flat, cases):
        local = config[(config.k == 20) & (config.n_b == 1000000) &
                       (config.size_ratio == ratio) & (config.target_jaccard == j)]
        for name in ALGORITHMS:
            g = local[local.algorithm == name]
            ax.scatter(g.payload_bytes, g.mae, color=colors[name], label=LABELS[name], s=26)
            front = g[g.error_space_frontier]
            ax.scatter(front.payload_bytes, front.mae, s=70, edgecolors=colors[name], facecolors='none')
        zero = local.ground_truth.iloc[0]
        ax.axhline(zero, color='gray', ls=':', lw=1, label='Always-zero baseline')
        ax.set_xscale('log'); ax.set_yscale('symlog', linthresh=1e-7)
        ax.set_ylim(bottom=0); ax.grid(alpha=.2)
        ax.set_xlabel('Mean measured pair payload (bytes)'); ax.set_ylabel('Jaccard MAE')
        ax.set_title(f'k=20, larger set=1e6, ratio={ratio}, J={j:g}')
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=4, fontsize=8)
    fig.suptitle('Accuracy-space tradeoffs; rings mark the registered two-objective frontier')
    save(fig, 'representative_frontiers')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    ps = sorted(hll_ss.parameter.unique())
    for ax, field, title in zip(axes, ['jaccard_mae_ratio_ss_to_hll', 'cardinality_mae_ratio_ss_to_hll'],
                              ['Jaccard error at equal payload', 'Cardinality error at equal payload']):
        ax.boxplot([hll_ss.loc[hll_ss.parameter == p, field] for p in ps], tick_labels=ps,
                   showfliers=False, widths=.6)
        ax.axhline(1, color='gray', ls='--'); ax.set_yscale('log')
        ax.set_xlabel('Register precision p'); ax.set_ylabel('SetSketch MAE / HLL MAE')
        ax.set_title(title); ax.grid(axis='y', alpha=.2)
    fig.suptitle('84 regimes per p; below 1 favors SetSketch; exploratory comparisons')
    save(fig, 'hll_setsketch')


def write_report(output, run, status, progress, summary, comparisons, hll_ss, coverage, configurations):
    zone = dt.timezone(dt.timedelta(hours=8))
    start = dt.datetime.fromisoformat(status['started_utc'])
    finished = dt.datetime.fromisoformat(status['finished_utc'])
    measured = dt.datetime.fromisoformat(progress['updated_utc'])
    total = (finished-start).total_seconds()/60
    lines = ['# RabbitSketch2：完整 synthetic-set pilot 深入分析', '',
             f"实验与自动分析于北京时间 **{finished.astimezone(zone):%Y-%m-%d %H:%M:%S}** 完成，"
             f'总耗时 {total:.1f} 分钟，其中测量约 {(measured-start).total_seconds()/60:.1f} 分钟，'
             f'自动分析及出图约 {(finished-measured).total_seconds():.0f} 秒。', '',
             '90,720/90,720 次测量完成，108 个批次均正常退出；7 个候选、54 个原生配置、'
             '84 个可行数据场景，每个配置/场景 20 次配对重复。MinHash 即当前 FastKMV；'
             'OrderMinHash 按用户要求不在实验和后续推荐范围内。', '',
             '**预注册结论：目前没有算法达到淘汰门槛。** 以下统一预算、相对误差和参数覆盖分析'
             '均为事后诊断，不修改原始筛选阈值，也不把事后排名作为已验证的推荐边界。', '',
             '## 1. 前沿与淘汰证据', '',
             '| 算法 | 四维前沿场景 /84 | 误差–字节前沿场景 /84 | 稳健支配场景 /84 | 空采样次数 |',
             '|---|---:|---:|---:|---:|']
    for r in summary.itertuples():
        lines.append(f'| {LABELS[r.algorithm]} | {r.four_objective_frontier_regimes} | '
                     f'{r.error_space_frontier_regimes} | {r.robustly_dominated_regimes} | {r.empty_sample_pairs} |')
    sampled = configurations[(configurations.algorithm == 'kssd') & (configurations.parameter > 0)]
    sampled_front = int(sampled.groupby('regime_id').error_space_frontier.any().sum())
    lines += ['', '四维前沿使用 MAE、成对 payload、构建时间和查询时间；'
              '“进入某场景前沿”表示该算法至少一个配置进入。稳健支配仍使用原协议的配对 bootstrap、'
              '实际资源增益与全部配置条件。不能用配置级支配比例替代算法级淘汰条件。', '',
              f'KSSD 的 `drlevel=0` 在全部 84 个场景中得到零误差，其最大成对 payload 为 '
              f'{configurations.loc[(configurations.algorithm == "kssd") & (configurations.parameter == 0), "payload_bytes"].max()/1e6:.0f} MB。'
              f'这一全量保留配置解释了 KSSD 的 84/84 前沿覆盖；仅看 `drlevel>0`，'
              f'它仍在 {sampled_front}/84 个场景进入误差–字节前沿。', '',
              'MinHash 与 unweighted ProbMinHash 的前沿覆盖较少，值得继续验证冗余假设；'
              '但当前证据不满足预注册的 ≥90% 场景要求，不能宣布淘汰。', '',
              '![前沿和相对精度覆盖](overview.png)', '',
              '## 2. 按真实 payload 预算比较', '',
              '预算为 A+B 两个 sketch 的总 payload，不是单个 sketch，也不含对象/缓存/构建临时内存。'
              '每个算法/场景选取在 20 次重复中**最大观测 payload 均不超过预算**、且平均 payload 最大的配置；'
              '选参不读取估计误差。这个观测上限不保证未来样本也满足同一内存上限。', '',
              '2、8、32 KiB 是本轮固定大小算法共同覆盖较完整的预算点：BinDash、MinHash、'
              'ProbMinHash、HLL、SetSketch 均有对应原生分辨率；短集合可能提前饱和。', '',
              '| 成对预算 | BinDash MAE 更低的场景 /84 | BinDash/MinHash MAE 比值中位数 | 配对 95% 区间完全低于 0 的场景 |',
              '|---|---:|---:|---:|']
    for budget in SHARED_BUDGETS:
        g = comparisons[(comparisons.budget_bytes == budget) & (comparisons.challenger == 'bindash')]
        lines.append(f'| {budget//1024} KiB | {int((g.mae_ratio < 1).sum())}/{len(g)} | '
                     f'{g.mae_ratio.median():.3f} | {int((g.paired_mae_difference_ci_upper < 0).sum())}/{len(g)} |')
    lines += ['', '这些配对区间是在单个已观测配置对上的探索性结果，未做多重比较校正；'
              '场景间也不等同于独立生物学样本。现阶段可以把 BinDash 作为紧凑 Jaccard 比较的重点验证候选，'
              '还不能据此宣布所有真实基因组任务的赢家。', '',
              '**参数覆盖会影响排名。** HLL/SetSketch 只扫了偶数 p，因此在 1、4、16、64 KiB 等预算下'
              '经常只用掉一半预算；BinDash 本轮最多 16,384 bins，即成对 64 KiB。'
              '128/256 KiB 下 BinDash 排名变化可能来自测试上限，不能解释为真实算法切换边界。'
              'KSSD 的 reduction 网格也很粗。完整的预算利用率见 `budget_coverage.csv`。', '',
              '![预算与数据场景](accuracy_budget_map.png)', '',
              '## 3. 低 Jaccard 需要额外检查分辨能力', '',
              '当真实 J=0.001 时，始终预测 0 的 MAE 就是 0.001，已满足绝对误差阈值 0.005 或 0.01。'
              '因此原始 query-reuse 表中的“最快且满足误差阈值”并不一定能分辨低相似度样本。', '',
              '补充诊断使用 `MAE 的 95% 上界 / J <= 25%`；'
              '`relative_guarded_reuse.csv` 进一步要求该上界同时 <=0.005。'
              '这是事后敏感性分析，保留原始绝对误差分析及淘汰规则。没有配置满足时明确标记不可行。', '',
              '| 成对预算 | MinHash 达标场景 /84 | BinDash 达标场景 /84 | SetSketch 达标场景 /84 |',
              '|---|---:|---:|---:|']
    for budget in SHARED_BUDGETS:
        g = coverage[coverage.budget_bytes == budget].set_index('algorithm')
        counts = [int(g.loc[name, 'relative_accuracy_pass_regimes']) for name in ('minhash', 'bindash', 'setsketch')]
        lines.append(f'| {budget//1024} KiB | {counts[0]} | {counts[1]} | {counts[2]} |')
    lines += ['', '![代表性误差–空间曲线](representative_frontiers.png)', '',
              '## 4. HLL 与 SetSketch 的相同精度参数对照', '']
    mid = hll_ss[hll_ss.parameter <= 14]
    hi = hll_ss[hll_ss.parameter == 16]
    lines += [f'相同 p、相同 payload 下，p=4–14 的 {len(mid)} 项配置/场景对照中，'
              f'SetSketch 有 {int((mid.jaccard_mae_ratio_ss_to_hll < 1).sum())} 项 Jaccard MAE 更低；'
              f'p=16 时为 {int((hi.jaccard_mae_ratio_ss_to_hll < 1).sum())}/{len(hi)}。'
              '趋势随精度变化，因此不能概括成一方全面支配另一方。', '',
              'p=16 的反转集中在较小集合：较大集合基数为 10^5、10^6、10^7 时，'
              'SetSketch 的 Jaccard MAE 更低的场景分别为 6/28、7/28、27/28。'
              '这需要按当前实现的估计器适用范围核查，不能直接解释为 SetSketch 范式的固有劣势。', '',
              f'基数误差采用每次重复中 A、B 两侧相对绝对误差的平均值。'
              f'全部同 p 对照中，SetSketch/HLL 基数 MAE 比值的中位数为 '
              f'{hll_ss.cardinality_mae_ratio_ss_to_hll.median():.3f}。'
              '独立的基数角色仍需要保留；本轮 Jaccard 筛选不能移除该角色。', '',
              '![HLL和SetSketch](hll_setsketch.png)', '']
    diagnostic = output/'setsketch_diagnostic/summary.json'
    if diagnostic.exists():
        data = json.loads(diagnostic.read_text())
        observed = next(d for d in data if d['precision']==16)
        diagnostic_lines = ['### SetSketch 小集合异常的独立复现', '',
            '使用原来的 main 入口、64 线程、相同输入种子及原生库，给冻结 worker 的副本增加诊断输出。'
            'p=12/14/16 各复现 20 个原始 pair，全部默认估计值、输入摘要和 payload 与原记录完全一致。'
            '诊断调用改变了查询顺序，因此这些新耗时不参与性能统计。', '',
            '固定 k=20、|A|=10,000、|B|=100,000，真实 J=0.001001001。', '',
            '| p | 默认 J 估计均值 | 默认估计器 MAE | 实例 inclusion–exclusion MAE | 双方均空的寄存器比例 |',
            '|---|---:|---:|---:|---:|']
        for d in data:
            diagnostic_lines.append(f"| {d['precision']} | {d['default_estimate_mean']:.8g} | "
                                    f"{d['default_mae']:.8g} | {d['inclexcl_mae']:.8g} | "
                                    f"{100*d['both_zero_fraction_mean']:.2f}% |")
        diagnostic_lines += ['',
            f"p=16 时默认估计器的 MAE 为 {observed['default_mae']:.6f}，"
            f"现有实例 inclusion–exclusion 路径为 {observed['inclexcl_mae']:.6f}；"
            '说明该异常与默认查询路径有关，不能简单归咎于寄存器容量不足。'
            'p=12/14 下 inclusion–exclusion 又不占优，因而不能未经验证就全局替换估计器。', '',
            '源码核查：当前基数估计已使用 partitioned register 的边际似然；'
            '默认 joint 内核仍基于寄存器相等/大于/小于计数，并在双方均空或饱和的寄存器超过 30% 时回退。'
            f"本例双方均空比例约 {100*observed['both_zero_fraction_mean']:.1f}%，20 次均未触发回退。"
            '这支持“低负载区间的 joint 模型或回退策略需要修正”的假设；尚未完成修复验证。', '',
            '不能仅把 30% 阈值调低：单独诊断现有静态 `jaccardFromCores` 路径时，'
            f"p=16 的 MAE 为 {observed['static_fallback_mae']:.3f}，输出还可能超出 [0,1]。"
            '该静态路径用 `factor/sum` 计算 union，与实例路径的基数处理不同。'
            '它在这组原始测量中并未实际触发；原始 Jaccard 估计仍全部位于 [0,1]。', '',
            '定位：[joint 内核](/home/bigssd/zt/rabbitsketch2_benchmark/source/RabbitSketch-main/src/SetSketch.cpp:490)，'
            '[实例 inclusion–exclusion](/home/bigssd/zt/rabbitsketch2_benchmark/source/RabbitSketch-main/src/SetSketch.cpp:369)，'
            '[静态 fallback](/home/bigssd/zt/rabbitsketch2_benchmark/source/RabbitSketch-main/src/SetSketch.cpp:976)。', '',
            '复现入口：`benchmarks/screening/diagnose_setsketch.py`；完整输出与编译/链接摘要在 '
            '`setsketch_diagnostic/`。原生库、原始记录和预注册结论均保留。', '']
        (output/'setsketch_anomaly.md').write_text('\n'.join(diagnostic_lines))
        lines += diagnostic_lines
    lines += ['## 5. 空采样与大小不平衡', '',
              '420 次测量有一侧或两侧采样为空，其中 FracMinHash 346 次、KSSD 74 次。'
              '对应 62 个配置/场景已按原协议排除出前沿；空采样不是进程失败。'
              '它们集中在小集合与较强 reduction/scaled 配置。', '',
              '对 FracMinHash，小集合的预期样本数约为 `|A|/scaled`。例如 |A|=1000、scaled=10000 时'
              '仅为 0.1；预算选择需要检查小集合保留的样本数，不能只看大集合大小。'
              '这一轮测量的是 Jaccard，尚未测 directed containment，不能从本轮给出 containment 的算法切换阈值。', '',
              '## 6. 性能结果的适用范围', '',
              '构建和查询耗时来自一个 main 进程的 64 线程样本列表负载；每个线程测量一个 pair。'
              '总批次时间包含合成集合生成、N 分隔序列构造、哈希、sketch 构建和重复查询。'
              '批次吞吐量不等于纯 genome 构建吞吐量，也不等于 all-vs-all 查询吞吐量。', '',
              'RSS 是整个共享进程的高水位，包含所有并行输入和分配器保留。'
              'KSSD 进程的最大值约 '
              f'{summary.loc[summary.algorithm == "kssd", "maximum_shared_process_rss_bytes"].iloc[0]/1024**3:.1f} GiB；'
              '它不能充当某一个 sketch 的内存占用。字节预算只使用实测 finalized query-array payload。', '',
              '## 7. 下一轮应验证的内容', '',
              '1. 优先核查 SetSketch 低负载 joint 估计器与静态 fallback；用小、中、大集合及不同 p 验证任何修复，'
              '为修复版本建立独立结果，避免与当前快照混合。',
              '2. 在预先指定的 2/8/32 KiB 成对预算下，对真实基因组的精确 canonical k-mer Jaccard 做独立验证，'
              '重点复核 BinDash、当前 MinHash 与 unweighted ProbMinHash 的差异。七个候选仍保留在注册初筛结果中。',
              '3. 补齐预算网格：HLL/SetSketch 加入所需奇数 p；如考察 128/256 KiB，扩展 BinDash bins，'
              '并核对其他算法的预算覆盖，再生成推荐边界。',
              '4. 单独开展 directed containment 与 abundance 实验；HLL/SetSketch 的 cardinality 角色独立评估。'
              'OrderMinHash 和专门的顺序实验已排除。',
              '5. 在完成生物学验证后，再对保留方法测量真实 FASTA 构建、重复查询和 all-vs-all 系统性能。', '',
              '## 可复现文件', '',
              f'- 原始完成实验（只读）：`{run}`。',
              '- `followup_manifest.json`：输入 SHA-256、原始行摘要核验、脚本版本和补充分析规则。',
              '- `algorithm_summary.csv`：预注册前沿、空采样与淘汰状态。',
              '- `budget_selected_configurations.csv`、`budget_coverage.csv`：按 payload 选参及覆盖。',
              '- `paired_budget_comparisons.csv`：与 MinHash 的配对误差对照及区间。',
              '- `hll_setsketch_comparisons.csv`：相同 p/相同 payload 的 Jaccard 和 cardinality 对照。',
              '- `relative_guarded_reuse.csv`：补充相对精度条件下的查询复用选择。',
              '- `empty_sampling.csv`：空采样配置/场景明细。',
              '- 所有图同时保存 PNG 和独立 PDF，可供审阅和导出。', '']
    (output/'report.md').write_text('\n'.join(lines))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, default=DEFAULT_ROOT/'runs/pilot')
    p.add_argument('--output', type=Path, default=HERE/'results/pilot_20261001')
    args = p.parse_args()
    run, output = args.run.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    plan = json.loads((run/'plan.json').read_text())
    status = json.loads((run/'status.json').read_text())
    progress = json.loads((run/'progress.json').read_text())
    if status['phase'] != 'completed' or status['failed'] or progress['completed'] != len(plan['jobs']):
        raise ValueError('this follow-up requires a complete, successful pilot')
    source_manifest = json.loads((run/'analysis/analysis_manifest.json').read_text())
    rows = [json.loads(line) for line in (run/'raw.jsonl').open()]
    rows.sort(key=lambda r:r['job_id'])
    if digest(rows) != source_manifest['raw_rows_sha256'] or digest(plan) != source_manifest['plan_sha256']:
        raise ValueError('raw measurements or plan differ from the registered analysis')
    assert len(rows) == len(plan['jobs']) == len({r['job_id'] for r in rows})
    assert all(r['status'] == 'ok' and r['batch_threads'] == 64 and
               r['measurement_profile'] == 'main_openmp_batch' for r in rows)
    configurations = pd.read_csv(run/'analysis/configurations.csv')
    decisions = pd.read_csv(run/'analysis/algorithm_decisions.csv')
    grouped = {}
    for row in rows:
        grouped.setdefault((row['regime_id'], row['config_id']), []).append(row)
    extra, errors, card_errors, bootstrap = [], {}, {}, {}
    for key, g in grouped.items():
        g.sort(key=lambda r:r['replicate'])
        assert [r['replicate'] for r in g] == list(range(plan['config']['replicates']))
        errors[key] = np.array([abs(r['estimate']-r['ground_truth']) for r in g])
        if g[0]['algorithm'] in ('hll', 'setsketch'):
            card_errors[key] = np.array([sum(abs(r[f'cardinality_estimate_{side}']-r[f'n_{side}'])/
                                            r[f'n_{side}'] for side in ('a','b'))/2 for r in g])
        extra.append(dict(regime_id=key[0], config_id=key[1],
                          max_observed_payload_bytes=max(r['sketch_payload_bytes'] for r in g),
                          empty_sample_a_count=sum(bool(r['empty_sample_a']) for r in g),
                          empty_sample_b_count=sum(bool(r['empty_sample_b']) for r in g),
                          checked_mae=float(errors[key].mean())))
    configurations = configurations.merge(pd.DataFrame(extra), on=['regime_id','config_id'], validate='one_to_one')
    assert np.allclose(configurations.mae, configurations.checked_mae, rtol=1e-10, atol=1e-15)
    for regime in plan['regimes']:
        seed = int(digest([plan['config']['analysis']['bootstrap_seed'], regime['regime_id']])[:16], 16)
        bootstrap[regime['regime_id']] = np.random.default_rng(seed).integers(
            0, plan['config']['replicates'], (plan['config']['analysis']['bootstrap_samples'], plan['config']['replicates']))
    summary = overview(configurations, decisions)
    selected = budget_selections(configurations, BUDGETS)
    comparisons = compare_budgets(selected, errors, bootstrap)
    hll_ss = hll_setsketch(configurations, errors, card_errors, bootstrap)
    coverage = (selected.assign(relative_accuracy_pass=selected.relative_mae_ci_upper <= RELATIVE_LIMIT)
                .groupby(['budget_bytes','budget_kib','algorithm']).agg(
                    available_regimes=('regime_id','size'),
                    median_budget_used_fraction=('budget_used_fraction','median'),
                    full_budget_regimes=('budget_used_fraction', lambda s:int((s >= .99).sum())),
                    median_relative_mae=('relative_mae','median'),
                    relative_accuracy_pass_regimes=('relative_accuracy_pass','sum')).reset_index())
    winners = selected.sort_values(['mae','payload_bytes','algorithm']).drop_duplicates(['regime_id','budget_bytes'])
    guarded = relative_guarded_reuse(configurations, plan['regimes'])
    empty = configurations[configurations.empty_samples > 0]
    for name, table in [('algorithm_summary',summary), ('budget_selected_configurations',selected),
                        ('paired_budget_comparisons',comparisons), ('budget_coverage',coverage),
                        ('hll_setsketch_comparisons',hll_ss), ('observed_accuracy_winners',winners),
                        ('relative_guarded_reuse',guarded), ('empty_sampling',empty)]:
        table.to_csv(output/(name+'.csv'),index=False)
    figures(output, configurations, selected, summary, hll_ss)
    write_report(output, run, status, progress, summary, comparisons, hll_ss, coverage, configurations)
    inputs = ['plan.json','raw.jsonl','status.json','progress.json','analysis/configurations.csv',
              'analysis/algorithm_decisions.csv','analysis/analysis_manifest.json','analysis/batches.csv']
    atomic_json(output/'followup_manifest.json', dict(
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), source_run=str(run),
        source_files_sha256={name:file_sha(run/name) for name in inputs},
        raw_rows_digest_verified=source_manifest['raw_rows_sha256'], script_sha256=file_sha(__file__),
        budgets_bytes=BUDGETS, fixed_size_common_budgets_bytes=SHARED_BUDGETS,
        parameter_selection='largest mean payload with maximum observed pair payload <= budget; no error tuning',
        relative_accuracy_limit=RELATIVE_LIMIT, additional_absolute_limit_for_reuse=.005,
        purpose='post_hoc_diagnostics_only; registered pruning thresholds unchanged',
        diagnostic_manifest_sha256=file_sha(output/'setsketch_diagnostic/manifest.json')
            if (output/'setsketch_diagnostic/manifest.json').exists() else None,
        output_files_sha256={f.name:file_sha(f) for f in output.iterdir()
                             if f.is_file() and f.name!='followup_manifest.json' and f.suffix!='.log'}))
    print(summary.to_string(index=False))
    print('Report:',output/'report.md')
    print('BinDash vs MinHash at shared budgets:')
    print(comparisons[(comparisons.challenger=='bindash') & comparisons.budget_bytes.isin(SHARED_BUDGETS)]
          .groupby('budget_kib').agg(cells=('regime_id','size'),wins=('mae_ratio',lambda s:int((s<1).sum())),
                                   median_ratio=('mae_ratio','median'),
                                   paired_upper_below_zero=('paired_mae_difference_ci_upper',lambda s:int((s<0).sum()))).to_string())


if __name__ == '__main__':
    main()
