#!/usr/bin/env python3
"""Six-candidate evidence tables. Keep source runs and timing populations separate."""
import argparse
import datetime as dt
import itertools
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from run import DEFAULT_ROOT, HERE, atomic_json, digest, file_sha, make_plan
from analyze import aggregate, decisions, dominates
from summarize_pilot import budget_selections, paired_interval, hll_setsketch

NAMES=['minhash','fracminhash','hll','setsketch','probminhash','bindash']
LABELS=dict(zip(NAMES,['MinHash','FracMinHash','HLL','SetSketch corrected','ProbMinHash unweighted','BinDash']))
BUDGETS=[2048,8192,32768]


def read_rows(run):
    return [json.loads(line) for line in (run/'raw.jsonl').open()]


def extras(plan,rows):
    groups={}; errors={};cards={};bootstrap={}; extra=[]
    for r in rows: groups.setdefault((r['regime_id'],r['config_id']),[]).append(r)
    n=plan['config']['replicates']
    for key,g in groups.items():
        g.sort(key=lambda r:r['replicate'])
        if [r['replicate'] for r in g]!=list(range(n)): raise ValueError('missing/duplicate replicate')
        errors[key]=np.array([abs(r['estimate']-r['ground_truth']) for r in g])
        if g[0]['algorithm'] in ('hll','setsketch'):
            cards[key]=np.array([sum(abs(r[f'cardinality_estimate_{s}']-r[f'n_{s}'])/r[f'n_{s}']
                                    for s in ('a','b'))/2 for r in g])
        extra.append(dict(regime_id=key[0],config_id=key[1],
                          max_observed_payload_bytes=max(r['sketch_payload_bytes'] for r in g)))
    for r in plan['regimes']:
        seed=int(digest([plan['config']['analysis']['bootstrap_seed'],r['regime_id']])[:16],16)
        bootstrap[r['regime_id']]=np.random.default_rng(seed).integers(0,n,(2000,n))
    return pd.DataFrame(extra),errors,cards,bootstrap


def paired_comparisons(selected,errors,bootstrap):
    result=[]
    for (budget,regime),g in selected.groupby(['budget_kib','regime_id']):
        g=g.assign(algorithm_order=g.algorithm.map({n:i for i,n in enumerate(NAMES)})).sort_values('algorithm_order')
        for a,b in itertools.combinations(g.itertuples(),2):
            lo,hi=paired_interval(errors[(regime,a.config_id)],errors[(regime,b.config_id)],bootstrap[regime])
            result.append(dict(budget_kib=budget,regime_id=regime,k=a.k,n_b=a.n_b,
                               size_ratio=a.size_ratio,target_jaccard=a.target_jaccard,
                               a=a.algorithm,b=b.algorithm,parameter_a=a.parameter,parameter_b=b.parameter,
                               mae_a=a.mae,mae_b=b.mae,ratio_a_to_b=a.mae/b.mae if b.mae else np.nan,
                               difference_ci_lower=lo,difference_ci_upper=hi))
    return pd.DataFrame(result)


def synthetic(original,correction,out,baseline_only=False):
    plan=make_plan(json.loads((HERE/'configs/pilot_six.json').read_text()))
    rows=[r for r in read_rows(original) if r['algorithm']!='kssd']
    if len(rows)!=84000 or any(r['status']!='ok' for r in rows): raise ValueError('incomplete original six')
    cache=out/'original_six_configurations.csv'
    cache_identity=dict(rows_digest=digest(sorted(rows,key=lambda r:r['job_id'])),plan_digest=digest(plan),
                        analyzer_sha256=file_sha(HERE/'analyze.py'))
    cache_manifest=out/'original_six_manifest.json'
    if cache.exists() and cache_manifest.exists() and json.loads(cache_manifest.read_text())==cache_identity:
        old=pd.read_csv(cache)
    else:
        old=pd.DataFrame(aggregate(plan,rows));old.to_csv(cache,index=False)
        pd.DataFrame(decisions(plan,old.to_dict('records'))).to_csv(out/'original_six_decisions.csv',index=False)
        atomic_json(cache_manifest,cache_identity)
    if baseline_only: return
    new=read_rows(correction)
    newplan=json.loads((correction/'plan.json').read_text())
    progress=json.loads((correction/'progress.json').read_text())
    if progress['phase']!='measured' or progress['failed'] or len(new)!=11760: raise ValueError('correction incomplete')
    old_ss={r['job_id']:r for r in rows if r['algorithm']=='setsketch'}
    for r in new:
        o=old_ss[r['job_id']]
        if r['status']!='ok' or any(r[k]!=o[k] for k in ('input_digest_a','input_digest_b','ground_truth',
                'payload_bytes_a','payload_bytes_b','cardinality_estimate_a','cardinality_estimate_b')):
            raise ValueError('correction changed data/state/cardinality, or failed')
    ss=pd.DataFrame(aggregate(newplan,new))
    config=pd.concat([old[old.algorithm!='setsketch'],ss],ignore_index=True)
    config['source_run']=np.where(config.algorithm=='setsketch',str(correction),str(original))
    extension=DEFAULT_ROOT/'runs/frac_budget_extension_20261001'
    extrows=read_rows(extension)
    extplan=json.loads((extension/'plan.json').read_text())
    extprogress=json.loads((extension/'progress.json').read_text())
    if extprogress['phase']!='measured' or extprogress['failed'] or len(extrows)!=2240:
        raise ValueError('FracMinHash budget extension incomplete')
    data_checks={(r['regime_id'],r['replicate']):(r['input_digest_a'],r['input_digest_b']) for r in rows}
    for r in extrows:
        if (r['input_digest_a'],r['input_digest_b'])!=data_checks[(r['regime_id'],r['replicate'])]:
            raise ValueError('extension input differs')
    extended=pd.DataFrame(aggregate(extplan,extrows));extended['source_run']=str(extension)
    config=pd.concat([config,extended],ignore_index=True)
    config['error_space_frontier']=False
    config['descriptive_cross_run_4d_frontier']=False
    # Recompute from all six, never retain flags computed with KSSD present.
    for regime,g in config[config.eligible].groupby('regime_id'):
        records=g.to_dict('records')
        for idx,r in zip(g.index,records):
            config.loc[idx,'error_space_frontier']=not any(dominates(x,r,('mae','payload_bytes')) for x in records)
            config.loc[idx,'descriptive_cross_run_4d_frontier']=not any(dominates(x,r) for x in records)
    config=config.drop(columns=['pareto_frontier','robustly_dominated','dominator_config_ids'],errors='ignore')
    current=[r for r in rows if r['algorithm']!='setsketch']+new+extrows
    extra,errors,cards,bootstrap=extras(plan,current)
    config=config.merge(extra,on=['regime_id','config_id'],validate='one_to_one')
    selected=budget_selections(config,BUDGETS)
    selected['relative_pass']=selected.relative_mae_ci_upper<=.25
    coverage=selected.groupby(['budget_kib','algorithm']).agg(
        available=('regime_id','size'),relative_pass=('relative_pass','sum'),
        median_relative_mae=('relative_mae','median'),median_absolute_mae=('mae','median'),
        median_used_fraction=('budget_used_fraction','median')).reset_index()
    comparisons=paired_comparisons(selected,errors,bootstrap)
    hll_ss=hll_setsketch(config,errors,cards,bootstrap)
    dimension=[]
    for name in ['target_jaccard','n_b','size_ratio','k']:
        x=selected.groupby(['budget_kib',name,'algorithm']).agg(
            cells=('regime_id','size'),relative_pass=('relative_pass','sum'),
            median_relative_mae=('relative_mae','median'),median_absolute_mae=('mae','median')).reset_index()
        x=x.rename(columns={name:'value'});x['dimension']=name;dimension.append(x)
    winners=selected.sort_values(['mae','payload_bytes','algorithm']).drop_duplicates(['budget_kib','regime_id']).copy()
    precision=selected.groupby(['budget_kib','regime_id']).relative_pass.any().rename('any_precision_pass').reset_index()
    winners=winners.merge(precision,on=['budget_kib','regime_id'],validate='one_to_one')
    winners['interpretation']='observed minimum MAE, exploratory; no automatic statistical winner'
    # Preserve all three dimensions when assessing imbalance; feasible J ranges differ.
    balanced=selected[selected.size_ratio==1]
    imbalanced=selected[selected.size_ratio!=1].merge(balanced,on=['algorithm','budget_kib','k','n_b','target_jaccard'],suffixes=('','_balanced'))
    imbalanced['relative_mae_ratio_to_balanced']=imbalanced.relative_mae/imbalanced.relative_mae_balanced
    changes=old[old.algorithm=='setsketch'].merge(ss,on=['regime_id','config_id'],suffixes=('_old','_corrected'))
    changes['mae_ratio_corrected_to_old']=changes.mae_corrected/changes.mae_old
    fronts=config.groupby(['algorithm','regime_id'])[['error_space_frontier','descriptive_cross_run_4d_frontier']].any().groupby('algorithm').sum().reset_index()
    for name,table in [('synthetic_configurations',config),('synthetic_budget_selected',selected),
                       ('synthetic_coverage',coverage),('synthetic_paired_comparisons',comparisons),
                       ('synthetic_dimensions',pd.concat(dimension)),('synthetic_winners',winners),
                       ('synthetic_hll_setsketch',hll_ss),('synthetic_imbalance_matched',imbalanced),
                       ('setsketch_before_after',changes),('synthetic_frontiers',fronts)]:
        table.to_csv(out/(name+'.csv'),index=False)
    atomic_json(out/'synthetic_manifest.json',dict(original=str(original),correction=str(correction),
        original_raw_sha256=file_sha(original/'raw.jsonl'),corrected_raw_sha256=file_sha(correction/'raw.jsonl'),
        frac_extension=str(extension),frac_extension_raw_sha256=file_sha(extension/'raw.jsonl'),
        selected_config='largest mean payload, maximum observed payload <= budget; no tuning on error',
        budgets_bytes=BUDGETS,rows=len(current),regimes=84,replicates=20,
        ci='paired pointwise percentile bootstrap; 2000 resamples; not familywise adjusted',
        timing='original run + SetSketch correction + FracMinHash grid extension: 4D frontier descriptive only; no revised pruning claim',
        excluded_by_user=['kssd','orderminhash'],minhash='FastKMV; no separate Legacy entry'))
    print(coverage.to_string(index=False),flush=True)


def real(run,out):
    progress=json.loads((run/'progress.json').read_text())
    if progress['phase']!='completed' or progress['failed']: raise ValueError('real validation incomplete')
    rows=read_rows(run)
    if len(rows)!=32400 or any(r['status']!='ok' for r in rows): raise ValueError('incomplete real rows')
    df=pd.DataFrame(rows)
    if set(df.measurement_profile)!= {'main_openmp_real_pair'}: raise ValueError('profile mismatch')
    df['absolute_error']=(df.estimate-df.ground_truth).abs()
    df['signed_error']=df.estimate-df.ground_truth
    df['relative_error']=df.absolute_error/df.ground_truth.where(df.ground_truth>0)
    group=['pair_id','tier','k','algorithm','parameter']
    config=df.groupby(group).agg(mae=('absolute_error','mean'),rmse=('signed_error',lambda x:np.sqrt((x*x).mean())),
        bias=('signed_error','mean'),relative_mae=('relative_error','mean'),ground_truth=('ground_truth','first'),
        payload_bytes=('sketch_payload_bytes','mean'),max_observed_payload_bytes=('sketch_payload_bytes','max'),
        build_seconds=('build_seconds','median'),query_seconds=('query_seconds','median'),
        cold_query_seconds=('query_cold_seconds','median'),n_a=('n_a','first'),n_b=('n_b','first'),
        seeds=('hash_seed','nunique'),empty_a=('empty_sample_a','sum'),empty_b=('empty_sample_b','sum')).reset_index()
    if set(config.seeds)!={3}: raise ValueError('not all three seeds')
    config['eligible']=(config.empty_a==0)&(config.empty_b==0)
    selections=[]
    for budget in BUDGETS:
        s=config[config.eligible & (config.max_observed_payload_bytes<=budget)].sort_values(
            ['payload_bytes','parameter'],ascending=[False,True]).drop_duplicates(['pair_id','k','algorithm']).copy()
        s['budget_kib']=budget/1024;selections.append(s)
    selected=pd.concat(selections,ignore_index=True)
    # Hold biological pairs as the grouping unit. Hash seeds are repeated measurements.
    summary=selected.groupby(['budget_kib','k','algorithm']).agg(pairs=('pair_id','size'),
        mean_absolute_error=('mae','mean'),median_absolute_error=('mae','median'),
        median_relative_mae=('relative_mae','median'),mean_bias=('bias','mean'),
        median_build_seconds=('build_seconds','median'),median_query_seconds=('query_seconds','median')).reset_index()
    comp=[]
    for (budget,k),g in selected.groupby(['budget_kib','k']):
        for a,b in itertools.combinations(NAMES,2):
            paired=g[g.algorithm==a].merge(g[g.algorithm==b],on='pair_id',suffixes=('_a','_b'))
            comp.append(dict(budget_kib=budget,k=k,a=a,b=b,pairs=len(paired),
                a_lower_mae=int((paired.mae_a<paired.mae_b).sum()),ties=int((paired.mae_a==paired.mae_b).sum()),
                pooled_mae_a=paired.mae_a.mean(),pooled_mae_b=paired.mae_b.mean(),
                pooled_mae_ratio=paired.mae_a.mean()/paired.mae_b.mean(),
                median_build_ratio=(paired.build_seconds_a/paired.build_seconds_b).median(),
                median_query_ratio=(paired.query_seconds_a/paired.query_seconds_b).median()))
    # The same assemblies occur in multiple pairs: no IID genome-pair CI is asserted.
    tiers=selected.groupby(['budget_kib','k','tier','algorithm']).agg(pairs=('pair_id','size'),
        mean_absolute_error=('mae','mean'),median_relative_mae=('relative_mae','median'),
        mean_jaccard=('ground_truth','mean')).reset_index()
    by_j=selected.copy();by_j['j_range']=pd.cut(by_j.ground_truth,[-1,0,.001,.01,.05,.1,.5,1],
        labels=['J=0','(0,.001]','(.001,.01]','(.01,.05]','(.05,.1]','(.1,.5]','(.5,1]'])
    j_summary=by_j.groupby(['budget_kib','k','j_range','algorithm'],observed=True).agg(
        pairs=('pair_id','size'),mean_absolute_error=('mae','mean'),median_relative_mae=('relative_mae','median')).reset_index()
    reuse=[];cross=[]
    for (budget,k,pair),g in selected.groupby(['budget_kib','k','pair_id']):
        g=g.assign(algorithm_order=g.algorithm.map({n:i for i,n in enumerate(NAMES)})).sort_values('algorithm_order')
        truth=g.ground_truth.iloc[0]
        pool=config[(config.pair_id==pair)&(config.k==k)&config.eligible&
                    (config.max_observed_payload_bytes<=budget*1024)]
        allowed=pool[pool.mae<=min(.005,.25*truth)] if truth>0 else pool.iloc[0:0]
        for q in [1,100,10000,1000000]:
            if allowed.empty:
                reuse.append(dict(budget_kib=budget,k=k,pair_id=pair,q=q,algorithm='none'))
            else:
                ranked=allowed.assign(cost=allowed.build_seconds+allowed.cold_query_seconds+(q-1)*allowed.query_seconds).sort_values('cost')
                r=ranked.iloc[0]
                reuse.append(dict(budget_kib=budget,k=k,pair_id=pair,q=q,algorithm=r.algorithm,
                                  parameter=r.parameter,seconds=r.cost,mae=r.mae,
                                  payload_bytes=r.payload_bytes))
        for a,b in itertools.combinations(g.itertuples(),2):
            denominator=a.query_seconds-b.query_seconds
            if denominator:
                # T(Q)=B+cold+(Q-1)*warm, Q>=1.
                q=1+(b.build_seconds+b.cold_query_seconds-a.build_seconds-a.cold_query_seconds)/denominator
                if q>=1 and np.isfinite(q):
                    cross.append(dict(budget_kib=budget,k=k,pair_id=pair,a=a.algorithm,b=b.algorithm,q_star=q,
                                      both_observed_precision_pass=truth>0 and max(a.mae,b.mae)<=min(.005,.25*truth),
                                      faster_at_q1=a.algorithm if a.build_seconds+a.cold_query_seconds < b.build_seconds+b.cold_query_seconds else b.algorithm,
                                      faster_at_large_q=a.algorithm if a.query_seconds<b.query_seconds else b.algorithm))
    for name,table in [('real_configurations',config),('real_budget_selected',selected),('real_summary',summary),
                       ('real_pairwise_comparisons',pd.DataFrame(comp)),('real_tiers',tiers),
                       ('real_jaccard_ranges',j_summary),('real_reuse_observed',pd.DataFrame(reuse)),
                       ('real_crossovers',pd.DataFrame(cross))]:
        table.to_csv(out/(name+'.csv'),index=False)
    atomic_json(out/'real_manifest.json',dict(source=str(run),raw_sha256=file_sha(run/'raw.jsonl'),
        rows=len(rows),pairs=200,unique_genomes=progress['unique_genomes'],k=[20,24],seeds=[42,101,202],
        confidence='descriptive real validation; assemblies reused between pairs, no IID pair significance claims',
        scope='native general-set Jaccard, exact canonical k-mer truth; no ANI/abundance/containment claims',
        reuse='all tested native configurations within budget, post-hoc observed MAE gate; warm same-pair cost model, not all-vs-all',
        crossover='fixed largest-payload configurations; both_observed_precision_pass marks acceptable accuracy; not optimal parameter switching'))
    print(summary.to_string(index=False),flush=True)


def figures(out):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/rabbitsketch_screening_mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap,BoundaryNorm
    from matplotlib.patches import Patch
    colors=dict(minhash='#777777',fracminhash='#E69F00',hll='#009E73',
                setsketch='#D55E00',probminhash='#CC79A7',bindash='#0072B2')
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=180,bbox_inches='tight')
        fig.savefig(out/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    s=pd.read_csv(out/'synthetic_budget_selected.csv')
    fig,ax=plt.subplots(figsize=(9,4.8))
    for name in NAMES:
        g=s[s.algorithm==name].groupby('budget_kib').relative_pass.sum()
        ax.plot(g.index,g.values,'-o',label=LABELS[name],color=colors[name])
    ax.set(xscale='log',xlabel='Combined A+B payload budget (KiB)',ylabel='Regimes passing upper MAE / J <= 25% (out of 84)',ylim=(0,86))
    ax.set_xticks([2,8,32],['2','8','32']);ax.grid(alpha=.2);ax.legend(fontsize=8)
    ax.set_title('Corrected six-candidate synthetic screening');save(fig,'synthetic_coverage')
    w=pd.read_csv(out/'synthetic_winners.csv')
    fig,axes=plt.subplots(2,3,figsize=(14,6),constrained_layout=True)
    cmap=ListedColormap([colors[n] for n in NAMES]);norm=BoundaryNorm(np.arange(-.5,6),6)
    for row,k in enumerate([20,24]):
        for col,ratio in enumerate([1,10,100]):
            # Display N_B=1e6, all budgets. Other lengths are in the companion CSV.
            g=w[(w.k==k)&(w.size_ratio==ratio)&(w.n_b==1000000)].copy()
            g['code']=g.algorithm.map({n:i for i,n in enumerate(NAMES)})
            m=g.pivot(index='budget_kib',columns='target_jaccard',values='code').sort_index(ascending=False)
            ax=axes[row,col];ax.imshow(m,cmap=cmap,norm=norm,aspect='auto')
            for r in g.itertuples():
                if not r.any_precision_pass:
                    ax.text(list(m.columns).index(r.target_jaccard),list(m.index).index(r.budget_kib),
                            'x',ha='center',va='center',fontsize=15,color='black')
            ax.set_xticks(range(len(m.columns)),[f'{j:g}' for j in m.columns]);ax.set_yticks(range(len(m.index)),[f'{b:g}' for b in m.index])
            ax.set(title=f'k={k}, size ratio {ratio}:1',xlabel='True Jaccard',ylabel='A+B KiB')
    fig.legend(handles=[Patch(color=colors[n],label=LABELS[n]) for n in NAMES],loc='outside lower center',ncol=3,fontsize=8)
    fig.suptitle('Observed minimum MAE at N_large=10^6; x = no method passes precision gate\nExploratory colors, not significance labels');save(fig,'selection_map')
    if (out/'real_summary.csv').exists():
        r=pd.read_csv(out/'real_summary.csv')
        fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
        for ax,k in zip(axes,[20,24]):
            for name in NAMES:
                g=r[(r.k==k)&(r.algorithm==name)]
                if len(g)!=3 or (g.pairs!=200).any(): continue
                ax.plot(g.budget_kib,g.mean_absolute_error,'-o',label=LABELS[name],color=colors[name])
            ax.set(xscale='log',yscale='log',xlabel='Combined A+B payload budget (KiB)',ylabel='Mean absolute Jaccard error',title=f'200 real genome pairs, k={k}, 3 seeds')
            ax.set_xticks([2,8,32],['2','8','32']);ax.grid(alpha=.2)
        axes[1].legend(fontsize=8);save(fig,'real_accuracy')
        fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
        r8=r[(r.k==20)&(r.budget_kib==8)]
        for row in r8.itertuples():
            for ax,x in zip(axes,[1000*row.median_build_seconds,1e6*row.median_query_seconds]):
                ax.scatter(x,row.mean_absolute_error,s=65,color=colors[row.algorithm],label=LABELS[row.algorithm])
        axes[0].set(xlabel='Median build time for A+B (ms)',ylabel='Mean absolute Jaccard error')
        axes[1].set(xscale='log',xlabel='Median warm query time (microseconds)',ylabel='Mean absolute Jaccard error')
        for ax in axes: ax.grid(alpha=.2)
        axes[1].legend(fontsize=7)
        fig.suptitle('Real genome tradeoffs: k=20, A+B budget 8 KiB, 64-thread batch load')
        save(fig,'real_tradeoffs')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,default=DEFAULT_ROOT/'runs/pilot')
    p.add_argument('--correction',type=Path,default=DEFAULT_ROOT/'runs/setsketch_validated_20261001')
    p.add_argument('--real',type=Path,default=DEFAULT_ROOT/'runs/real_validation_20261001')
    p.add_argument('--output',type=Path,default=HERE/'results/systematic_20261001')
    p.add_argument('--baseline-only',action='store_true')
    p.add_argument('--synthetic-only',action='store_true')
    p.add_argument('--real-only',action='store_true')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    if not args.real_only: synthetic(args.original,args.correction,args.output,args.baseline_only)
    if args.baseline_only: return
    if not args.synthetic_only: real(args.real,args.output)
    figures(args.output)
    atomic_json(args.output/'analysis_version.json',dict(script_sha256=file_sha(__file__),
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),scope=NAMES,
        artifacts={f.name:file_sha(f) for f in args.output.iterdir() if f.suffix in ('.csv','.pdf','.png')}))

if __name__=='__main__': main()
