#!/usr/bin/env python3
"""Evidence-only ANI advisor prototype; never invent unmeasured configurations."""
import argparse,json,os
from pathlib import Path
import pandas as pd
_bundled=Path(__file__).resolve().parents[1]/'results/regimes_20261001'
ROOT=Path(os.environ.get('RABBITSKETCH_ADVISOR_DATA',str(_bundled if _bundled.exists() else '/home/bigssd/zt/rabbitsketch2_benchmark/analysis/regimes_20261001')))
def choose(family,ani_bin,budget,target,incident_reuse,cost_n=10000):
    pair_file=ROOT/'ani_pair_config.csv'
    if not pair_file.exists():pair_file=ROOT/'ani_pair_config.csv.gz'
    scores=pd.read_csv(ROOT/'ani_config_accuracy.csv');pairs=pd.read_csv(pair_file,low_memory=False);costs=pd.read_csv(ROOT/'workload_measured.csv')
    candidates=[]
    for c in costs[(costs.k==24)&(costs.objects==cost_n)].to_dict('records'):
        mask=(scores.family==family)&(scores.ani_bin==ani_bin)&(scores.k==24)&(scores.algorithm==c['algorithm'])&(scores.parameter==c['parameter'])
        cal=scores[mask&(scores.split=='calibration')];val=scores[mask&(scores.split=='validation')]
        if cal.empty or val.empty:continue
        cal=cal.iloc[0];val=val.iloc[0]
        p=pairs[(pairs.family==family)&(pairs.ani_bin==ani_bin)&(pairs.k==24)&(pairs.algorithm==c['algorithm'])&(pairs.parameter==c['parameter'])]
        calbytes=p[p.split=='calibration'].max_pair_bytes.max();valbytes=p[p.split=='validation'].max_pair_bytes.max()
        if cal.mae>target or cal.p95>2*target or calbytes>budget:continue
        candidates.append(dict(algorithm=c['algorithm'],parameter=int(c['parameter']),k=24,calibration_mae=cal.mae,calibration_p95=cal.p95,validation_mae=val.mae,validation_p95=val.p95,max_observed_pair_bytes=float(max(calbytes,valbytes)),validation_pass=bool(val.mae<=target and val.p95<=2*target and valbytes<=budget),modeled_seconds_per_object=c['B']+incident_reuse*.5*c['t_pair']))
    result=dict(status='no_calibration_candidate',reference_family=family,ani_bin=ani_bin,target_mae_fraction=target,target_p95_fraction=2*target,pair_budget_bytes=budget,incident_queries_per_object=incident_reuse,cost_calibration_N=cost_n,cost_scope='64-thread mixed-bacterial all-vs-all throughput, excluding input parsing and matrix serialization; not isolated-pair latency')
    if not candidates:return result
    winner=min(candidates,key=lambda c:c['modeled_seconds_per_object']);result.update(status='validated_in_this_suite' if winner['validation_pass'] else 'calibration_choice_failed_holdout',choice=winner)
    if winner['algorithm']=='setsketch':result['required_variant']='isolated partitioned-joint SetSketch patch'
    return result
def main():
    global ROOT
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--data-dir',type=Path,help='Directory containing the frozen advisor CSV tables');ap.add_argument('--family',choices=['substitution','mash142'],default='substitution');ap.add_argument('--ani-bin',choices=['80–90','90–95','95–99','99–100'],required=True);ap.add_argument('--pair-budget',type=int,default=8192);ap.add_argument('--target-error',type=float,default=.001);ap.add_argument('--queries-per-sketch',type=float,default=1);args=ap.parse_args()
    if args.data_dir is not None:ROOT=args.data_dir
    if args.pair_budget<=0 or args.target_error<=0 or args.queries_per_sketch<0:raise ValueError('invalid request')
    print(json.dumps(choose(args.family,args.ani_bin,args.pair_budget,args.target_error,args.queries_per_sketch),indent=2,ensure_ascii=False))
if __name__=='__main__':main()
