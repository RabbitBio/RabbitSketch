#!/usr/bin/env python3
"""Pointwise paired intervals resampling genomes, not correlated pair rows."""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark/analysis/regimes_20261001')
d=pd.read_csv(ROOT/'ani_selected_rows.csv.gz',low_memory=False)
d=d[(d.split=='validation')&d.family.isin(['substitution','mash142'])&d.ani_bin.isin(['80–90','90–95','95–99','99–100'])]
keys=['pair_id','a','b','family','reference','k','budget','algorithm','ani_bin']
d=d.groupby(keys,dropna=False,observed=True)[['error_native','sampling_error_native']].mean().reset_index()
rng=np.random.default_rng(20261003);rows=[]
for group,g in d.groupby(['family','k','budget','ani_bin']):
    for left,right in [('bindash','minhash'),('setsketch','bindash')]:
        a=g[g.algorithm==left];b=g[g.algorithm==right]
        z=a.merge(b,on='pair_id',suffixes=('_left','_right'),validate='one_to_one')
        if z.empty:continue
        if group[0]=='substitution':
            ids=sorted(z.a_left.unique());idx={n:i for i,n in enumerate(ids)};counts=rng.multinomial(len(ids),np.ones(len(ids))/len(ids),size=1000);weights=counts[:,z.a_left.map(idx)]
        else:
            ids=sorted(set(z.a_left)|set(z.b_left));idx={n:i for i,n in enumerate(ids)};counts=rng.multinomial(len(ids),np.ones(len(ids))/len(ids),size=1000);weights=counts[:,z.a_left.map(idx)]*counts[:,z.b_left.map(idx)]
        sums=weights.sum(axis=1);weights=weights[sums>0];sums=sums[sums>0]
        for metric in ['error_native','sampling_error_native']:
            delta=(z[metric+'_left']-z[metric+'_right']).to_numpy();boots=weights@delta/sums
            lo,hi=np.quantile(boots,[.025,.975]);rows.append(dict(family=group[0],k=group[1],budget=group[2],ani_bin=group[3],left=left,right=right,metric=metric,pairs=len(z),genomes=len(ids),difference=delta.mean(),ci_low=lo,ci_high=hi,usable_bootstraps=len(boots),scope='pointwise genome-resampling interval conditional on three tested hash seeds; no multiplicity adjustment'))
pd.DataFrame(rows).to_csv(ROOT/'paired_ani_intervals.csv',index=False)
print('Paired intervals',len(rows))
