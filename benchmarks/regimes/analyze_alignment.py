#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark');OUT=ROOT/'analysis/regimes_20261001'
raw=pd.DataFrame(json.loads((ROOT/'runs/regimes_20261001/alignment_audit/results.json').read_text()))
raw['nucmer_ani']=(raw.nucmer_ref_identity+raw.nucmer_query_identity)/2
raw['eligible']=(raw.one_to_one_length_fraction_ref>=.2)&(raw.one_to_one_length_fraction_query>=.2)
raw['reference_difference']=raw.fastani-raw.nucmer_ani
raw.to_csv(OUT/'alignment_references.csv',index=False)
d=pd.read_csv(OUT/'ani_selected_rows.csv.gz',low_memory=False)
z=d.merge(raw[['pair_id','nucmer_ani','eligible']],on='pair_id',validate='many_to_one');z=z[z.eligible]
z['nucmer_bin']=pd.cut(z.nucmer_ani,[0,.8,.9,.95,.99,1.00001],right=False,labels=['below80','80–90','90–95','95–99','99–100'])
z['error_vs_nucmer']=abs(z.native_ani-z.nucmer_ani);z['exact_model_error_vs_nucmer']=abs(z.exact_native_ani-z.nucmer_ani)
summary=z.groupby(['family','k','budget','algorithm','nucmer_bin'],observed=True).agg(n=('pair_id','size'),pairs=('pair_id','nunique'),mae_vs_nucmer=('error_vs_nucmer','mean'),p95_vs_nucmer=('error_vs_nucmer',lambda x:x.quantile(.95)),mae_vs_fastani=('error_native','mean'),model_mae_vs_nucmer=('exact_model_error_vs_nucmer','mean')).reset_index();summary.to_csv(OUT/'alignment_sketch_accuracy.csv',index=False)
base=z.drop_duplicates(['pair_id','k']);base[['pair_id','family','k','nucmer_ani','identity','nucmer_bin','exact_model_error_vs_nucmer']].to_csv(OUT/'alignment_exact_model_errors.csv',index=False)
result=dict(audited_pairs=len(raw),eligible_pairs=int(raw.eligible.sum()),mean_absolute_reference_difference=float(abs(raw.reference_difference).mean()),p95_absolute_reference_difference=float(abs(raw.reference_difference).quantile(.95)),max_absolute_reference_difference=float(abs(raw.reference_difference).max()),rounding_percentage_points=.0001)
(OUT/'alignment_audit_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
