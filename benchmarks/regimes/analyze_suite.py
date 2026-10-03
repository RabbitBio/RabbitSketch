#!/usr/bin/env python3
"""Task-specific summaries, independent-reference errors and held-out selection."""
from pathlib import Path
import argparse, csv, json, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark');DATA=ROOT/'data/regimes_20261001';RUN=ROOT/'runs/regimes_20261001';OUT=ROOT/'analysis/regimes_20261001'
OUT.mkdir(parents=True,exist_ok=True)
NAMES=dict(minhash='MinHash',bindash='BinDash',bindash_hll='BinDash+HLL',setsketch='SetSketch-fixed',fracminhash='FracMinHash',kssd='KSSD',hll='HLL',probminhash='ProbMinHash')
COLORS=dict(minhash='#0072B2',bindash='#D55E00',bindash_hll='#882255',setsketch='#009E73',fracminhash='#CC79A7',kssd='#E69F00',hll='#56B4E9',probminhash='#332288')
INCLUDE_SUPPLEMENT=True
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':130,'pdf.fonttype':42})
def read(name):
    directory=ROOT/'data/containment_holdout_20261001' if name.startswith('contain_') else DATA
    return pd.read_csv(directory/(name+'.tsv'),sep='\t')
def save(name,frame):frame.to_csv(OUT/(name+'.csv'),index=False)
def figure(name,fig):fig.tight_layout();fig.savefig(OUT/(name+'.png'),dpi=180);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
def native(j,k):
    j=np.asarray(j);return np.power(np.divide(2*j,1+j),1/np.asarray(k))
def mash(j,k):
    j=np.asarray(j);safe=np.maximum(j,np.finfo(float).tiny);return np.where(j>0,np.maximum(0,1+np.log(2*safe/(1+safe))/np.asarray(k)),0)
def group_stats(d,keys,value):
    return d.groupby(keys,dropna=False,observed=True).agg(n=(value,'size'),valid=(value,'count'),mae=(value,'mean'),p95=(value,lambda x:x.quantile(.95)),rmse=(value,lambda x:np.sqrt((x*x).mean()))).reset_index()
def budget_select(d):
    keys=['pair_id','k','algorithm','parameter']
    sizes=d.assign(bytes=d.payload_bytes_a+d.payload_bytes_b).groupby(keys).bytes.agg(['mean','max']).reset_index()
    result=[]
    for budget in [2048,8192,32768]:
        selected=sizes[sizes['max']<=budget].sort_values(['mean','parameter'],ascending=[False,True]).drop_duplicates(['pair_id','k','algorithm'])
        result.append(d.merge(selected[keys],on=keys,validate='many_to_one').assign(budget=budget))
    return pd.concat(result,ignore_index=True)
def load_phase(phase,allow_partial=False):
    roots=[RUN]
    supplement=ROOT/'runs/regimes_supplement_20261001'
    if INCLUDE_SUPPLEMENT and (supplement/'suite_manifest.json').exists():roots.append(supplement)
    compact=ROOT/'runs/compressed_validation_20261001'
    if INCLUDE_SUPPLEMENT and (compact/'suite_manifest.json').exists():roots.append(compact)
    if phase=='containment':roots=[ROOT/'runs/containment_holdout_20261001']
    jobs=[(root,j) for root in roots for j in json.loads((root/'suite_manifest.json').read_text())['jobs'] if j['phase']==phase]
    pending=[j for root,j in jobs if not (root/'jobs'/j['job_id']/'completion.json').exists()]
    if pending and not allow_partial:raise RuntimeError(f'{phase}: {len(pending)} jobs pending')
    rows={};expected={};counts={}
    for root,j in jobs:
        path=root/'jobs'/j['job_id']/'raw.jsonl'
        if not path.exists():continue
        for line in path.open():
            if not line.startswith('{'):continue
            r=json.loads(line);typ=r.pop('type');r.update(job_id=j['job_id'],phase=phase)
            rows.setdefault(typ,[]).append(r)
        counts[j['job_id']]=json.loads((path.parent/'completion.json').read_text())['rows']
    frames={k:pd.DataFrame(v) for k,v in rows.items()}
    obj=frames.get('object',pd.DataFrame())
    if len(obj) and phase not in ['weighted','workload']:
        exact=obj[obj.algorithm=='exact'].set_index(['object_id','k']).input_digest
        for _,r in obj[obj.algorithm!='exact'][['object_id','k','input_digest']].iterrows():
            if int(r.input_digest)!=int(exact.loc[r.object_id,r.k]):raise AssertionError('normalized FASTA digest mismatch')
    (OUT/(phase+'_audit.json')).write_text(json.dumps(dict(planned_jobs=len(jobs),completed_jobs=len(counts),pending=len(pending),row_counts={k:len(v) for k,v in frames.items()}),indent=2)+'\n')
    return frames
def ani(frames):
    d=frames['pair'];truth=d[d.algorithm=='exact'][['pair_id','k','n_a','n_b','intersection']].copy();truth['j_true']=truth.intersection/(truth.n_a+truth.n_b-truth.intersection)
    meta=read('ani_pairs');objects=read('ani_objects').set_index('object_id').path.to_dict()
    fa={}
    for line in (RUN/'fastani/references.tsv').open():
        a,b,ident,m,n=line.strip().split('\t');fa[a,b]=(float(ident)/100,int(m)/int(n))
    identities=[]
    for r in meta.to_dict('records'):
        if r['family']=='substitution':r.update(identity_source='known_positional_substitution',reference_status='available',ani_ab=r['identity'],ani_ba=r['identity'],af_ab=1.,af_ba=1.)
        else:
            a,b=objects[r['a']],objects[r['b']];x=fa.get((a,b),(np.nan,np.nan));y=fa.get((b,a),(np.nan,np.nan));ok=np.isfinite(x[0]) and np.isfinite(y[0]) and min(x[1],y[1])>=.2
            r.update(identity_source='FastANI_1.34_bidirectional_mean',reference_status='available' if ok else 'missing_or_low_AF',identity=(x[0]+y[0])/2 if ok else np.nan,ani_ab=x[0],ani_ba=y[0],af_ab=x[1],af_ba=y[1])
        identities.append(r)
    meta=pd.DataFrame(identities);save('ani_references',meta)
    meta.loc[meta.identity<.8,'reference_status']='below_80_outside_primary_range'
    meta['ani_bin']=pd.cut(meta.identity,[0,.8-1e-9,.9,.95,.99,1+1e-9],right=False,labels=['below80','80–90','90–95','95–99','99–100'])
    save('ani_references',meta)
    base=truth.merge(meta,on='pair_id',validate='many_to_one');base['exact_native_ani']=native(base.j_true,base.k);base['exact_mash_ani']=mash(base.j_true,base.k)
    base['model_error_native']=abs(base.exact_native_ani-base.identity);base['model_error_mash']=abs(base.exact_mash_ani-base.identity)
    save('ani_exact_model_errors',base)
    z=d[d.algorithm!='exact'].merge(base.drop(columns=['a','b']),on=['pair_id','k'],validate='many_to_one')
    z['native_ani']=native(z.jaccard,z.k);z['mash_ani']=mash(z.jaccard,z.k)
    z['error_native']=abs(z.native_ani-z.identity);z['error_mash']=abs(z.mash_ani-z.identity);z['sampling_error_native']=abs(z.native_ani-z.exact_native_ani)
    z['bias_native']=z.native_ani-z.identity;z['zero_jaccard']=(z.jaccard==0).astype(float)
    selected=budget_select(z)
    selected[['pair_id','a','b','family','reference','split','k','budget','algorithm','parameter','ani_bin','identity','native_ani','mash_ani','exact_native_ani','exact_mash_ani','error_native','sampling_error_native']].to_csv(OUT/'ani_selected_rows.csv.gz',index=False)
    keys=['family','split','k','budget','algorithm','ani_bin']
    summaries=[]
    for metric in ['error_native','error_mash','sampling_error_native']:
        summaries.append(group_stats(selected,keys,metric).assign(metric=metric))
    summary=pd.concat(summaries,ignore_index=True);save('ani_budget_accuracy',summary)
    save('ani_config_accuracy',group_stats(z,['family','split','k','algorithm','parameter','ani_bin'],'error_native'))
    save('ani_zero_rate',selected.groupby(keys,observed=True).zero_jaccard.mean().reset_index())
    threshold=[]
    for tau,width in [(.95,.01),(.99,.005)]:
        s=selected[(selected.identity>=tau-width)&(selected.identity<=tau+width)].copy();s['truth_positive']=s.identity>=tau-1e-12;s['prediction_positive']=s.native_ani>=tau
        s['misclassified']=(s.truth_positive!=s.prediction_positive).astype(float);s['false_positive']=((~s.truth_positive)&s.prediction_positive).astype(float);s['false_negative']=(s.truth_positive&(~s.prediction_positive)).astype(float)
        for label,subset in [('all_in_window',s),('exclude_exact_threshold',s[abs(s.identity-tau)>1e-6])]:
            threshold.append(subset.groupby(['family','split','k','budget','algorithm']).agg(n=('pair_id','size'),misclassification_rate=('misclassified','mean'),false_positive_fraction=('false_positive','mean'),false_negative_fraction=('false_negative','mean')).reset_index().assign(threshold=tau,window=width,evaluation=label))
    save('ani_threshold_errors',pd.concat(threshold,ignore_index=True))
    # Store pair means compactly for later held-out, task-constrained decisions.
    agg=z.groupby(['pair_id','family','reference','split','k','algorithm','parameter','ani_bin'],dropna=False,observed=True).agg(error=('error_native','mean'),sampling_error=('sampling_error_native','mean'),bytes_a=('payload_bytes_a','max'),bytes_b=('payload_bytes_b','max')).reset_index()
    agg['max_pair_bytes']=agg.bytes_a+agg.bytes_b;save('ani_pair_config',agg)
    fig,axs=plt.subplots(1,3,figsize=(14,4))
    for ax,family in zip(axs,['substitution','mash142','diverse200']):
        s=summary[(summary.family==family)&(summary.split=='validation')&(summary.k==24)&(summary.budget==8192)&(summary.metric=='error_native')&summary.ani_bin.isin(['80–90','90–95','95–99','99–100'])]
        for alg,g in s.groupby('algorithm',observed=True):
            x=[['80–90','90–95','95–99','99–100'].index(b) for b in g.ani_bin];ax.plot(x,100*g.mae,'o-',label=NAMES[alg],color=COLORS[alg])
        b=base[(base.family==family)&(base.split=='validation')&(base.k==24)&base.ani_bin.isin(['80–90','90–95','95–99','99–100'])].groupby('ani_bin',observed=True).model_error_native.mean()
        ax.plot([['80–90','90–95','95–99','99–100'].index(i) for i in b.index],100*b.values,'k--',label='Exact-J model error')
        ax.set_xticks(range(4),['80–90','90–95','95–99','99–100']);ax.set_xlabel('Reference identity (%)');ax.set_ylabel('ANI MAE (percentage points)');ax.set_title(family);ax.set_yscale('log');ax.grid(alpha=.2)
    axs[0].legend(fontsize=8);figure('ani_error_decomposition',fig)
    save('ani_model_error_summary',pd.concat([group_stats(base,['family','split','k','ani_bin'],m).assign(metric=m) for m in ['model_error_native','model_error_mash']],ignore_index=True))
    return selected
def cardinality(frames):
    o=frames['object'];truth=o[o.algorithm=='exact'][['object_id','k','n']];z=o[o.algorithm!='exact'].drop(columns='n').merge(truth,on=['object_id','k'],validate='many_to_one').merge(read('card_meta'),on='object_id',validate='many_to_one')
    z['relative_error']=abs(z.cardinality-z.n)/z.n;z['relative_bias']=(z.cardinality-z.n)/z.n
    save('cardinality_objects',z);summary=group_stats(z,['algorithm','parameter','k','length','kind'],'relative_error');save('cardinality_accuracy',summary)
    save('cardinality_build',z.groupby(['algorithm','parameter','k','length','kind']).agg(median_build_seconds=('build_seconds','median'),median_cardinality_seconds=('cardinality_seconds','median'),payload_bytes=('payload_bytes','median')).reset_index())
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for ax,length in zip(axs,[5_000_000,50_000_000]):
        s=summary[(summary.k==24)&(summary.length==length)]
        for alg,g in s.groupby('algorithm'):ax.plot(2.**g.parameter,100*g.mae,'o-',color=COLORS[alg],label=NAMES[alg])
        ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xlabel('Single-object register bytes');ax.set_ylabel('Mean cardinality error (%)');ax.set_title(f'{length/1e6:g} Mb');ax.grid(alpha=.2);ax.legend()
    figure('cardinality_tradeoff',fig)
def containment(frames):
    d=frames['pair'];truth=d[d.algorithm=='exact'][['pair_id','k','n_a','n_b','intersection']].copy();truth['truth_ab']=truth.intersection/truth.n_a;truth['truth_ba']=truth.intersection/truth.n_b
    z=d[d.algorithm!='exact'].merge(truth,on=['pair_id','k'],validate='many_to_one').merge(read('contain_pairs').drop(columns=['a','b']),on='pair_id',validate='many_to_one')
    z['available']=(z.empty_a==0)&(z.empty_b==0)&np.isfinite(z.containment_derived_ab)
    # Match UnifiedQuery's common intersection cap, including the reverse
    # direction: independently clipping each C would be inconsistent.
    cap=np.maximum(1,np.maximum(z.containment_derived_ab,z.containment_derived_ba))
    z['derived_bounded']=np.maximum(0,z.containment_derived_ab/cap)
    z['derived_bounded_reverse']=np.maximum(0,z.containment_derived_ba/cap)
    z['error']=np.where(z.available,abs(z.derived_bounded-z.truth_ab),np.nan)
    z['native_error']=np.where(z.available,abs(z.containment_native_ab-z.truth_ab),np.nan)
    # KSSD has no public native directed-containment query in this revision.
    z.loc[z.algorithm=='kssd','native_error']=np.nan
    z['native_out_of_range']=(~np.isfinite(z.containment_native_ab))|(z.containment_native_ab<0)|(z.containment_native_ab>1)
    selected=budget_select(z);keys=['split','k','budget','algorithm','nominal_ratio','nominal_containment']
    summ=group_stats(selected,keys,'error');avail=selected.groupby(keys).agg(availability=('available','mean'),pair_count=('pair_id','nunique'),native_invalid_rate=('native_out_of_range','mean')).reset_index();summ=summ.merge(avail,on=keys)
    save('containment_budget_accuracy',summ);save('containment_native_accuracy',group_stats(selected,keys,'native_error'));save('containment_exact_truth',truth.merge(read('contain_pairs'),on='pair_id'))
    save('containment_config_accuracy',group_stats(z,['split','k','algorithm','parameter','nominal_ratio','nominal_containment'],'error'))
    save('containment_query_cost',selected.groupby(keys).agg(median_joint_query_seconds=('query_seconds','median')).reset_index())
    reverse=selected.copy();reverse['error_reverse']=np.where(reverse.available,abs(reverse.derived_bounded_reverse-reverse.truth_ba),np.nan);save('containment_reverse_accuracy',group_stats(reverse,keys,'error_reverse'))
    # Empirical uncertainty among seeds is nested in biological reference IDs.
    selected[['pair_id','reference','split','k','algorithm','parameter','budget','nominal_ratio','nominal_containment','error','native_error','available']].to_csv(OUT/'containment_selected_rows.csv.gz',index=False)
    plot_containment_summary(summ)
def plot_containment_summary(summ):
    fig,axs=plt.subplots(2,3,figsize=(14,7),sharex='col',height_ratios=[2,1])
    for col,c in enumerate([1.,.5,.1]):
        s=summ[(summ.split=='validation')&(summ.k==24)&(summ.budget==8192)&(summ.nominal_containment==c)]
        for alg,g in s.groupby('algorithm'):
            g=g.sort_values('nominal_ratio')
            label=NAMES[alg]+(' (experimental)' if alg=='kssd' else '')
            style='o--' if alg=='kssd' else 'o-'
            axs[0,col].plot(g.nominal_ratio,g.mae,style,color=COLORS[alg],label=label)
            axs[1,col].plot(g.nominal_ratio,g.availability,style,color=COLORS[alg])
        axs[0,col].set_ylabel('MAE among available outputs')
        axs[0,col].set_title(f'Nominal containment = {c:g}')
        axs[1,col].set_ylabel('Available output fraction')
        axs[1,col].set_ylim(-.05,1.08)
        axs[1,col].set_xlabel('Reference / query length')
        for ax in axs[:,col]:ax.set_xscale('log');ax.grid(alpha=.2)
    axs[0,0].legend(fontsize=8)
    fig.suptitle('Held-out containment: k=24, 8 KiB pair payload budget; empty samples are unavailable',fontsize=11)
    figure('containment_boundaries',fig)
def workload(frames):
    b=frames['build'].copy();q=frames['allpairs'].copy();input_time=frames['input'].set_index('job_id').wall_seconds
    setup=frames.get('shared_setup',pd.DataFrame(columns=['job_id','wall_seconds'])).set_index('job_id').wall_seconds
    qstat=q.groupby('job_id').agg(query_wall_median=('wall_seconds','median'),query_wall_min=('wall_seconds','min'),query_wall_max=('wall_seconds','max'),query_pairs=('pairs','first'),checksum_min=('checksum','min'),checksum_max=('checksum','max')).reset_index()
    z=b.merge(qstat,on='job_id',validate='one_to_one');z['setup_seconds']=z.job_id.map(setup).fillna(0);z['input_seconds']=z.job_id.map(input_time);z['build_total']=z.wall_seconds+z.prepare_query_wall_seconds+z.setup_seconds;z['total_seconds']=z.build_total+z.query_wall_median;z['total_with_input']=z.total_seconds+z.input_seconds
    z['B']=z.build_total/z.objects;z['t_pair']=z.query_wall_median/z.query_pairs;z['query_pairs_per_second']=z.query_pairs/z.query_wall_median;save('workload_measured',z)
    cross=[]
    for n,g in z.groupby('objects'):
        for pm,pb in [(128,512),(512,2048),(2048,8192)]:
            if g[(g.algorithm=='minhash')&(g.parameter==pm)].empty or g[(g.algorithm=='bindash')&(g.parameter==pb)].empty:continue
            ma=g[(g.algorithm=='minhash')&(g.parameter==pm)].iloc[0];bd=g[(g.algorithm=='bindash')&(g.parameter==pb)].iloc[0];den=ma.t_pair-bd.t_pair
            qs=(bd.B-ma.B)/den if den else np.nan;valid=np.isfinite(qs) and qs>0
            cross.append(dict(objects=n,minhash_parameter=pm,bindash_parameter=pb,B_minhash=ma.B,B_bindash=bd.B,t_minhash=ma.t_pair,t_bindash=bd.t_pair,Q_charged_crossover=qs if valid else np.nan,incident_reuse_crossover=2*qs if valid else np.nan,allpairs_N_crossover=2*qs+1 if valid else np.nan,positive_crossover=valid))
    save('workload_crossovers',pd.DataFrame(cross))
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    for (alg,p),g in z.groupby(['algorithm','parameter']):
        if (alg,int(p)) not in [('minhash',512),('bindash',2048),('setsketch',12),('fracminhash',10000),('kssd',3),('kssd',4)]:continue
        g=g.sort_values('objects');label=f'{NAMES[alg]} ({p})';style='--' if alg=='kssd' and p==4 else '-';axs[0].plot(g.objects,g.total_seconds,marker='o',label=label,color=COLORS[alg],linestyle=style)
        axs[1].plot(g.objects,g.query_pairs_per_second,marker='o',linestyle=style,label=label,color=COLORS[alg])
    for ax in axs:ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Distinct genomes');ax.grid(alpha=.2)
    axs[0].set_ylabel('Build + prepare + all-pairs seconds');axs[1].set_ylabel('Pair comparisons / second');axs[0].legend(fontsize=7);figure('workload_scaling',fig)
    return z
def weighted(frames):
    z=frames['pair'].merge(read('weighted_pairs').drop(columns=['a','b']),on='pair_id',validate='many_to_one');z['error']=abs(z.jaccard-z.probability_jaccard);z['error_wrong_metric']=abs(z.jaccard-z.minmax_jaccard)
    keys=['family','model','perturbation','split','parameter'];summary=group_stats(z,keys,'error');save('weighted_accuracy',summary);save('weighted_wrong_metric',group_stats(z,keys,'error_wrong_metric'))
    save('weighted_pair_accuracy',z.groupby(['pair_id','family','model','perturbation','parameter']).agg(mae=('error','mean'),p95=('error',lambda x:x.quantile(.95)),estimate=('jaccard','mean'),truth=('probability_jaccard','first'),support_jaccard=('support_jaccard','first'),minmax_jaccard=('minmax_jaccard','first')).reset_index())
    fig,axs=plt.subplots(1,2,figsize=(10,4));s=z[(z.parameter==2048)&(z.family=='synthetic_weights')&(z.perturbation!=10)].groupby('pair_id').agg(truth=('probability_jaccard','first'),estimate=('jaccard','mean'))
    axs[0].scatter(s.truth,s.estimate,s=8,alpha=.5,color=COLORS['probminhash']);axs[0].plot([0,1],[0,1],'k--');axs[0].set_xlabel('Exact probability Jaccard');axs[0].set_ylabel('Mean sketch estimate (20 seeds)');axs[0].set_title('Synthetic abundance, m=2048')
    s=z[z.family=='cami_marine_profile'].groupby('parameter').error.agg(['mean',lambda x:x.quantile(.95)]);axs[1].plot(s.index*16,s['mean'],'o-',label='MAE');axs[1].plot(s.index*16,s.iloc[:,1],'s--',label='P95 error');axs[1].set_xscale('log',base=2);axs[1].set_yscale('log');axs[1].set_xlabel('Pair-total winner bytes');axs[1].set_ylabel('Absolute probability-Jaccard error');axs[1].set_title('CAMI II Marine profiles');axs[1].legend();figure('weighted_validation',fig)
def selection():
    if not (OUT/'ani_pair_config.csv').exists() or not (OUT/'workload_measured.csv').exists():return
    d=pd.read_csv(OUT/'ani_pair_config.csv',low_memory=False);accuracy=pd.read_csv(OUT/'ani_config_accuracy.csv');cost=pd.read_csv(OUT/'workload_measured.csv');cost=cost[cost.objects==10000]
    rng=np.random.default_rng(95199);rows=[]
    # Select on controlled calibration references; verify on disjoint references.
    for family in ['substitution','mash142']:
      s=d[(d.family==family)&(d.k==24)]
      for abin in ['80–90','90–95','95–99','99–100']:
       for budget in [8192,32768]:
        for eps in [.001,.005]:
         for Q in [.5,5000.,100000.]:
          candidates=[]
          for (alg,p),g in s[s.ani_bin==abin].groupby(['algorithm','parameter']):
           cc=cost[(cost.algorithm==alg)&(cost.parameter==p)]
           cal=g[g.split=='calibration'];val=g[g.split=='validation']
           if cc.empty or cal.empty or val.empty or cal.max_pair_bytes.max()>budget:continue
           scores=accuracy[(accuracy.family==family)&(accuracy.k==24)&(accuracy.ani_bin==abin)&(accuracy.algorithm==alg)&(accuracy.parameter==p)]
           cg=scores[scores.split=='calibration'].iloc[0];vg=scores[scores.split=='validation'].iloc[0]
           e=cg.mae;p95=cg.p95
           # No independence claim for pairwise natural genomes: descriptive gate.
           if e>eps or p95>2*eps:continue
           c=cc.iloc[0];candidates.append(dict(algorithm=alg,parameter=p,calibration_mae=e,calibration_p95=p95,validation_mae=vg.mae,validation_p95=vg.p95,validation_budget_coverage=float((val.max_pair_bytes<=budget).mean()),validation_pass=bool(vg.mae<=eps and vg.p95<=2*eps and val.max_pair_bytes.max()<=budget),cost=c.B+Q*c.t_pair))
          base=dict(family=family,ani_bin=abin,budget=budget,target_mae=eps,Q_charged=Q,eligible_configurations=len(candidates))
          rows.append(dict(**base,**min(candidates,key=lambda r:r['cost'])) if candidates else dict(**base,algorithm='infeasible'))
    save('heldout_selection',pd.DataFrame(rows))
    result=pd.DataFrame(rows)
    from matplotlib.colors import ListedColormap
    labels=['infeasible','minhash','bindash','setsketch','fracminhash','kssd'];colors=['#EEEEEE']+[COLORS[a] for a in labels[1:]]
    for family in ['substitution','mash142']:
        fig,axs=plt.subplots(2,2,figsize=(11,6))
        for ax,(budget,eps) in zip(axs.flat,[(8192,.001),(32768,.001),(8192,.005),(32768,.005)]):
            s=result[(result.family==family)&(result.budget==budget)&(result.target_mae==eps)]
            arr=np.zeros((3,4));cells={}
            for iy,Q in enumerate([.5,5000.,100000.]):
                for ix,abin in enumerate(['80–90','90–95','95–99','99–100']):
                    r=s[(s.Q_charged==Q)&(s.ani_bin==abin)].iloc[0];arr[iy,ix]=labels.index(r.algorithm);cells[iy,ix]=r
            ax.imshow(arr,cmap=ListedColormap(colors),vmin=0,vmax=len(labels)-1,aspect='auto')
            for (iy,ix),r in cells.items():
                text='No fit' if r.algorithm=='infeasible' else NAMES[r.algorithm]+f'\n{r.parameter:g}'
                if r.algorithm!='infeasible' and not bool(r.validation_pass):text+=' ×'
                ax.text(ix,iy,text,ha='center',va='center',fontsize=8,color='black' if r.algorithm in ['infeasible','kssd'] else 'white')
            ax.set_xticks(range(4),['80–90','90–95','95–99','99–100']);ax.set_yticks(range(3),['1','10,000','200,000']);ax.set_xlabel('Reference ANI (%)');ax.set_ylabel('Incident reuse per sketch');ax.set_title(f'{budget//1024} KiB / pair; MAE ≤ {100*eps:g} pp; P95 ≤ {200*eps:g} pp')
        fig.suptitle(f'{family}: calibration choice, independent holdout check (× = failed)',fontsize=11);figure('selection_'+family,fig)
def main():
    global INCLUDE_SUPPLEMENT
    ap=argparse.ArgumentParser();ap.add_argument('--phases',nargs='+',default=['ani','cardinality','containment','workload','weighted']);ap.add_argument('--primary-only',action='store_true');args=ap.parse_args();INCLUDE_SUPPLEMENT=not args.primary_only
    for phase in args.phases:
        frames=load_phase(phase);print(phase,{k:len(v) for k,v in frames.items()},flush=True)
        if 'build' in frames:save(phase+'_build_batches',frames['build'])
        globals()[phase](frames)
    selection()
if __name__=='__main__':main()
