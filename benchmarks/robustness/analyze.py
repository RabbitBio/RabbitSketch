#!/usr/bin/env python3
"""Frozen-configuration validation, paired cluster intervals, no selector refitting."""
import collections
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *

NAMES={'minhash':'MinHash','bindash':'BinDash','setsketch':'SetSketch-fixed','fracminhash':'FracMinHash','kssd':'KSSD'}
COLORS=dict(zip(NAMES,['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00']))
METRICS=['j_error','sampling_error','fastani_error','column_error','homologous_error']
def save(name,df):df.to_csv(ANALYSIS/(name+'.csv'),index=False)
def read(name):return pd.read_csv(DATA/(name+'.tsv'),sep='\t',low_memory=False)
def transform(j,k):return np.power(2*np.asarray(j)/(1+np.asarray(j)),1/np.asarray(k))

def load():
    manifest=json.loads((RUN/'suite_manifest.json').read_text());frames=collections.defaultdict(list)
    jobs=manifest['jobs'];audit=[]
    assert len(jobs)==94
    for j in jobs:
        d=RUN/'jobs'/j['job_id'];completion=json.loads((d/'completion.json').read_text())
        assert sha(d/'raw.jsonl')==completion['sha256']
        counts=collections.Counter()
        for line in (d/'raw.jsonl').open():
            if not line.startswith('{'):continue
            r=json.loads(line);typ=r.pop('type');r.update(panel=j['phase'],job_id=j['job_id'],half_subk=j.get('half_subk',6));frames[typ].append(r);counts[typ]+=1
        assert dict(counts)==completion['rows'];audit.append(dict(job_id=j['job_id'],**counts))
    frames={k:pd.DataFrame(v) for k,v in frames.items()}
    o=frames['object'];exact=o[o.algorithm=='exact'].set_index(['panel','k','object_id']).input_digest
    for r in o[o.algorithm!='exact'].itertuples():assert int(r.input_digest)==int(exact.loc[r.panel,r.k,r.object_id])
    atomic(ANALYSIS/'raw_audit.json',dict(jobs=len(jobs),rows={k:len(v) for k,v in frames.items()},suite_sha256=sha(RUN/'suite_manifest.json'),all_raw_hashes_verified=True,normalized_input_digests_match=True))
    save('job_counts',pd.DataFrame(audit));return frames

def prepare(frames):
    truth=frames['pair'];truth=truth[truth.algorithm=='exact'][['panel','pair_id','k','n_a','n_b','intersection']].copy()
    truth['j_true']=truth.intersection/(truth.n_a+truth.n_b-truth.intersection);truth['exact_ani']=transform(truth.j_true,truth.k)
    natural=read('natural_pairs');natural['panel']='natural';natural['cluster_a']=natural.genus_a.astype(str);natural['cluster_b']=natural.genus_b.astype(str)
    mixed=read('mixed_references');mixed['panel']='mixed';mixed['cluster_a']=mixed.reference;mixed['cluster_b']=mixed.reference
    # Metadata owns a/b; exact truth has no endpoints. No many-to-many merges.
    meta=pd.concat([natural,mixed],ignore_index=True)
    base=truth.merge(meta,on=['panel','pair_id'],validate='many_to_one')
    for name,col in [('fastani','fastani'),('column','column_identity'),('homologous','homologous_ungapped_identity')]:base[name+'_model_error']=abs(base.exact_ani-base[col])
    save('exact_truth_and_models',base)
    native=frames['pair'];native=native[native.algorithm!='exact'].drop(columns=['n_a','n_b','intersection'],errors='ignore')
    z=native.merge(base.drop(columns=['a','b']),on=['panel','pair_id','k'],validate='many_to_one')
    z['estimated_ani']=transform(z.jaccard,z.k);z['j_error']=abs(z.jaccard-z.j_true);z['sampling_error']=abs(z.estimated_ani-z.exact_ani)
    for name,col in [('fastani','fastani'),('column','column_identity'),('homologous','homologous_ungapped_identity')]:z[name+'_error']=abs(z.estimated_ani-z[col])
    z['pair_bytes']=z.payload_bytes_a+z.payload_bytes_b;z['zero_j']=(z.jaccard==0).astype(float);z['empty']=((z.empty_a!=0)|(z.empty_b!=0)).astype(float)
    assert z.jaccard.between(0,1).all();assert np.isfinite(z.j_error).all()
    z.to_csv(ANALYSIS/'all_native_rows.csv.gz',index=False)
    # Budget selection uses bytes alone. Require EVERY observed seed to fit.
    sizes=z.groupby(['panel','pair_id','k','algorithm','parameter']).pair_bytes.agg(['mean','max']).reset_index();selected=[]
    for budget in [2048,8192,32768]:
        s=sizes[sizes['max']<=budget].sort_values(['mean','parameter'],ascending=[False,True]).drop_duplicates(['panel','pair_id','k','algorithm'])
        selected.append(z.merge(s[['panel','pair_id','k','algorithm','parameter']],on=['panel','pair_id','k','algorithm','parameter'],validate='many_to_one').assign(budget=budget))
    selected=pd.concat(selected,ignore_index=True);selected.to_csv(ANALYSIS/'budget_selected_rows.csv.gz',index=False)
    return base,z,selected

def metrics(rows,keys,name):
    result=[]
    for group,g in rows.groupby(keys,dropna=False,observed=True):
        if not isinstance(group,tuple):group=(group,)
        base=dict(zip(keys,group));base.update(pairs=g.pair_id.nunique(),seed_observations=len(g),clusters=len(set(g.cluster_a)|set(g.cluster_b)),mean_pair_bytes=g.pair_bytes.mean(),max_pair_bytes=g.pair_bytes.max(),zero_rate=g.zero_j.mean(),empty_rate=g['empty'].mean())
        for m in METRICS:
            valid=g[g[m].notna()]
            if valid.empty:continue
            v=valid[m];perpair=valid.groupby(['pair_id','cluster_a','cluster_b'])[m].mean().reset_index()
            # Within-genus macro only; distant pairs are retained in pair-weighted endpoints.
            within=perpair[perpair.cluster_a==perpair.cluster_b];macro=within.groupby('cluster_a')[m].mean().mean()
            result.append(dict(base,metric=m,valid_pairs=valid.pair_id.nunique(),mae=v.mean(),p95=v.quantile(.95),rmse=np.sqrt((v*v).mean()),within_cluster_macro_mae=macro))
    out=pd.DataFrame(result);save(name,out);return out

def strata(z):
    n=z[z.panel=='natural'];m=z[z.panel=='mixed']
    prior=set(json.loads((ANALYSIS/'data_provenance.json').read_text())['prior_accuracy_genus_taxids'])
    unseen=~n.cluster_a.isin(prior)&~n.cluster_b.isin(prior)
    n=pd.concat([n.assign(stratum='all'),n.assign(stratum=n.ani_bin),n[n.j_true>.5].assign(stratum='J>0.5'),n[unseen].assign(stratum='unseen_accuracy_genera')],ignore_index=True)
    m=m.assign(stratum=m.apply(lambda r:f'id{int(round(r.target_identity*100))}_indel{int(round(r.indel_share*100))}',axis=1))
    return pd.concat([n,m],ignore_index=True)

def intervals(selected):
    rng=np.random.default_rng(202610011);data=strata(selected);out=[]
    keys=['panel','k','budget','stratum']
    for group,g in data.groupby(keys):
        for left,right in [('bindash','minhash'),('setsketch','bindash'),('fracminhash','minhash'),('kssd','bindash')]:
            keep=['pair_id','cluster_a','cluster_b']
            a=g[g.algorithm==left].groupby(keep)[METRICS].mean().reset_index();b=g[g.algorithm==right].groupby(keep)[METRICS].mean().reset_index()
            z=a.merge(b,on=keep,suffixes=('_left','_right'),validate='one_to_one')
            if z.empty:continue
            ids=sorted(set(z.cluster_a)|set(z.cluster_b));idx={v:i for i,v in enumerate(ids)}
            c=rng.multinomial(len(ids),np.ones(len(ids))/len(ids),size=2000)
            ia=z.cluster_a.map(idx).to_numpy();ib=z.cluster_b.map(idx).to_numpy()
            # Same-genus edges get one cluster weight, cross-genus edges a dyadic weight.
            weights=c[:,ia]*np.where(ia==ib,1,c[:,ib])
            for metric in METRICS:
                delta=(z[metric+'_left']-z[metric+'_right']).to_numpy();valid=np.isfinite(delta)
                if not valid.any():continue
                w=weights[:,valid];den=w.sum(axis=1);good=den>0;boots=w[good]@delta[valid]/den[good]
                lo,hi=np.quantile(boots,[.025,.975]);out.append(dict(zip(keys,group),left=left,right=right,metric=metric,pairs=int(valid.sum()),clusters=len(set(z.loc[valid,'cluster_a'])|set(z.loc[valid,'cluster_b'])),difference=delta[valid].mean(),ci_low=lo,ci_high=hi,bootstrap_replicates=len(boots)))
    out=pd.DataFrame(out);save('paired_cluster_intervals',out);return out

def coverage(base,selected):
    nat=read('natural_candidates');names=json.loads((DATA/'taxonomy_names.json').read_text())
    for rank in ['species','genus','family','phylum']:nat[rank+'_name']=nat[rank].astype(str).map(names)
    save('natural_taxonomy',nat)
    save('natural_composition',nat.groupby(['phylum_name','family_name','genus_name'],dropna=False).agg(species=('species','nunique'),genomes=('object_id','size')).reset_index())
    rows=[]
    for panel in ['natural','mixed']:
        for k in [20,24]:
            b=base[(base.panel==panel)&(base.k==k)]
            for budget in [2048,8192,32768]:
                for alg in NAMES:
                    s=selected[(selected.panel==panel)&(selected.k==k)&(selected.budget==budget)&(selected.algorithm==alg)]
                    rows.append(dict(panel=panel,k=k,budget=budget,algorithm=alg,total_pairs=len(b),budget_feasible_pairs=s.pair_id.nunique(),budget_coverage=s.pair_id.nunique()/len(b),empty_observations=s['empty'].sum(),seed_observations=len(s)))
    save('budget_coverage',pd.DataFrame(rows))
    n=base[(base.panel=='natural')&(base.k==24)]
    save('natural_reference_coverage',n.groupby(['ani_bin','pair_class'],dropna=False).agg(pairs=('pair_id','size'),genera_a=('genus_a','nunique'),median_j=('j_true','median'),median_fastani=('fastani','median'),median_af_a=('af_a','median'),median_af_b=('af_b','median')).reset_index())

def costs(frames):
    b=frames['build'].copy();setup=frames.get('shared_setup',pd.DataFrame())
    if len(setup):
        s=setup[['job_id','wall_seconds']].rename(columns={'wall_seconds':'setup_wall_seconds'});b=b.merge(s,on='job_id',how='left',validate='one_to_one')
    else:b['setup_wall_seconds']=0.
    b['setup_wall_seconds']=b.setup_wall_seconds.fillna(0);b['build_with_setup']=b.wall_seconds+b.setup_wall_seconds
    b['ready_wall_seconds']=b.build_with_setup+b.prepare_query_wall_seconds
    inp=frames['input'][['job_id','wall_seconds']].rename(columns={'wall_seconds':'input_wall_seconds'});b=b.merge(inp,on='job_id',validate='one_to_one')
    save('build_batches',b)
    save('build_summary',b.groupby(['panel','k','algorithm','parameter']).agg(replicates=('seed','size'),median_build=('build_with_setup','median'),min_build=('build_with_setup','min'),max_build=('build_with_setup','max'),median_ready=('ready_wall_seconds','median'),median_input=('input_wall_seconds','median'),median_payload=('payload_bytes','median'),peak_rss=('process_peak_rss_bytes','max')).reset_index())
    objects=frames['object'];objects[objects.algorithm!='exact'].to_csv(ANALYSIS/'native_object_costs.csv.gz',index=False)
    p=frames['pair'];p=p[p.algorithm!='exact']
    save('pair_query_cost',p.groupby(['panel','k','algorithm','parameter']).agg(pair_observations=('pair_id','size'),median_j_query_seconds=('query_seconds','median'),mean_j_query_seconds=('query_seconds','mean')).reset_index())

def external_advisor(z):
    source=ROOT/'analysis/regimes_20261001/heldout_selection.csv';old=pd.read_csv(source);out=[]
    for i,r in old.iterrows():
        label=r.ani_bin.replace('–','-')
        # Mash advice transfers to natural genomes; substitution advice to each mixed composition.
        if r.family=='mash142':groups=[('natural',None,z[(z.panel=='natural')&(z.ani_bin==label)])]
        else:
            # Nominal, lower-inclusive controlled bins avoid moving a boundary pair
            # because of one-base rounding; errors still use the achieved identity.
            groups=[]
            for q in [0,.1,.25,.5]:
                lo,hi=map(float,label.split('-'));g=z[(z.panel=='mixed')&(z.indel_share==q)&(z.target_identity>=lo/100)&(z.target_identity<(hi/100 if hi<100 else 1.000001))]
                groups.append(('mixed',q,g))
        for panel,q,g in groups:
            common=dict(request_index=i,source_family=r.family,panel=panel,indel_share=q,ani_bin=label,budget=r.budget,target_mae=r.target_mae,incident_reuse=2*r.Q_charged,algorithm=r.algorithm,parameter=r.parameter,old_validation_pass=r.validation_pass)
            if r.algorithm=='infeasible':out.append(dict(common,status='no_frozen_candidate',pairs=g[g.k==24].pair_id.nunique()));continue
            s=g[(g.k==24)&(g.algorithm==r.algorithm)&(g.parameter==r.parameter)]
            if s.empty:out.append(dict(common,status='no_external_pairs',pairs=0));continue
            error=s.fastani_error if panel=='natural' else s.column_error
            mae=error.mean();p95=error.quantile(.95);budgetpass=s.pair_bytes.max()<=r.budget;passed=mae<=r.target_mae and p95<=2*r.target_mae and budgetpass
            model=s.drop_duplicates('pair_id')['fastani_model_error' if panel=='natural' else 'column_model_error']
            out.append(dict(common,status='pass' if passed else 'fail',pairs=s.pair_id.nunique(),mae=mae,p95=p95,max_pair_bytes=s.pair_bytes.max(),budget_pass=budgetpass,exact_j_model_mae=model.mean(),exact_j_model_p95=model.quantile(.95),exact_j_model_satisfies=int(model.mean()<=r.target_mae and model.quantile(.95)<=2*r.target_mae)))
    save('frozen_advisor_external',pd.DataFrame(out));atomic(ANALYSIS/'frozen_advisor_provenance.json',dict(source=str(source),sha256=sha(source),refit=False,failed_choice_replaced=False,cost_model='unchanged Stage II 10000-genome 64-thread CPU calibration',mixed_truth='achieved alignment-column identity; bin membership uses nominal targets to avoid one-base rounding at boundaries; different estimand from FastANI, explicitly scoped'))

def alignment(z,base):
    result=RUN/'alignment/results.json'
    if not result.exists():return
    a=pd.DataFrame(json.loads(result.read_text()));save('alignment_results',a)
    good=a[(a.status=='ok')&(a.qualified==1)][['panel','pair_id','nucmer']]
    b=base.merge(good,on=['panel','pair_id'],validate='many_to_one');b['nucmer_model_error']=abs(b.exact_ani-b.nucmer);b['fastani_nucmer_difference']=b.fastani-b.nucmer;save('alignment_reference_comparison',b)
    mb=b[b.panel=='mixed']
    save('alignment_mixed_summary',mb.groupby(['k','target_identity','indel_share']).agg(pairs=('pair_id','size'),fastani=('fastani','mean'),nucmer=('nucmer','mean'),column_identity=('column_identity','mean'),homologous_identity=('homologous_ungapped_identity','mean'),nucmer_model_mae=('nucmer_model_error','mean')).reset_index())
    x=z.merge(good,on=['panel','pair_id'],validate='many_to_one');x['nucmer_error']=abs(x.estimated_ani-x.nucmer)
    rows=[]
    for group,g in x.groupby(['panel','k','algorithm','parameter','ani_bin'],dropna=False):
        for metric in ['nucmer_error','fastani_error','column_error','homologous_error','sampling_error']:
            s=g[g[metric].notna()]
            if len(s):rows.append(dict(zip(['panel','k','algorithm','parameter','ani_bin'],group),metric=metric,pairs=s.pair_id.nunique(),mae=s[metric].mean(),p95=s[metric].quantile(.95)))
    save('alignment_sketch_accuracy',pd.DataFrame(rows))

def figures(summary,base):
    plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    def figsave(name,fig):
        fig.tight_layout();fig.savefig(ANALYSIS/(name+'.png'),dpi=180);fig.savefig(ANALYSIS/(name+'.pdf'));plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for ax,metric,title in zip(axes,['j_error','sampling_error','fastani_error'],['Exact Jaccard error','ANI sampling error','ANI error vs FastANI']):
        d=summary[(summary.panel=='natural')&(summary.k==24)&(summary.budget==8192)&(summary.metric==metric)&summary.stratum.isin(['90-95','95-99','99-100'])]
        for alg,g in d.groupby('algorithm'):
            g=g.set_index('stratum').reindex(['90-95','95-99','99-100']);ax.plot(range(3),g.mae*(1 if metric=='j_error' else 100),'o-',label=NAMES[alg],color=COLORS[alg])
        if metric=='fastani_error':
            b=base[(base.panel=='natural')&(base.k==24)].groupby('ani_bin').fastani_model_error.mean().reindex(['90-95','95-99','99-100'])
            ax.plot(range(3),b.values*100,'k--',label='Exact-J model residual');ax.legend(fontsize=8)
        ax.set_xticks(range(3),['90–95','95–99','99–100']);ax.set_xlabel('Bidirectional FastANI (%)');ax.set_title(title);ax.set_ylabel('MAE' if metric=='j_error' else 'MAE (pp)');ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);fig.suptitle('Independent genomes · k=24 · pair payload ≤ 8 KiB');figsave('independent_genome_accuracy',fig)
    fig,axes=plt.subplots(1,4,figsize=(15,3.8),sharey=True)
    for ax,t in zip(axes,[90,95,97,99]):
        d=summary[(summary.panel=='mixed')&(summary.k==24)&(summary.budget==8192)&(summary.metric=='sampling_error')]
        for alg,g in d.groupby('algorithm'):
            g=g.set_index('stratum').reindex([f'id{t}_indel{q}' for q in [0,10,25,50]]);ax.plot([0,10,25,50],g.mae*100,'o-',label=NAMES[alg],color=COLORS[alg])
        ax.set_title(f'Column identity {t}%');ax.set_xlabel('Indel share of edited bases (%)');ax.grid(alpha=.2)
    axes[0].set_ylabel('ANI sampling MAE (pp)');axes[0].legend(fontsize=8);fig.suptitle('Mixed mutation validation · k=24 · pair payload ≤ 8 KiB');figsave('mixed_mutation_sampling',fig)
    m=base[(base.panel=='mixed')&(base.k==24)];fig,axes=plt.subplots(1,3,figsize=(12,4))
    for ax,metric,title in zip(axes,['column_model_error','homologous_model_error','fastani_model_error'],['Known alignment-column truth','Known ungapped homologous truth','FastANI reference']):
        for q,g in m.groupby('indel_share'):
            r=g.groupby('target_identity')[metric].mean();ax.plot(r.index*100,r.values*100,'o-',label=f'indel {q:.0%}')
        ax.set_title(title);ax.set_xlabel('Target column identity (%)');ax.set_ylabel('Exact-J model MAE (pp)');ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);figsave('mixed_mutation_model_residual',fig)

def main():
    frames=load();base,z,selected=prepare(frames)
    metrics(strata(z),['panel','k','algorithm','parameter','stratum'],'configuration_accuracy')
    summary=metrics(strata(selected),['panel','k','budget','algorithm','stratum'],'budget_accuracy')
    coverage(base,selected);costs(frames);intervals(selected);external_advisor(z);alignment(z,base);figures(summary,base)
    m=selected[selected.panel=='mixed'].copy();m['j_band']=pd.cut(m.j_true,[0,.05,.1,.25,.5,1.000001],include_lowest=True,right=False).astype(str)
    metrics(m,['k','budget','algorithm','indel_share','j_band'],'mixed_observed_j_accuracy')
    save('mixed_model_summary',base[base.panel=='mixed'].groupby(['k','target_identity','indel_share']).agg(pairs=('pair_id','size'),mean_j=('j_true','mean'),column_identity=('column_identity','mean'),homologous_identity=('homologous_ungapped_identity','mean'),fastani=('fastani','mean'),column_model_mae=('column_model_error','mean'),column_model_p95=('column_model_error',lambda v:v.quantile(.95)),homologous_model_mae=('homologous_model_error','mean'),fastani_model_mae=('fastani_model_error','mean')).reset_index())
    complete=(RUN/'alignment/completion.json').exists()
    atomic(ANALYSIS/'analysis_status.json',dict(native_pair_observations=len(z),native_statistics_complete=True,alignment_complete=complete))
    print('Native statistical analysis complete',len(z),'observations; alignment complete:',complete,flush=True)

if __name__=='__main__':main()
