#!/usr/bin/env python3
"""Standalone scientific plots for observable selection and its failure modes."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

from scenario_common import *

SHORT = dict(minhash='FKMV',fracminhash='Frac',hll='HLL',setsketch='SS*',kssd='KSSD',
             probminhash_unweighted='PMH',bindash='BD')
COLORS = dict(minhash='#0072B2',fracminhash='#CC79A7',hll='#56B4E9',setsketch='#009E73',
              kssd='#E69F00',probminhash_unweighted='#332288',bindash='#D55E00')
POLICY_COLORS = dict(pilot='#0072B2',metadata='#CC79A7',fixed='#333333')
POLICY_NAMES = dict(pilot='Pilot-based rule',metadata='Metadata-only rule',fixed='Development-fixed best')
METRICS = ['jaccard','containment','cardinality','neighbor','species','clustering']
METRIC_NAMES = dict(jaccard='Jaccard',containment='Query containment',cardinality='Query cardinality',
                    neighbor='Exact-J neighbor recovery',species='Source-species retrieval',clustering='Query clustering')
CONDITIONS = CONFIG['conditions'] + json.loads((HERE/'validation_extra_conditions.json').read_text())
ORDER = [c['name'] for c in CONDITIONS]
LABELS = dict(full='Full genomes',fragment_100k='100 kb target fragment',fragment_10k='10 kb fragment',
    fragment_1k='1 kb fragment',incomplete_10pct='10% assembly',
    incomplete_10pct_contaminated='10% assembly + 20% contaminant',reads_cov01_clean='150 bp / 0.1x / no error',
    reads_cov01_error01='150 bp / 0.1x / 1% error',reads_cov1_error01='150 bp / 1x / 1% error',
    reads_cov1_error03='150 bp / 1x / 3% error',unseen_fragment_5k='NEW: 5 kb fragment',
    unseen_incomplete_30pct_contam10='NEW: 30% assembly + 10% contaminant',
    unseen_reads_cov05_error02='NEW: 150 bp / 0.5x / 2% error',
    unseen_reads250_cov05_error01='NEW: 250 bp / 0.5x / 1% error')


def save(fig,name):
    directory=OUT/'figures';directory.mkdir(parents=True,exist_ok=True)
    fig.savefig(directory/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(directory/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)


def candidates():
    d=pd.read_csv(OUT/'scenario_candidates.csv').fillna('')
    d=d[(d.k==24)&(d.budget==8192)].set_index(['condition','metric'])
    fig,ax=plt.subplots(figsize=(13.5,8.8))
    for i,condition in enumerate(ORDER):
        for j,metric in enumerate(METRICS):
            row=d.loc[(condition,metric),:]; algs=row.candidates.split(';') if row.candidates else []
            no_feasible=row.feasible_alternative_rate==0
            color='#888888' if no_feasible else COLORS.get(row.best_observed,'#dddddd')
            ax.add_patch(plt.Rectangle((j-.5,i-.5),1,1,facecolor=color,alpha=.16,edgecolor='white'))
            labels=[SHORT[a] for a in algs]
            if len(labels)>4:label='/'.join(labels[:3])+f'\n+ {len(labels)-3} candidates'
            else:label=' / '.join(labels[:2])+ ('\n'+' / '.join(labels[2:]) if len(labels)>2 else '')
            if no_feasible:label='None meet target\nLowest loss: '+SHORT.get(row.best_observed,'-')
            elif row.feasible_alternative_rate<1:label+=' †'
            ax.text(j,i,label or 'Unavailable',ha='center',va='center',fontsize=8)
    ax.set(xlim=(-.5,5.5),ylim=(len(ORDER)-.5,-.5),xticks=range(6),yticks=range(len(ORDER)))
    ax.set_xticklabels(['Jaccard','Containment','Cardinality','NN recovery','Species','Clustering'])
    ax.set_yticklabels([LABELS[c] for c in ORDER],fontsize=9)
    ax.xaxis.tick_top();ax.tick_params(length=0)
    ax.axhline(9.5,color='black',lw=1)
    ax.set_title('Held-out descriptive candidates | k=24 | final pair payload budget: 8 KiB',pad=36)
    fig.text(.02,-.025,'30 genera held out from development. Contrasts are simultaneous within each cell only. SS*: corrected SetSketch.\n'
        'Grey: no measured configuration meets the standard target in any workload. †: some workloads have no feasible configuration.\n'
        'Candidate membership describes relative loss, not target attainment. Full candidates and feasibility fractions are in CSV.\n'
        'Targets: J MAE/P95 .02/.05; C .05/.10; relative cardinality .10/.25; NN >= .90; species >= .95; ARI >= .80; validity >= .99.',fontsize=9)
    save(fig,'scenario_candidates_k24_8k')


def primary_intervals():
    d=pd.read_csv(OUT/'advisor_vs_fixed_intervals.csv')
    d=d[(d.condition=='all')&(d.policy=='pilot')]
    fig,axes=plt.subplots(2,3,figsize=(14,8))
    colors=['#0072B2','#009E73','#D55E00']
    for ax,metric in zip(axes.flat,METRICS):
        for offset,(k,color) in enumerate(zip(CONFIG['k'],colors)):
            g=d[(d.metric==metric)&(d.k==k)].set_index('budget').loc[CONFIG['pair_budgets']]
            y=np.arange(4)+(offset-1)*.20
            ax.errorbar(g.difference,y,xerr=[g.difference-g.simultaneous_low,g.simultaneous_high-g.difference],
                        fmt='o',mfc='white',color=color,ms=4,capsize=2,label=f'k={k}')
            substantial=g.material_improvement.to_numpy(dtype=bool)
            ax.scatter(g.difference.to_numpy()[substantial],y[substantial],s=20,color=color,zorder=4)
        ax.axvline(0,color='black',lw=.8);ax.set_yticks(range(4),['512 B','2 KiB','8 KiB','32 KiB'])
        ax.set_title(METRIC_NAMES[metric],fontsize=11);ax.grid(axis='x',alpha=.2)
        ax.set_xlabel('Loss(rule) - loss(fixed); left is better',fontsize=9)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.ticklabel_format(axis='x',style='sci',scilimits=(-3,3),useMathText=True)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Frozen pilot rule vs development-fixed best | all 14 conditions | 30 held-out genera')
    fig.tight_layout(rect=(0,.025,1,.95))
    fig.text(.02,.005,'95% genus bootstrap max-t intervals over 72 contrasts. Filled points also meet the 10% and 0.0001 loss-reduction thresholds.\n'
        'Open points do not meet all criteria. Different objectives use different loss definitions.',fontsize=9)
    save(fig,'advisor_vs_fixed_simultaneous')


def budget_curves():
    d=pd.read_csv(OUT/'policy_summary.csv');d=d[(d.k==24)&(d.condition=='all')]
    fig,axes=plt.subplots(2,3,figsize=(14,8))
    for ax,metric in zip(axes.flat,METRICS):
        for policy in ['pilot','metadata','fixed']:
            g=d[(d.metric==metric)&(d.policy==policy)].set_index('budget').loc[CONFIG['pair_budgets']]
            ax.plot(range(4),g.loss,'o-',color=POLICY_COLORS[policy],label=POLICY_NAMES[policy])
        g=d[(d.metric==metric)&(d.policy=='pilot')].set_index('budget').loc[CONFIG['pair_budgets']]
        ax.plot(range(4),g.loss-g.mean_regret,'--',color='#888888',label='Offline penalized-loss oracle')
        ax.set_xticks(range(4),['512 B','2 KiB','8 KiB','32 KiB'])
        ax.set_title(METRIC_NAMES[metric],fontsize=11);ax.set_ylabel('Failure-aware loss; lower is better')
        ax.grid(alpha=.2)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=2,fontsize=9)
    fig.suptitle('Budget and selection error | k=24 | 14 validation conditions')
    fig.tight_layout(rect=(0,.09,1,.95));save(fig,'budget_selection_curves')


def target_attainment():
    d=pd.read_csv(OUT/'precision_profile_summary.csv')
    d=d[(d.k==24)&(d.budget==8192)&(d.profile=='standard')]
    fig,axes=plt.subplots(1,3,figsize=(16,9),sharey=True)
    for ax,policy in zip(axes,['fixed','metadata','pilot']):
        g=d[d.policy==policy].pivot(index='condition',columns='metric',values='attainment_rate').loc[ORDER,METRICS]
        im=ax.imshow(g.to_numpy(),vmin=0,vmax=1,cmap='YlGnBu',aspect='auto')
        for i in range(len(ORDER)):
            for j in range(6):
                v=g.iloc[i,j];ax.text(j,i,f'{v:.0%}',ha='center',va='center',fontsize=7,color='white' if v>.65 else 'black')
        ax.set_xticks(range(6),['J','C','Card','NN','Species','Cluster'],rotation=45,ha='right')
        ax.set_yticks(range(len(ORDER)),[LABELS[c] for c in ORDER],fontsize=8)
        ax.axhline(9.5,color='white',lw=1.5);ax.set_title(POLICY_NAMES[policy],fontsize=11)
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.6,label='Fraction of workloads meeting both error targets and availability')
    fig.suptitle('Actual target attainment, not predicted feasibility | k=24 | 8 KiB')
    fig.subplots_adjust(left=.26,right=.85,bottom=.13,top=.92,wspace=.12)
    save(fig,'precision_attainment_k24_8k')


def availability():
    d=pd.read_csv(OUT/'policy_summary.csv');d=d[(d.k==24)&(d.metric=='jaccard')&d.policy.isin(ALGORITHMS)]
    fig,axes=plt.subplots(1,4,figsize=(16,9),sharey=True)
    for ax,budget in zip(axes,CONFIG['pair_budgets']):
        g=d[d.budget==budget].pivot(index='condition',columns='policy',values='availability').loc[ORDER,ALGORITHMS]
        im=ax.imshow(g.to_numpy(),vmin=0,vmax=1,cmap='RdYlBu',aspect='auto')
        for i in range(len(ORDER)):
            for j in range(7):
                v=g.iloc[i,j]
                if v<.99:ax.text(j,i,f'{v:.0%}',ha='center',va='center',fontsize=6,color='black')
        ax.set_xticks(range(7),[SHORT[a] for a in ALGORITHMS],rotation=60,ha='right')
        ax.set_yticks(range(len(ORDER)),[LABELS[c] for c in ORDER],fontsize=8)
        ax.set_title(f'{budget/1024:g} KiB per pair')
        ax.axhline(9.5,color='black',lw=.8)
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.6,label='Valid and within budget')
    fig.suptitle('Empty/invalid estimates and budget failures remain visible | Jaccard | k=24')
    fig.subplots_adjust(left=.27,right=.85,bottom=.16,top=.92,wspace=.13)
    save(fig,'availability_and_budget_failures')


def downstream():
    d=pd.read_csv(OUT/'policy_summary.csv');d=d[(d.k==24)&(d.budget==8192)]
    exact=pd.read_csv(OUT/'exact_downstream_summary.csv');exact=exact[(exact.k==24)&(exact.budget==8192)].set_index('condition')
    fig,axes=plt.subplots(3,1,figsize=(14,11),sharex=True)
    for ax,metric in zip(axes,['neighbor','species','clustering']):
        for policy in ['pilot','metadata','fixed']:
            g=d[(d.metric==metric)&(d.policy==policy)].set_index('condition').loc[ORDER]
            y=g.loss if metric=='clustering' else 1-g.loss
            ax.plot(range(len(ORDER)),y,'o-',color=POLICY_COLORS[policy],label=POLICY_NAMES[policy],ms=4)
        y=(1-exact.loc[ORDER,'ari'])/2 if metric=='clustering' else exact.loc[ORDER,'species_accuracy'] if metric=='species' else np.ones(len(ORDER))
        ax.plot(range(len(ORDER)),y,'--',color='#999999',label='Exact k-mer reference')
        ax.set_ylabel('Clustering loss' if metric=='clustering' else 'Recovery' if metric=='neighbor' else 'Accuracy')
        ax.axvline(9.5,color='black',lw=.8);ax.grid(alpha=.2);ax.set_title(METRIC_NAMES[metric],fontsize=11)
    axes[0].legend(fontsize=9,ncol=2)
    axes[-1].set_xticks(range(len(ORDER)),[LABELS[c] for c in ORDER],rotation=45,ha='right',fontsize=8)
    fig.suptitle('Harder downstream tasks | k=24 | 8 KiB | six-reference known-genus galleries')
    fig.tight_layout(rect=(0,.01,1,.97));save(fig,'downstream_tasks_k24_8k')


def pilot_diagnostics():
    workloads=json.loads((DATA/'validation_workloads.json').read_text());rows=[]
    features=pd.read_csv(OUT/'validation_features.csv').set_index(['workload_id','k'])
    for k in CONFIG['k']:
        p=pd.read_csv(OUT/'validation'/f'pilot_k{k}_pairs.csv.gz').set_index('pair_id')
        o=pd.read_csv(OUT/'validation'/f'pilot_k{k}_objects.csv.gz').set_index('object_id')
        t=pd.read_csv(OUT/'validation'/f'exact_k{k}_pairs.csv.gz').set_index('pair_id')
        for w in workloads:
            q=p.loc[w['retrieval_pairs']];truth=t.loc[w['retrieval_pairs']]
            c=bounded_containment(q.jaccard,o.loc[q.a,'cardinality'].to_numpy(),o.loc[q.b,'cardinality'].to_numpy())
            f=features.loc[(w['workload_id'],k),:]
            rows.append(dict(workload_id=w['workload_id'],genus=w['genus'],condition=w['condition'],k=k,
                pilot_ratio_p90=10**f.log_pilot_ratio_p90,pilot_j_mae=float(abs(q.jaccard-truth.truth_j).mean()),
                pilot_c_mae=float(np.mean(abs(c-truth.truth_c.to_numpy()))),
                pilot_zero_fraction=float(np.mean(q.jaccard==0)),truth_j_mean=float(truth.truth_j.mean()),
                truth_c_mean=float(truth.truth_c.mean()),pilot_payload_bytes=int(f.pilot_payload_bytes)))
    frame=pd.DataFrame(rows);save_csv(OUT/'pilot_diagnostics.csv',frame)
    g=frame[frame.k==24].groupby('condition').mean(numeric_only=True).loc[ORDER]
    fig,axes=plt.subplots(1,3,figsize=(15,8),sharey=True)
    axes[0].barh(range(len(ORDER)),g.pilot_zero_fraction,color='#0072B2')
    axes[1].barh(range(len(ORDER)),g.pilot_c_mae,color='#D55E00')
    axes[2].barh(range(len(ORDER)),g.truth_c_mean,color='#009E73')
    for ax,title in zip(axes,['Pilot Jaccard = 0 fraction','Pilot containment MAE','Exact query containment (mean)']):
        ax.set_title(title,fontsize=10);ax.set_xlim(0,1);ax.grid(axis='x',alpha=.2);ax.axhline(9.5,color='black',lw=.8)
    axes[0].set_yticks(range(len(ORDER)),[LABELS[c] for c in ORDER],fontsize=8);axes[0].invert_yaxis()
    fig.suptitle('Pilot resolution limits | K=32 | independent seed 919 | k=24')
    fig.tight_layout(rect=(0,.01,1,.95));save(fig,'pilot_resolution_limits')


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    candidates();primary_intervals();budget_curves();target_attainment();availability();downstream();pilot_diagnostics()
    atomic(OUT/'figures_complete.json',dict(figures=7,png=7,pdf=7,code_sha256=sha(__file__)))


if __name__=='__main__':main()
