#!/usr/bin/env python3
"""Audit the public artifacts, frozen choices, alternative scores and provenance."""
from advisor import Advisor, PILOT_FEATURES
from scenario_common import *


def main():
    manifest=json.loads((PUBLISHED/'artifact_sha256.json').read_text())
    for name,digest in manifest.items():
        if sha(PUBLISHED/name)!=digest:raise ValueError('published artifact changed: '+name)
    files={str(p.relative_to(PUBLISHED)) for p in PUBLISHED.rglob('*') if p.is_file()
           and p.name not in ['artifact_sha256.json','publication_verification.json']}
    if files!=set(manifest):raise ValueError('unmanifested publication files')
    for name,digest in json.loads((PUBLISHED/'workflow_sha256.json').read_text()).items():
        if sha(HERE/name)!=digest:raise ValueError('published workflow source changed: '+name)
    dev=pd.read_csv(PUBLISHED/'development_genomes.tsv',sep='\t')
    val=pd.read_csv(PUBLISHED/'validation_genomes.tsv',sep='\t')
    for column in ['object_id','biosample','fingerprint','genus']:
        if set(dev[column])&set(val[column]):raise ValueError('development/validation overlap: '+column)
    if val.fingerprint.nunique()!=len(val) or val.biosample.nunique()!=len(val):raise ValueError('validation duplicates')
    if not (val.groupby(['genus','species']).size()==5).all():raise ValueError('unbalanced species')
    if not (val.groupby('genus').species.nunique()==3).all():raise ValueError('unbalanced genera')
    decisions=pd.read_csv(PUBLISHED/'validation_decisions.csv.gz');decisions['config']=decisions.config.fillna('')
    outcomes=pd.read_csv(PUBLISHED/'decision_outcomes.csv.gz');outcomes['config']=outcomes.config.fillna('')
    native=pd.read_csv(PUBLISHED/'validation_native_scores.csv.gz')
    decision_keys=['workload_id','k','budget','metric','policy']
    native_keys=['workload_id','k','budget','metric','config']
    if decisions.duplicated(decision_keys).any() or outcomes.duplicated(decision_keys).any():raise ValueError('duplicate decision')
    if native.duplicated(native_keys).any():raise ValueError('duplicate alternative')
    joined=decisions[decision_keys+['config']].merge(outcomes[decision_keys+['config']],on=decision_keys,validate='one_to_one',suffixes=('_frozen','_scored'))
    if len(joined)!=len(decisions) or not (joined.config_frozen==joined.config_scored).all():raise ValueError('choice changed after freeze')
    selected=outcomes[outcomes.config.ne('')].merge(native[native_keys+['loss','availability','budget_violation']],on=native_keys,validate='many_to_one',suffixes=('_scored','_native'))
    if len(selected)!=outcomes.config.ne('').sum():raise ValueError('selected native score missing')
    for column in ['loss','availability','budget_violation']:
        if not np.allclose(selected[column+'_scored'],selected[column+'_native'],rtol=0,atol=1e-12,equal_nan=True):
            raise ValueError('selected score changed: '+column)
    unsupported=outcomes[outcomes.config.eq('')]
    if not ((unsupported.loss==1)&(unsupported.availability==0)).all():raise ValueError('unsupported output hidden')
    real=native[native.algorithm!='exact']
    oracle=real.groupby(native_keys[:-1]).loss.min().rename('fresh_oracle').reset_index()
    check=outcomes.merge(oracle,on=native_keys[:-1],validate='many_to_one')
    if not np.allclose(check.oracle_loss,check.fresh_oracle,rtol=0,atol=1e-12):raise ValueError('offline oracle mismatch')
    if ((outcomes.loss<0)|(outcomes.loss>1+1e-10)|~np.isfinite(outcomes.loss)).any():raise ValueError('invalid failure-aware loss')
    features=pd.read_csv(PUBLISHED/'validation_features.csv')
    model=Advisor(PUBLISHED/'model')
    if model.metadata['scalings']['pilot']['fields']!=PILOT_FEATURES:raise ValueError('unexpected model features')
    forbidden=['truth','species','genus','nominal','coverage','contamination','error_rate']
    if any(any(word in f for word in forbidden) for f in PILOT_FEATURES):raise ValueError('forbidden model feature')
    # Re-evaluate one deterministic workload from every input condition and k,
    # across every task/budget. No native measurements or exact truth are used.
    chosen=features.sort_values('workload_id').groupby(['condition','k'],sort=True).first().reset_index()
    lookup=decisions.set_index(decision_keys)
    recomputed=0
    for f in chosen.to_dict('records'):
        for policy in ['pilot','metadata']:
            neighbors=model.neighbors(f,int(f['k']),f['kind'],policy)
            for budget in CONFIG['pair_budgets']:
                for metric in model.metadata['metrics']:
                    result=model.choose(f,int(f['k']),f['kind'],metric,budget,policy,neighbors)
                    original=lookup.loc[(f['workload_id'],f['k'],budget,metric,policy),:]
                    if result['config']!=original.config:raise ValueError('portable frozen decision mismatch')
                    if abs(result['predicted_loss']-original.predicted_loss)>1e-7:raise ValueError('portable prediction mismatch')
                    recomputed+=1
    inventory=pd.read_csv(PUBLISHED/'native_archive_manifest.csv')
    if inventory.job_id.duplicated().any() or inventory.binary_sha256.nunique()!=1:raise ValueError('native identity inconsistency')
    prep={split:json.loads((PUBLISHED/(split+'_preparation.json')).read_text()) for split in ['development','validation']}
    for split,group in inventory.groupby('split'):
        if not (group.pairs==prep[split]['pairs']).all() or not (group.objects==prep[split]['objects']).all():
            raise ValueError('incomplete native input panel')
    evaluation=json.loads((PUBLISHED/'evaluation_complete.json').read_text())
    if not evaluation['validation_outcomes_after_decisions']:raise ValueError('outcomes preceded decisions')
    inputs=json.loads((PUBLISHED/'inputs_verified.json').read_text())
    if not inputs['all_match'] or inputs['manifest_sha256']!=sha(PUBLISHED/'input_sha256.json'):
        raise ValueError('input verification mismatch')
    for name in ['model_frozen.json','decisions_frozen.json']:
        identity=json.loads((RUN/name).read_text())
        for filename,digest in identity['files'].items():
            if sha(filename)!=digest:raise ValueError('original frozen artifact changed')
    result=dict(artifacts=len(manifest),all_hashes_match=True,validation_genomes=len(val),
        validation_genera=int(val.genus.nunique()),frozen_decisions=len(decisions),
        selected_scores_checked=len(selected),unsupported_decisions_preserved=len(unsupported),
        portable_decisions_recomputed=recomputed,native_archives=len(inventory),
        full_input_files_verified=inputs['files'],code_sha256=sha(__file__))
    atomic(PUBLISHED/'publication_verification.json',result)
    atomic(RUN/'complete.json',dict(status='complete',publication=str(PUBLISHED),
        artifact_manifest_sha256=sha(PUBLISHED/'artifact_sha256.json'),
        verification_sha256=sha(PUBLISHED/'publication_verification.json'),**result))
    atomic(RUN/'pipeline_status.json',dict(stage='complete',publication=str(PUBLISHED)))
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
