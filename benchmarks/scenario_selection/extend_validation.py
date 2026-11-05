#!/usr/bin/env python3
"""Register unseen conditions before any validation pilot/outcome, without refitting."""
from concurrent.futures import ThreadPoolExecutor
import datetime
import shutil

from prepare import derive
from scenario_common import *


def main():
    done = DATA / 'validation_extension_complete.json'
    if done.exists(): print('Validation extension already prepared'); return
    if (RUN / 'validation-pilot_plan.json').exists() or (RUN / 'model_frozen.json').exists():
        raise ValueError('register held-out conditions before model freeze and validation')
    conditions_path = HERE / 'validation_extra_conditions.json'
    conditions = json.loads(conditions_path.read_text())
    identity = dict(registered_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        code_sha256=sha(__file__), extra_conditions_sha256=sha(conditions_path),
        addendum_sha256=sha(HERE / 'ANALYSIS_ADDENDUM.md'),
        original_preparation_sha256=sha(DATA / 'validation_preparation.json'),
        no_validation_pilot_or_outcomes=True)
    marker = DATA / 'validation_extension_registration.json'
    if not marker.exists(): atomic(marker, identity)
    else:
        prior = json.loads(marker.read_text())
        for k in ['code_sha256','extra_conditions_sha256','addendum_sha256']:
            if prior[k] != identity[k]: raise ValueError('registered extension changed')
    objects = read_tsv(DATA / 'validation_objects.tsv')
    genomes = {r['object_id']: r for r in read_tsv(DATA / 'validation_genomes.tsv')}
    sources = [r for r in objects if r['condition'] == 'full' and r['role'] == 'query']
    donor_for = {r['source']: r['donor'] for r in objects if r['condition'] != 'full'}
    old_workloads = json.loads((DATA / 'validation_workloads.json').read_text())
    references = {w['genus']: w['references'] for w in old_workloads}
    def create(source):
        sid = source['source']; donor_id = donor_for[sid]
        records = fasta_records(source['path']); donor = fasta_records(genomes[donor_id]['path'])
        result = []
        for condition in conditions:
            oid = sid + '__' + condition['name']; path = DATA / 'queries/validation' / (oid + '.fa')
            meta = path.with_suffix('.json')
            if meta.exists():
                row = json.loads(meta.read_text())
                if row['sha256'] != sha(path): raise ValueError('extended query changed')
            else:
                seqs, provenance = derive(records, donor, condition, sid)
                with path.open('wb') as f:
                    for i, s in enumerate(seqs): f.write(b'>' + str(i).encode() + b'\n' + s + b'\n')
                lengths = list(map(len, seqs))
                row = dict(object_id=oid,path=str(path),source=sid,donor=donor_id,split='validation',
                    genus=source['genus'],species=source['species'],kind=condition['kind'],condition=condition['name'],
                    bases=sum(lengths),records=len(lengths),median_record_length=float(np.median(lengths)),
                    sha256=sha(path),source_sha256=sha(source['path']),donor_sha256=sha(genomes[donor_id]['path']),
                    provenance=provenance,role='query')
                atomic(meta,row)
            result.append({k:v for k,v in row.items() if k!='provenance'})
        return result
    new_objects = []
    with ThreadPoolExecutor(max_workers=16) as pool:
        for i, rows in enumerate(pool.map(create,sources),1):
            new_objects.extend(rows)
            if i%30==0:print('Unseen-condition sources',i,'/',len(sources),flush=True)
    archive = DATA / 'before_validation_extension'; archive.mkdir(exist_ok=True)
    for suffix in ['objects.tsv','pairs.tsv','workloads.json','preparation.json']:
        p=DATA/('validation_'+suffix)
        if not (archive/p.name).exists():shutil.copyfile(p,archive/p.name)
    objects += new_objects
    if len({r['object_id'] for r in objects})!=len(objects):raise ValueError('duplicate extended object')
    pairs=read_tsv(DATA/'validation_pairs.tsv'); workloads=list(old_workloads)
    for genus,refs in sorted(references.items()):
        for condition in conditions:
            queries=sorted(r['object_id'] for r in new_objects if r['genus']==genus and r['condition']==condition['name'])
            if len(queries)!=9:raise ValueError('unbalanced held-out workload')
            retrieval=[];cluster=[]
            for q in queries:
                for r in refs:
                    pid=key((q,r))[:24];pairs.append(dict(pair_id=pid,a=q,b=r));retrieval.append(pid)
            for i,a in enumerate(queries):
                for b in queries[i+1:]:
                    pid=key((a,b))[:24];pairs.append(dict(pair_id=pid,a=a,b=b));cluster.append(pid)
            workloads.append(dict(workload_id=f'validation_{genus}_{condition["name"]}',split='validation',
                genus=genus,condition=condition['name'],kind=condition['kind'],queries=queries,references=refs,
                retrieval_pairs=retrieval,cluster_pairs=cluster))
    if len({p['pair_id'] for p in pairs})!=len(pairs):raise ValueError('duplicate extended pair')
    tsv(DATA/'validation_objects.tsv',objects);tsv(DATA/'validation_pairs.tsv',pairs)
    atomic(DATA/'validation_workloads.json',workloads)
    prep=json.loads((DATA/'validation_preparation.json').read_text())
    prep.update(objects=len(objects),pairs=len(pairs),workloads=len(workloads),
        derived_query_files=len(objects)-prep['genomes'],original_conditions=10,unseen_conditions=4,
        extension_registration_sha256=sha(marker),
        files={str(DATA/('validation'+s)):sha(DATA/('validation'+s)) for s in ['_objects.tsv','_pairs.tsv','_workloads.json']})
    atomic(DATA/'validation_preparation.json',prep)
    inputs=json.loads((RUN/'input_sha256.json').read_text());inputs.update({r['path']:r['sha256'] for r in new_objects})
    atomic(RUN/'input_sha256.json',inputs)
    atomic(done,dict(new_conditions=len(conditions),new_objects=len(new_objects),workloads=len(workloads),
        preparation_sha256=sha(DATA/'validation_preparation.json'),input_manifest_sha256=sha(RUN/'input_sha256.json'),
        registration_sha256=sha(marker)))
    overall=json.loads((DATA/'complete.json').read_text());overall['validation_extension_sha256']=sha(done)
    overall['input_manifest_sha256']=sha(RUN/'input_sha256.json');atomic(DATA/'complete.json',overall)
    print(json.dumps(prep,indent=2),flush=True)


if __name__=='__main__':main()
