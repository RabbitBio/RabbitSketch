#!/usr/bin/env python3
"""Correct donor leakage: both constituent genomes stay inside one split."""
from pathlib import Path
import csv, json, re
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark');OLD=ROOT/'data/regimes_20261001';OUT=ROOT/'data/containment_holdout_20261001';OUT.mkdir(exist_ok=True);GEN=OUT/'genomes';GEN.mkdir(exist_ok=True)
prov=list(csv.DictReader((OLD/'controlled_provenance.tsv').open(),delimiter='\t'))
def load(p):return ''.join(s.strip() for s in p.read_text().splitlines() if not s.startswith('>'))
def save(name,seq):
    p=GEN/(name+'.fna')
    if not p.exists():p.write_text('>'+name+'\n'+'\n'.join(seq[i:i+100] for i in range(0,len(seq),100))+'\n')
    return str(p)
refs={r['reference']:load(OLD/'controlled'/(r['reference']+'.fna')) for r in prov}
groups={s:[r['reference'] for r in prov if r['split']==s] for s in ['calibration','validation']}
objects=[];pairs=[]
for r in prov:
    rid=r['reference'];seq=refs[rid];group=groups[r['split']];donor_id=group[(group.index(rid)+1)%len(group)];donor=refs[donor_id]
    objects.append(dict(object_id=rid,path=str(OLD/'controlled'/(rid+'.fna'))))
    for length in [5_000_000,2_000_000,1_000_000,500_000,100_000,50_000,10_000]:
        for overlap in [1.,.8,.5,.2,.1]:
            name=f'{rid}_frag{length}_c{overlap:.1f}';shared=int(length*overlap);query=seq[:shared]+'N'+donor[:length-shared]
            objects.append(dict(object_id=name,path=save(name,query)))
            pairs.append(dict(pair_id=name,a=name,b=rid,family='containment',reference=rid,split=r['split'],query_length=length,reference_length=len(seq),nominal_ratio=len(seq)/length,nominal_containment=overlap,donor=donor_id))
for name,rows in [('contain_objects',objects),('contain_pairs',pairs)]:
    with (OUT/(name+'.tsv')).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
for s in ['calibration','validation']:
    used={p[k] for p in pairs if p['split']==s for k in ['reference','donor']};assert used==set(groups[s])
assert set(groups['calibration']).isdisjoint(groups['validation'])
(OUT/'provenance.json').write_text(json.dumps(dict(parent=str(OLD),change='replace cross-split donor by next reference within the same split; all 700 pairs rebuilt and exact truth recomputed',split_genomes=groups),indent=2)+'\n')
print('Prepared',len(pairs),'containment pairs with disjoint source-genome splits',flush=True)
