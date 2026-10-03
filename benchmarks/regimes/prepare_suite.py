#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, json
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark');OUT=ROOT/'data/regimes_20261001'
def read(name):return list(csv.DictReader((OUT/(name+'.tsv')).open(),delimiter='\t'))
def write(name,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/(name+'.tsv')).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
objects=read('mash_objects')+read('ani_control_objects');pairs=read('mash_pairs')+read('ani_control_pairs')
mashsplit={o['object_id']:('calibration' if i<71 else 'validation') for i,o in enumerate(sorted(read('mash_objects'),key=lambda r:hashlib.sha256(r['object_id'].encode()).hexdigest()))}
for p in pairs:
    if p['family']=='mash142':p['split']=mashsplit[p['a']] if mashsplit[p['a']]==mashsplit[p['b']] else 'cross_split'
old=list(csv.DictReader((ROOT/'data/real_genome_200_20261001/pairs.tsv').open(),delimiter='\t'))
found={}
for p in old:
    for s in ['A','B']:
        oid='real_'+p['acc'+s]
        if oid not in found:objects.append(dict(object_id=oid,path=p['file'+s]));found[oid]=1
    pairs.append(dict(pair_id='real_'+p['pair_id'],a='real_'+p['accA'],b='real_'+p['accB'],family='diverse200',reference=p['accA'],split='validation',tier=p['tier']))
write('ani_objects',objects);write('ani_pairs',pairs)
allpaths=sorted(Path('/home/bigssd/zt/ref/bacteria').glob('GC[AF]*.fna.gz'),key=lambda p:hashlib.sha256(p.name.encode()).hexdigest())
workload=[]
for p in allpaths:
    if p.stat().st_size>10_000_000:continue
    workload.append(dict(object_id=p.name.split('_genomic')[0],path=str(p)))
    if len(workload)==10_000:break
write('workload_objects',workload)
print('ANI objects',len(objects),'pairs',len(pairs),'workload',len(workload),flush=True)
