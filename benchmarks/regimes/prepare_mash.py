#!/usr/bin/env python3
"""Resolve exact local accession versions and recompute FastANI references."""
from pathlib import Path
import csv, gzip, hashlib, json, re, shutil
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark')
OUT=ROOT/'data/regimes_20261001'
OUT.mkdir(parents=True,exist_ok=True)
(OUT/'mash_genomes').mkdir(exist_ok=True)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def write_tsv(p,rows):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
paths={re.match(r'GC[AF]_\d+\.\d+',p.name).group():p for p in Path('/home/bigssd/zt/ref/bacteria').glob('GC[AF]*.fna.gz')}
source=list(csv.DictReader((ROOT/'sources/mash500/genome_manifest.tsv').open(),delimiter='\t'))
objects=[]; provenance=[]; missing=[]
for r in source:
    a=r['assembly_accession']
    if a not in paths:missing.append(r);continue
    gz=paths[a]; dest=OUT/'mash_genomes'/gz.name.removesuffix('.gz')
    if not dest.exists():
        tmp=dest.with_suffix('.tmp')
        with gzip.open(gz,'rb') as f,tmp.open('wb') as o:shutil.copyfileobj(f,o)
        tmp.replace(dest)
    oid='mash_'+r['genome_id']
    objects.append(dict(object_id=oid,path=str(dest)))
    provenance.append(dict(**r,path=str(dest),sha256=sha(dest),source_gz=str(gz),source_sha256=sha(gz)))
write_tsv(OUT/'mash_objects.tsv',objects)
write_tsv(OUT/'mash_unresolved.tsv',missing)
write_tsv(OUT/'mash_provenance.tsv',provenance)
pairs=[]
for i,a in enumerate(objects):
    for b in objects[i+1:]:pairs.append(dict(pair_id=a['object_id']+'__'+b['object_id'],a=a['object_id'],b=b['object_id'],family='mash142'))
write_tsv(OUT/'mash_pairs.tsv',pairs)
real=list(csv.DictReader((ROOT/'data/real_genome_200_20261001/pairs.tsv').open(),delimiter='\t'))
all_paths=sorted({o['path'] for o in objects}|{p[c] for p in real for c in ('fileA','fileB')})
(OUT/'fastani_genomes.txt').write_text(''.join(p+'\n' for p in all_paths))
(OUT/'mash_resolution.json').write_text(json.dumps(dict(original_genomes=500,matched_exact_versions=len(objects),pairs=len(pairs),fastani_genomes=len(all_paths),historical_ani_reused=False),indent=2)+'\n')
print((OUT/'mash_resolution.json').read_text(),flush=True)
