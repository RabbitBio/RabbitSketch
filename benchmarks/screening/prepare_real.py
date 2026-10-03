#!/usr/bin/env python3
"""Resolve the historical pair manifest to identical-accession local gzip sources."""
import csv
import gzip
import json
from pathlib import Path
import shutil

from run import DEFAULT_ROOT,HERE,atomic_json,file_sha

original=HERE.parent.parent/'examples/real_genome_eval/pairs.tsv'
out=DEFAULT_ROOT/'data/real_genome_200_20261001'
out.mkdir(parents=True,exist_ok=True)
genomes=out/'genomes';genomes.mkdir(exist_ok=True)
rows=list(csv.DictReader(original.open(),delimiter='\t'))
oldfiles=sorted({r[k] for r in rows for k in ('fileA','fileB')})
resolved={};sources=[]
for i,name in enumerate(oldfiles):
    source=Path('/home/bigssd/zt/ref/bacteria')/(Path(name).name+'.gz')
    target=genomes/Path(name).name
    if not target.exists():
        tmp=target.with_suffix('.fna.tmp')
        with gzip.open(source,'rb') as src, tmp.open('wb') as dest: shutil.copyfileobj(src,dest)
        tmp.replace(target)
    resolved[name]=str(target)
    sources.append(dict(historical_path=name,compressed_source=str(source),
                        compressed_sha256=file_sha(source),path=str(target),sha256=file_sha(target)))
    if (i+1)%50==0:print('Resolved',i+1,'/',len(oldfiles),flush=True)
for row in rows:
    for key in ('fileA','fileB'):row[key]=resolved[row[key]]
with (out/'pairs.tsv').open('w') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
atomic_json(out/'manifest.json',dict(original_pairs_path=str(original),original_pairs_sha256=file_sha(original),
    resolved_pairs_sha256=file_sha(out/'pairs.tsv'),sources=sources,
    note='Same named accessions and local source filenames; historic uncompressed file hashes unavailable.'))
print('Prepared',len(rows),'pairs,',len(oldfiles),'genomes',flush=True)
