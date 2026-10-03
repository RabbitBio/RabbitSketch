#!/usr/bin/env python3
"""Deterministic real-sequence mutations and fragments; no benchmark outcomes used."""
from pathlib import Path
import csv, hashlib, json, re
import numpy as np
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark')
OUT=ROOT/'data/regimes_20261001';GEN=OUT/'controlled';GEN.mkdir(exist_ok=True)
def tsv(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def seqs(path):
    return [re.sub(r'\s+','',s.split('\n',1)[1]).upper() for s in path.read_text().split('>')[1:]]
def save(name,seq):
    path=GEN/(name+'.fna')
    if not path.exists():
        with path.open('w') as f:f.write('>'+name+'\n');f.write('\n'.join(seq[i:i+100] for i in range(0,len(seq),100))+'\n')
    return str(path)
def object_row(oid,seq):return dict(object_id=oid,path=save(oid,seq))
old=list(csv.DictReader((ROOT/'data/real_genome_200_20261001/pairs.tsv').open(),delimiter='\t'))
candidates={r['file'+s]:r['org'+s].split()[0] for r in old for s in ('A','B')}
refs=[];seen=set()
for path in sorted(candidates,key=lambda p:hashlib.sha256(Path(p).name.encode()).hexdigest()):
    genus=candidates[path]
    if genus in seen:continue
    p=Path(path)
    if p.stat().st_size<5_050_000:continue
    sequences=seqs(p);seq=max(sequences,key=len)
    if len(seq)<5_000_000:continue
    seq=seq[:5_000_000]
    if sum(seq.count(c) for c in 'ACGT')<4_950_000:continue
    refs.append((path,genus,seq));seen.add(genus)
    if len(refs)==20:break
if len(refs)<10:raise RuntimeError('not enough independent genera with 5 Mb contigs')
print('references',len(refs),flush=True)
objects=[];pairs=[];provenance=[];contain_objects=[];contain_pairs=[];card_objects=[];card_meta=[]
target=[.80,.85,.90,.94,.945,.95,.955,.96,.97,.985,.99,.995,.999,1.]
for i,(path,genus,seq) in enumerate(refs):
    rid=f'ref{i:02d}';base=object_row(rid,seq);objects.append(base);contain_objects.append(base)
    split='calibration' if i<len(refs)//2 else 'validation'
    provenance.append(dict(reference=rid,path=path,genus=genus,split=split,source_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),bases=len(seq)))
    a=np.frombuffer(seq.encode(),dtype=np.uint8);valid=np.flatnonzero(np.isin(a,[65,67,71,84]));lookup=np.zeros(256,dtype=np.uint8);lookup[[65,67,71,84]]=[0,1,2,3];letters=np.array([65,67,71,84],dtype=np.uint8)
    rng=np.random.default_rng(2026100100+i);perm=rng.permutation(valid);delta=rng.integers(1,4,size=len(valid),dtype=np.uint8)
    for j,identity in enumerate(target):
        oid=f'{rid}_ani{identity:.4f}';n=int(round((1-identity)*len(valid)));b=a.copy();pos=perm[:n];b[pos]=letters[(lookup[b[pos]]+delta[:n])%4]
        objects.append(object_row(oid,b.tobytes().decode()))
        pairs.append(dict(pair_id=oid,a=rid,b=oid,family='substitution',reference=rid,split=split,target_identity=identity,identity=1-n/len(valid),valid_positions=len(valid),mutations=n))
    donor=refs[(i+len(refs)//2)%len(refs)][2]
    for length in [5_000_000,2_000_000,1_000_000,500_000,100_000,50_000,10_000]:
        for overlap in [1.,.8,.5,.2,.1]:
            oid=f'{rid}_frag{length}_c{overlap:.1f}';shared=int(length*overlap)
            # Two records prevent an artificial junction k-mer.
            query=seq[:shared]+'N'+donor[:length-shared]
            contain_objects.append(object_row(oid,query));contain_pairs.append(dict(pair_id=oid,a=oid,b=rid,family='containment',reference=rid,split=split,query_length=length,reference_length=len(seq),nominal_ratio=len(seq)/length,nominal_containment=overlap,donor=refs[(i+len(refs)//2)%len(refs)][1]))
    for length in [10_000,100_000,1_000_000,5_000_000]:
        oid=f'{rid}_card{length}';card_objects.append(object_row(oid,seq[:length]));card_meta.append(dict(object_id=oid,length=length,reference=rid,kind='bacterial',split=split))
    print(rid,genus,'prepared',flush=True)
# Five nonoverlapping human windows; read only the first two chromosomes.
human=Path('/home/bigssd/zt/human.fna');hseq=[];chrom='';hindex=0
def human_windows(chrom,sequence):
    global hindex
    for offset in range(0,len(sequence)-50_000_000+1,50_000_000):
        if hindex==5:return
        oid=f'human_card50m_{hindex}';card_objects.append(object_row(oid,sequence[offset:offset+50_000_000]));card_meta.append(dict(object_id=oid,length=50_000_000,reference=chrom+':'+str(offset),kind='human',split='validation'));hindex+=1
with human.open() as f:
    for line in f:
        if line.startswith('>'):
            if hseq:human_windows(chrom,''.join(hseq));hseq=[]
            if hindex==5:break
            chrom=line[1:].split()[0]
        else:hseq.append(line.strip().upper())
if hindex<5:human_windows(chrom,''.join(hseq))
for name,rows in [('ani_control_objects',objects),('ani_control_pairs',pairs),('contain_objects',contain_objects),('contain_pairs',contain_pairs),('card_objects',card_objects),('card_meta',card_meta),('controlled_provenance',provenance)]:tsv(OUT/(name+'.tsv'),rows)
(OUT/'controlled_preparation.json').write_text(json.dumps(dict(reference_count=len(refs),ani_pairs=len(pairs),containment_pairs=len(contain_pairs),cardinality_objects=len(card_objects),data_seed=2026100100),indent=2)+'\n')
print((OUT/'controlled_preparation.json').read_text(),flush=True)
