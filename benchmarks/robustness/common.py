from pathlib import Path
import csv,hashlib,json,os,re
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark')
DATA=ROOT/'data/robustness_20261001'
RUN=ROOT/'runs/robustness_20261001'
ANALYSIS=ROOT/'analysis/robustness_20261001'
HERE=Path(__file__).resolve().parent
for p in [DATA,RUN,ANALYSIS]:p.mkdir(parents=True,exist_ok=True)
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def key(value):return hashlib.sha256(('rs2-robustness-20261001:'+str(value)).encode()).hexdigest()
def read_tsv(path):return list(csv.DictReader(open(path),delimiter='\t'))
def tsv(path,rows):
    if not rows:raise ValueError('empty table '+str(path))
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with open(path,'w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def atomic(path,obj):
    p=Path(path);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n');tmp.replace(p)
def accroot(text):
    m=re.search(r'GC[AF]_(\d+)',str(text));return m[1] if m else None
def fasta_records(raw):
    records=[]
    for block in raw.split(b'>')[1:]:
        head,sep,seq=block.partition(b'\n')
        if not sep:continue
        seq=re.sub(rb'\s+',b'',seq).upper();seq=re.sub(rb'[^ACGT]',b'N',seq)
        records.append((head.split()[0],seq))
    if not records:raise ValueError('empty FASTA')
    return records
RC=bytes.maketrans(b'ACGTN',b'TGCAN')
def fingerprint(records):
    parts=sorted(hashlib.sha256(min(seq,seq.translate(RC)[::-1])).digest() for _,seq in records)
    return hashlib.sha256(b''.join(parts)).hexdigest()
