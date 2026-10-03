#!/usr/bin/env python3
"""Independent NUCmer/dnadiff check on a predeclared stratified natural subset."""
from pathlib import Path
import csv, fcntl, hashlib, json, os, re, shutil, subprocess, tarfile, time, urllib.request
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark');DATA=ROOT/'data/regimes_20261001';OUT=ROOT/'runs/regimes_20261001/alignment_audit';OUT.mkdir(exist_ok=True)
TOOLS=ROOT/'tools';TOOLS.mkdir(exist_ok=True);SRC=TOOLS/'mummer-4.0.1';INSTALL=TOOLS/'mummer-4.0.1-install'
URL='https://github.com/mummer4/mummer/releases/download/v4.0.1/mummer-4.0.1.tar.gz';archive=TOOLS/'mummer-4.0.1.tar.gz'
if not archive.exists():
    tmp=archive.with_suffix('.tmp');urllib.request.urlretrieve(URL,tmp);tmp.replace(archive)
if not SRC.exists():
    with tarfile.open(archive) as t:
        for m in t.getmembers():
            if not (m.name=='mummer-4.0.1' or m.name.startswith('mummer-4.0.1/')) or '..' in Path(m.name).parts:raise ValueError('unsafe archive path')
        t.extractall(TOOLS,filter='data')
print('Waiting for single benchmark slot',flush=True)
lock=(ROOT/'runs/.active_benchmark.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
def command(cmd,log,cwd=None,timeout=None):
    with log.open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,cwd=cwd,check=True,timeout=timeout)
if not (INSTALL/'bin/nucmer').exists():
    command(['./configure','--prefix='+str(INSTALL)],OUT/'configure.log',SRC)
    command(['make','-j64'],OUT/'make.log',SRC)
    command(['make','install'],OUT/'install.log',SRC)
objects={r['object_id']:r['path'] for r in csv.DictReader((DATA/'ani_objects.tsv').open(),delimiter='\t')}
pairs=list(csv.DictReader((DATA/'ani_pairs.tsv').open(),delimiter='\t'));obs={}
for line in (ROOT/'runs/regimes_20261001/fastani/references.tsv').open():
    a,b,x,m,n=line.rstrip().split('\t');obs[a,b]=(float(x)/100,int(m)/int(n))
natural=[]
for p in pairs:
    if p['family']=='substitution':continue
    a,b=objects[p['a']],objects[p['b']];x=obs.get((a,b));y=obs.get((b,a))
    if not x or not y or min(x[1],y[1])<.2:continue
    ident=(x[0]+y[0])/2
    if ident<.8:continue
    natural.append(dict(**p,file_a=a,file_b=b,fastani=ident))
# Eight per family/bin, plus eight near each threshold; deterministic hashes.
selected={}
for family in ['mash142','diverse200']:
    for lo,hi in [(.8,.9),(.9,.95),(.95,.99),(.99,1.000001),(.94,.96),(.985,.995)]:
        pool=[p for p in natural if p['family']==family and lo<=p['fastani']<hi]
        pool.sort(key=lambda p:hashlib.sha256(('nucmer-audit-v1'+p['pair_id']).encode()).hexdigest())
        for p in pool[:8]:selected[p['pair_id']]=p
plan=list(selected.values());fields=list(dict.fromkeys(k for p in plan for k in p))
with (OUT/'pairs.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(plan)
manifest=dict(source=URL,archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),version='4.0.1',threads=64,pairs=len(plan),selection='eight hash-selected pairs per family / ANI bin plus threshold windows; union; no sketch error used',nucmer_options=['--maxmatch','-c','65','-l','20','-t','64'])
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
results=[]
for i,p in enumerate(plan):
    d=OUT/p['pair_id'];d.mkdir(exist_ok=True);prefix=d/'align';completed=d/'result.json'
    if completed.exists():results.append(json.loads(completed.read_text()));continue
    print(f'{i+1}/{len(plan)} {p["pair_id"]}',flush=True);start=time.monotonic()
    cmd=[str(INSTALL/'bin/nucmer'),'--maxmatch','-c','65','-l','20','-t','64','-p',str(prefix),p['file_a'],p['file_b']]
    command(cmd,d/'nucmer.log',timeout=600)
    command([str(INSTALL/'bin/dnadiff'),'-d',str(prefix)+'.delta','-p',str(prefix)],d/'dnadiff.log',timeout=600)
    report=Path(str(prefix)+'.report').read_text();block=re.split(r'\n1-to-1\s',report)[1].split('\nM-to-M')[0]
    match=re.search(r'AvgIdentity\s+([\d.]+)\s+([\d.]+)',block)
    if not match:raise ValueError('missing 1-to-1 identity')
    af=re.search(r'AlignedBases\s+\d+\(([\d.]+)%\)\s+\d+\(([\d.]+)%\)',report)
    if not af:raise ValueError('missing aligned fraction')
    lengths=re.search(r'TotalLength\s+(\d+)\s+(\d+)',block);totals=re.search(r'TotalBases\s+(\d+)\s+(\d+)',report)
    result=dict(pair_id=p['pair_id'],family=p['family'],fastani=p['fastani'],nucmer_ref_identity=float(match[1])/100,nucmer_query_identity=float(match[2])/100,af_ref=float(af[1])/100,af_query=float(af[2])/100,one_to_one_length_fraction_ref=int(lengths[1])/int(totals[1]),one_to_one_length_fraction_query=int(lengths[2])/int(totals[2]),wall_seconds=time.monotonic()-start,command=cmd,report_sha256=hashlib.sha256(report.encode()).hexdigest())
    completed.write_text(json.dumps(result,indent=2)+'\n');results.append(result)
with (OUT/'results.json').open('w') as f:json.dump(results,f,indent=2)
(OUT/'completion.json').write_text(json.dumps(dict(completed=len(results),failed=0),indent=2)+'\n');print('Completed',len(results),'alignment audits',flush=True)
