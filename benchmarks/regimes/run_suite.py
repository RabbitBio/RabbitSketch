#!/usr/bin/env python3
"""Serial benchmark scheduler; each job is the existing main with 64 threads."""
from pathlib import Path
import argparse, csv, datetime, fcntl, hashlib, json, os, shutil, subprocess, time
from concurrent.futures import ThreadPoolExecutor
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark');DATA=ROOT/'data/regimes_20261001';HERE=Path(__file__).resolve().parent
OUT=ROOT/'runs/regimes_20261001';BINARY=ROOT/'build_validated/rabbitsketch_pilot'
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def atomic(path,value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)
def plan():
    jobs=[]
    def add(phase,algorithm,p,k,seed,objects,pairs,mode='accuracy',**extra):
        j=dict(phase=phase,algorithm=algorithm,parameter=p,k=k,seed=seed,objects=str(DATA/(objects+'.tsv')),pairs=str(DATA/(pairs+'.tsv')) if pairs else '-',mode=mode,**extra)
        j['job_id']=hashlib.sha256(json.dumps(j,sort_keys=True).encode()).hexdigest()[:24];jobs.append(j)
    grid=dict(minhash=[128,512,2048],bindash=[512,2048,8192],setsketch=[10,12,14],fracminhash=[500,1000,2000,5000,10000,20000,50000],kssd=[2,3])
    for phase,objects,pairs in [('ani','ani_objects','ani_pairs'),('cardinality','card_objects',None),('containment','contain_objects','contain_pairs')]:
        for k in [20,24]:add(phase,'exact',1,k,42,objects,pairs)
        configs=dict(hll=[8,10,12,14,16],setsketch=[8,10,12,14,16]) if phase=='cardinality' else grid
        if phase=='containment':configs=dict(configs,hll=[10,12,14])
        for algorithm,parameters in configs.items():
            for p in parameters:
                for k in [20,24]:
                    for seed in ([42] if algorithm=='kssd' else [42,101,202]):add(phase,algorithm,p,k,seed,objects,pairs)
    # Actual unique-genome all-vs-all, no diagonal and no replicated sketches.
    for n in [256,1024,4096,10000]:
        for algorithm,p in [('minhash',512),('bindash',2048),('kssd',3),('setsketch',12),('fracminhash',10000),('minhash',2048),('bindash',8192),('kssd',2)]:
            add('workload',algorithm,p,24,42,'workload_objects',None,mode='allpairs',limit=n,passes=3)
    for p in [128,512,2048,8192]:
        for seed in range(42,62):add('weighted','probminhash',p,24,seed,'weighted_objects','weighted_pairs',mode='weighted')
    return jobs
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phases',nargs='+',default=['ani','cardinality','containment','workload','weighted']);args=ap.parse_args()
    OUT.mkdir(exist_ok=True);(OUT/'jobs').mkdir(exist_ok=True)
    lock=(ROOT/'runs/.active_benchmark.lock').open('a');print('Waiting for benchmark lock',flush=True);fcntl.flock(lock,fcntl.LOCK_EX)
    jobs=plan();snapshot=OUT/'code_snapshot'
    if not snapshot.exists():
        snapshot.mkdir();shutil.copy2(BINARY,snapshot/'rabbitsketch_pilot')
        for p in HERE.glob('*'):
            if p.is_file():shutil.copy2(p,snapshot/p.name)
        for name in ['main_batch.cpp','real_input.h','CMakeLists.txt'] :shutil.copy2(HERE.parent/'screening'/name,snapshot/name)
        for p in (HERE.parent/'screening/patches').glob('*.patch'):shutil.copy2(p,snapshot/p.name)
    binary=snapshot/'rabbitsketch_pilot'
    file_inventory=OUT/'input_sha256.json'
    if not file_inventory.exists():
        paths=set()
        for manifest in {j['objects'] for j in jobs}:
            paths.update(r['path'] for r in csv.DictReader(open(manifest),delimiter='\t'))
        ordered=sorted(paths)
        print('Hashing',len(ordered),'input files',flush=True)
        with ThreadPoolExecutor(max_workers=16) as pool: hashes=dict(zip(ordered,pool.map(sha,ordered)))
        atomic(file_inventory,hashes)
    source=ROOT/'source/RabbitSketch-validated'
    source_commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    source_diff=subprocess.check_output(['git','-C',str(source),'diff','--binary'])
    (snapshot/'native_source.patch').write_bytes(source_diff)
    identity=dict(jobs=jobs,binary_sha256=sha(binary),threads=64,data_manifests={str(p):sha(p) for p in sorted(DATA.glob('*.tsv'))},input_sha_manifest=sha(file_inventory),source_commit=source_commit,source_diff_sha256=hashlib.sha256(source_diff).hexdigest(),script_sha256=sha(__file__))
    if (OUT/'suite_manifest.json').exists():
        if json.loads((OUT/'suite_manifest.json').read_text())!=identity:raise RuntimeError('suite identity changed')
    else:atomic(OUT/'suite_manifest.json',identity)
    env=dict(os.environ,OMP_NUM_THREADS='64',OMP_THREAD_LIMIT='64',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='spread',OMP_PLACES='threads',OPENBLAS_NUM_THREADS='1')
    selected=[j for j in jobs if j['phase'] in args.phases]
    for idx,j in enumerate(selected):
        d=OUT/'jobs'/j['job_id'];d.mkdir(exist_ok=True)
        if (d/'completion.json').exists():continue
        conf=d/'config.tsv';conf.write_text('#regimes-v1\n'+''.join(f'{k}\t{v}\n' for k,v in j.items() if k not in ['phase','job_id']))
        atomic(OUT/'progress.json',dict(status='running',job=j,position=idx+1,selected_jobs=len(selected),utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        print(f'{idx+1}/{len(selected)} {j["phase"]} {j["algorithm"]} p={j["parameter"]} k={j["k"]} seed={j["seed"]} N={j.get("limit","")}',flush=True)
        start=time.monotonic();cmd=[str(binary),'--screening-manifest',str(conf),'-t','64']
        with (d/'raw.tmp').open('w') as raw,(d/'stderr.log').open('w') as err:
            ret=subprocess.run(cmd,stdout=raw,stderr=err,env=env).returncode
        if ret:atomic(OUT/'progress.json',dict(status='failed',job=j,returncode=ret));raise RuntimeError(f'job failed: {d}')
        counts={};invalid=0
        with (d/'raw.tmp').open() as f:
            for line in f:
                if not line.startswith('{'):continue
                r=json.loads(line);counts[r['type']]=counts.get(r['type'],0)+1
                if r['type']=='allpairs':invalid+=r['invalid']
        expected_objects=int(j.get('limit',sum(1 for _ in open(j['objects']))-1));expected_pairs=0 if j['pairs']=='-' else sum(1 for _ in open(j['pairs']))-1
        if counts.get('object',0)!=expected_objects or counts.get('pair',0)!=expected_pairs or invalid:raise RuntimeError(f'invalid row counts {counts}')
        if j['mode']=='allpairs' and counts.get('allpairs')!=j['passes']:raise RuntimeError('missing all-pairs pass')
        (d/'raw.tmp').replace(d/'raw.jsonl')
        atomic(d/'completion.json',dict(job=j,command=cmd,wall_seconds=time.monotonic()-start,rows=counts,sha256=sha(d/'raw.jsonl')))
    atomic(OUT/'progress.json',dict(status='completed_selected_phases',phases=args.phases,jobs=len(selected),utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print('Completed phases',args.phases,flush=True)
if __name__=='__main__':main()
