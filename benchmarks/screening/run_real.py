#!/usr/bin/env python3
"""Independent real-genome validation through the existing main, one batch at a time."""
import argparse
import csv
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from run import atomic_json, canonical, digest, file_sha, source_identity, DEFAULT_ROOT, HERE
from run_main import FIELDS

PARAMETERS = {
    'minhash': [128, 512, 2048], 'fracminhash': [50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000],
    'hll': [10,12,14], 'setsketch': [10,12,14],
    'probminhash': [128,512,2048], 'bindash': [512,2048,8192],
}
KS=[20,24]
SEEDS=[42,101,202]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pairs',type=Path,required=True)
    ap.add_argument('--binary',type=Path,required=True)
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--threads',type=int,default=64)
    args=ap.parse_args()
    if not 1<=args.threads<=len(os.sched_getaffinity(0)): raise ValueError('invalid threads')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    lock=(DEFAULT_ROOT/'runs/.active_benchmark.lock').open('a')
    print('Waiting for the single benchmark process slot',flush=True)
    fcntl.flock(lock,fcntl.LOCK_EX)
    pairs=list(csv.DictReader(args.pairs.open(),delimiter='\t'))
    if len({p['pair_id'] for p in pairs})!=len(pairs): raise ValueError('duplicate pair IDs')
    paths=sorted({p[k] for p in pairs for k in ('fileA','fileB')})
    identity=dict(binary_sha256=file_sha(args.binary),source=source_identity(args.source),
                  pairs_sha256=file_sha(args.pairs),genomes={p:file_sha(p) for p in paths},
                  script_sha256=file_sha(__file__),threads=args.threads,ks=KS,seeds=SEEDS,
                  parameters=PARAMETERS,measurement_profile='main_openmp_real_pair',
                  truth='exact distinct canonical 2-bit k-mers; record boundaries separated by N',
                  reuse='each pair builds its own two sketches; no cross-pair build cache')
    if (out/'manifest.json').exists():
        if json.loads((out/'manifest.json').read_text())!=identity: raise ValueError('identity changed')
    else: atomic_json(out/'manifest.json',identity)
    env=dict(os.environ,OMP_NUM_THREADS=str(args.threads),OMP_DYNAMIC='FALSE',
             OMP_THREAD_LIMIT=str(args.threads),OMP_PROC_BIND='spread',OMP_PLACES='threads',
             OPENBLAS_NUM_THREADS='1')
    batch_dir=out/'batches'; batch_dir.mkdir(exist_ok=True)
    count=0
    def execute(jobs, phase):
        nonlocal count
        todo=[]
        for job in jobs:
            rowpath=out/phase/(job['job_id']+'.json')
            if rowpath.exists():
                r=json.loads(rowpath.read_text())
                if r['status']!='ok' or any(r[k]!=v for k,v in job.items()):
                    raise ValueError('invalid saved real row')
            else: todo.append(job)
        if not todo: return
        bid=digest([phase,[j['job_id'] for j in todo]])[:24]
        tsv=batch_dir/(bid+'.tsv')
        tsv.write_text(''.join('\t'.join(str(j[k]) for k in FIELDS+('fasta_a','fasta_b'))+'\n' for j in todo))
        expected={j['job_id']:j for j in todo};received=set();start=time.monotonic()
        print(f'{phase}: {todo[0]["algorithm"]} p={todo[0]["parameter"]} k={todo[0]["k"]} seed={todo[0]["hash_seed"]}, {len(todo)} pairs',flush=True)
        cmd=[str(args.binary.resolve()),'--screening-manifest',str(tsv),'-t',str(args.threads)]
        with (batch_dir/(bid+'.stderr')).open('w') as err:
            proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err,text=True,env=env)
            try:
                for line in proc.stdout:
                    if not line.startswith('{'): err.write(line);continue
                    r=json.loads(line)
                    if 'batch_event' in r:
                        if r['threads']!=args.threads: raise ValueError('thread mismatch')
                        continue
                    jid=r['job_id']
                    if jid not in expected or jid in received: raise ValueError('unknown/duplicate row')
                    job=expected[jid]
                    if r['status']!='ok' or r['measurement_profile']!='main_openmp_real_pair':
                        raise ValueError(str(r))
                    if phase=='truth':
                        row=dict(job)
                        row.update(exact_n_a=r.pop('n_a'),exact_n_b=r.pop('n_b'),
                                   exact_intersection=r.pop('intersection'))
                        row.update(r)
                    else:
                        if abs(r['ground_truth']-job['ground_truth'])>1e-15: raise ValueError('truth mismatch')
                        if any(r['input_digest_'+s]!=job['input_digest_'+s] for s in ('a','b')):
                            raise ValueError('FASTA content differs from exact truth')
                        row=dict(job,**{k:v for k,v in r.items() if k not in job})
                    atomic_json(out/phase/(jid+'.json'),row);received.add(jid);count+=1
                    atomic_json(out/'progress.json',dict(phase=phase,new_completed=count,
                        current_algorithm=job['algorithm'],current_parameter=job['parameter'],
                        current_k=job['k'],current_seed=job['hash_seed'],batch_pid=proc.pid,
                        planned_measurements=len(pairs)*len(KS)*len(SEEDS)*sum(map(len,PARAMETERS.values()))))
                if proc.wait()!=0 or len(received)!=len(todo): raise RuntimeError('incomplete batch')
            finally:
                if proc.poll() is None: proc.kill();proc.wait()
        atomic_json(batch_dir/(bid+'.json'),dict(command=cmd,wall_seconds=time.monotonic()-start,
                                                pairs=len(todo),phase=phase,threads=args.threads))
    truths={}
    for k in KS:
        jobs=[]
        for p in pairs:
            j=dict(pair_id=p['pair_id'],tier=p['tier'],algorithm='exact',parameter=1,k=k,n_a=1,n_b=1,
                   intersection=0,data_seed=0,hash_seed=42,query_min_ms=5,half_subk=6,bbits=16,
                   fasta_a=p['fileA'],fasta_b=p['fileB'])
            j['job_id']=digest(j)[:24];jobs.append(j)
        execute(jobs,'truth')
        for j in jobs: truths[(j['pair_id'],k)]=json.loads((out/'truth'/(j['job_id']+'.json')).read_text())
    all_jobs=[]
    # Declared order; no parameters selected from measured errors.
    for algorithm,parameters in PARAMETERS.items():
        for parameter in parameters:
            for k in KS:
                for seed in SEEDS:
                    jobs=[]
                    for p in pairs:
                        t=truths[(p['pair_id'],k)]
                        j=dict(pair_id=p['pair_id'],tier=p['tier'],algorithm=algorithm,parameter=parameter,
                               k=k,n_a=t['exact_n_a'],n_b=t['exact_n_b'],intersection=t['exact_intersection'],
                               ground_truth=t['ground_truth'],data_seed=0,hash_seed=seed,query_min_ms=5,
                               half_subk=6,bbits=16,fasta_a=p['fileA'],fasta_b=p['fileB'],
                               input_digest_a=t['input_digest_a'],input_digest_b=t['input_digest_b'])
                        j['job_id']=digest(j)[:24];jobs.append(j)
                    all_jobs+=jobs
                    execute(jobs,'rows')
    atomic_json(out/'plan.json',dict(parameters=PARAMETERS,ks=KS,seeds=SEEDS,jobs=all_jobs))
    with (out/'raw.jsonl').open('w') as stream:
        for j in all_jobs: stream.write(canonical(json.loads((out/'rows'/(j['job_id']+'.json')).read_text()))+'\n')
    atomic_json(out/'progress.json',dict(phase='completed',completed=len(all_jobs),failed=0,
                unique_genomes=len(paths),pairs=len(pairs),threads=args.threads,
                finished_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
    print('Completed',len(all_jobs),'real-genome measurements',flush=True)

if __name__=='__main__': main()
