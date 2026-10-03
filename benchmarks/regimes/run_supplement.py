#!/usr/bin/env python3
"""Isolated task-specific supplements; leave the primary suite frozen."""
import hashlib,json,sys
import run_suite as run
def plan():
    jobs=[]
    def add(phase,algorithm,p,k,seed,objects,pairs,mode='accuracy',**extra):
        j=dict(phase=phase,algorithm=algorithm,parameter=p,k=k,seed=seed,objects=str(run.DATA/(objects+'.tsv')),pairs=str(run.DATA/(pairs+'.tsv')) if pairs else '-',mode=mode,**extra)
        j['job_id']=hashlib.sha256(json.dumps(j,sort_keys=True).encode()).hexdigest()[:24];jobs.append(j)
    for p in [384,1536,6144]:
        for k in [20,24]:
            for seed in [42,101,202]:add('containment','bindash_hll',p,k,seed,'contain_objects','contain_pairs')
    for phase,objects,pairs in [('ani','ani_objects','ani_pairs'),('containment','contain_objects','contain_pairs')]:
        for k in [20,24]:add(phase,'kssd',4,k,42,objects,pairs,half_subk=7)
    for n in [256,1024,4096,10000]:add('workload','kssd',4,24,42,'workload_objects',None,mode='allpairs',limit=n,passes=3,half_subk=7)
    return jobs
run.OUT=run.ROOT/'runs/regimes_supplement_20261001'
run.plan=plan
if __name__=='__main__':run.main()
