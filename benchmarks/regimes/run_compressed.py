#!/usr/bin/env python3
"""Match KSSD's small payload with compact BinDash; no general grid expansion."""
import hashlib,json
import run_suite as run
def plan():
    jobs=[]
    def add(phase,a,p,k,seed,mode='accuracy',**extra):
        j=dict(phase=phase,algorithm=a,parameter=p,k=k,seed=seed,objects=str(run.DATA/('workload_objects.tsv' if phase=='workload' else 'ani_objects.tsv')),pairs='-' if phase=='workload' else str(run.DATA/'ani_pairs.tsv'),mode=mode,**extra)
        j['job_id']=hashlib.sha256(json.dumps(j,sort_keys=True).encode()).hexdigest()[:24];jobs.append(j)
    for k in [20,24]:
        for seed in [42,101,202]:add('ani','bindash',128,k,seed)
    for n in [1024,10000]:
        for a,p in [('minhash',128),('bindash',512),('bindash',128)]:add('workload',a,p,24,42,'allpairs',limit=n,passes=3)
    return jobs
run.OUT=run.ROOT/'runs/compressed_validation_20261001'
run.plan=plan
if __name__=='__main__':run.main()
