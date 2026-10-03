#!/usr/bin/env python3
import hashlib,json
import run_suite as run
def plan():
    jobs=[]
    def add(algorithm,p,k,seed,**extra):
        j=dict(phase='containment',algorithm=algorithm,parameter=p,k=k,seed=seed,objects=str(run.DATA/'contain_objects.tsv'),pairs=str(run.DATA/'contain_pairs.tsv'),mode='accuracy',**extra)
        j['job_id']=hashlib.sha256(json.dumps(j,sort_keys=True).encode()).hexdigest()[:24];jobs.append(j)
    for k in [20,24]:add('exact',1,k,42)
    grid=dict(minhash=[128,512,2048],bindash=[512,2048,8192],setsketch=[10,12,14],fracminhash=[500,1000,2000,5000,10000,20000,50000],kssd=[2,3,4],hll=[10,12,14],bindash_hll=[384,1536,6144])
    for a,ps in grid.items():
        for p in ps:
            for k in [20,24]:
                for seed in ([42] if a=='kssd' else [42,101,202]):add(a,p,k,seed,**({'half_subk':7} if a=='kssd' and p==4 else {}))
    return jobs
run.DATA=run.ROOT/'data/containment_holdout_20261001'
run.OUT=run.ROOT/'runs/containment_holdout_20261001'
run.plan=plan
if __name__=='__main__':run.main()
