#!/usr/bin/env python3
"""Reuse the existing main adapter, frozen native binary and sequential scheduler."""
import importlib.util,sys,shutil
from common import *

def plan():
    jobs=[]
    for phase in ['natural','mixed']:
        for k in [24,20]:
            configs=[('exact',1,{})]
            grid=([('minhash',p,{}) for p in [128,512,2048]]+
                  [('bindash',p,{}) for p in [128,512,2048,8192]]+
                  [('setsketch',12,{})]+[('fracminhash',p,{}) for p in [5000,10000]]+
                  [('kssd',3,{'half_subk':6}),('kssd',4,{'half_subk':7})]) if k==24 else [
                  ('minhash',512,{}),('bindash',2048,{}),('setsketch',12,{}),('fracminhash',10000,{}),('kssd',3,{'half_subk':6})]
            configs+=grid
            for algo,p,extra in configs:
                for seed in ([42] if algo in ['exact','kssd'] else [42,101,202]):
                    j=dict(phase=phase,algorithm=algo,parameter=p,k=k,seed=seed,objects=str(DATA/(phase+'_objects.tsv')),pairs=str(DATA/(phase+'_pairs.tsv')),mode='accuracy',**extra)
                    j['job_id']=hashlib.sha256(json.dumps(j,sort_keys=True).encode()).hexdigest()[:24];jobs.append(j)
    assert len(jobs)==94
    return jobs

def main():
    spec=importlib.util.spec_from_file_location('regimes_runner',HERE.parent/'regimes/run_suite.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    srcnames=['PROTOCOL.md','common.py','prepare_natural.py','prepare_mixed.py','run_references.py','run_native.py','run_alignment.py','run_pipeline.py','test_robustness.py']
    identity={n:sha(HERE/n) for n in srcnames};f=RUN/'preparation_code_manifest.json'
    if f.exists():assert json.loads(f.read_text())==identity,'preparation code changed'
    else:
        atomic(f,identity);dest=RUN/'preparation_code_snapshot';dest.mkdir(exist_ok=True)
        for n in srcnames:shutil.copy2(HERE/n,dest/n)
    module.DATA=DATA;module.OUT=RUN;module.plan=plan
    sys.argv=[__file__,'--phases','natural','mixed'];module.main()

if __name__=='__main__':main()
