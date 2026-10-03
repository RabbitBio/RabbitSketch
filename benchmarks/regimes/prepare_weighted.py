#!/usr/bin/env python3
"""Exact probability-Jaccard truth and reproducible positive-weight inputs."""
from pathlib import Path
import csv, hashlib, json
import numpy as np
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark')

def probability_jaccard(a,b):
    a=np.asarray(a,dtype=np.longdouble);b=np.asarray(b,dtype=np.longdouble)
    good=(a>0)&(b>0)
    ratio=np.divide(a,b,out=np.full_like(a,np.inf),where=b>0)
    order=np.argsort(ratio,kind='stable');ar=a[order];br=b[order];r=ratio[order]
    prefix=np.r_[np.longdouble(0),np.cumsum(br)]
    suffix=np.r_[np.cumsum(ar[::-1])[::-1],np.longdouble(0)]
    idx=np.searchsorted(r,ratio[good],side='left')
    den=b[good]*suffix[idx]+a[good]*prefix[idx]
    return float(np.sum(a[good]*b[good]/den)) if good.any() else 0.

def main():
    out=ROOT/'data/regimes_20261001';wd=out/'weights';wd.mkdir(exist_ok=True)
    objects=[];pairs=[];values={};provenance=[]
    def add(name,keys,w):
        p=wd/(name+'.tsv')
        with p.open('w') as f:
            for key,x in zip(keys,w):
                if x>0:f.write(f'{int(key)}\t{float(x):.17g}\n')
        objects.append(dict(object_id=name,path=str(p)));values[name]=np.array(w)
    def pair(a,b,**meta):
        x,y=values[a],values[b];support=float(((x>0)&(y>0)).sum()/((x>0)|(y>0)).sum())
        pairs.append(dict(pair_id=a+'__'+b,a=a,b=b,**meta,probability_jaccard=probability_jaccard(x,y),minmax_jaccard=float(np.minimum(x,y).sum()/np.maximum(x,y).sum()),support_jaccard=support))
    keys=np.arange(1,1025,dtype=np.uint64)*np.uint64(6364136223846793005)
    for replicate in range(20):
        rng=np.random.default_rng(2026100200+replicate)
        for model in ['uniform','lognormal','zipf0.5','zipf1.0','zipf1.5','zipf2.0']:
            if model=='uniform':a=np.ones(1024);target=rng.lognormal(0,1,1024)
            elif model=='lognormal':a=rng.lognormal(0,1,1024);target=rng.permutation(a)
            else:a=(np.arange(1,1025,dtype=float)**(-float(model[4:])))[rng.permutation(1024)];target=rng.permutation(a)
            a/=a.sum();target/=target.sum();name=f'w{replicate:02d}_{model}';add(name,keys,a)
            for change in [0.,.25,.5,.75,1.,10.]:
                b=a*10 if change==10 else (1-change)*a+change*target
                other=name+f'_{change:g}';add(other,keys,b)
                pair(name,other,family='synthetic_weights',model=model,replicate=replicate,perturbation=change,split='calibration' if replicate<10 else 'validation')
    source=Path('/home/bigssd/zt/cami2/marine/simulation_short_read')
    profiles=[]
    for i in range(10):
        p=source/f'abundance{i}.tsv';d={}
        for row in csv.reader(p.open(),delimiter='\t'):
            if row[0] in d:raise ValueError('duplicate CAMI key')
            d[row[0]]=float(row[1])
        profiles.append(d);provenance.append(dict(sample=i,path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),positive=sum(v>0 for v in d.values()),total_weight=sum(d.values())))
    names=sorted(set().union(*profiles));camikeys=[int.from_bytes(hashlib.sha256(n.encode()).digest()[:8],'little') for n in names]
    if len(set(camikeys))!=len(names):raise ValueError('profile key collision')
    for i,d in enumerate(profiles):add(f'cami{i}',camikeys,[d.get(n,0) for n in names])
    for i in range(10):
        for j in range(i+1,10):pair(f'cami{i}',f'cami{j}',family='cami_marine_profile',model='official_abundance',replicate=i,perturbation=-1,split='validation')
    for name,rows in [('weighted_objects',objects),('weighted_pairs',pairs),('cami_profile_provenance',provenance)]:
        with (out/(name+'.tsv')).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (out/'cami_key_map.json').write_text(json.dumps(dict(zip(names,camikeys)),indent=2)+'\n')
    print('weighted objects',len(objects),'pairs',len(pairs),flush=True)
if __name__=='__main__':main()
