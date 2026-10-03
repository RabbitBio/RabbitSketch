#!/usr/bin/env python3
"""Reference-only natural selection and directional FastANI, before native scoring."""
import argparse,collections,datetime,fcntl,shutil,subprocess,time
from common import *

def run_fastani(name,ql,rl):
    out=RUN/name;out.mkdir(exist_ok=True);exe=shutil.which('fastANI')
    cmd=[exe,'--ql',str(ql),'--rl',str(rl),'-t','64','--minFraction','0.2','-o',str(out/'references.tsv')]
    ident=dict(command=cmd,query_list_sha256=sha(ql),reference_list_sha256=sha(rl),binary_sha256=sha(exe))
    if (out/'completion.json').exists():
        done=json.loads((out/'completion.json').read_text());assert done['identity']==ident and done['sha256']==sha(out/'references.tsv');return
    lock=(ROOT/'runs/.active_benchmark.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
    atomic(out/'manifest.json',dict(identity=ident,version=subprocess.run([exe,'--version'],capture_output=True,text=True).stderr.strip(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print('FastANI',name,flush=True);start=time.monotonic()
    with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
        subprocess.run(cmd,stdout=stdout,stderr=stderr,check=True,env=dict(os.environ,OMP_NUM_THREADS='64',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1'))
    atomic(out/'completion.json',dict(identity=ident,wall_seconds=time.monotonic()-start,rows=sum(1 for _ in (out/'references.tsv').open()),sha256=sha(out/'references.tsv')))
    lock.close()

def read_fastani(names):
    obs={}
    for name in names:
        for line in (RUN/name/'references.tsv').open():
            a,b,x,m,n=line.rstrip().split('\t');obs[a,b]=(float(x)/100,int(m)/int(n))
    return obs

def add_refs(pairs,objects,obs):
    result=[]
    for p in pairs:
        a,b=objects[p['a']],objects[p['b']];x=obs.get((a,b));y=obs.get((b,a))
        ok=bool(x and y and min(x[1],y[1])>=.2)
        ident=(x[0]+y[0])/2 if ok else None
        label=ani_bin(ident) if ok else 'missing'
        result.append(dict(p,fastani=ident if ok else '',fastani_a=x[0] if x else '',fastani_b=y[0] if y else '',af_a=x[1] if x else '',af_b=y[1] if y else '',reference_available=int(ok),ani_bin=label))
    return result

def ani_bin(x):
    if x is None:return 'missing'
    for lo,hi,label in [(.90,.95,'90-95'),(.95,.99,'95-99'),(.99,1.000001,'99-100')]:
        if lo<=x<hi:return label
    return 'below90'

def balanced_take(rows,n):
    groups=collections.defaultdict(list)
    for r in sorted(rows,key=lambda r:key(r['pair_id'])):groups[r['genus_a']].append(r)
    groups={g:iter(v) for g,v in groups.items()};out=[]
    while groups and len(out)<n:
        for g in sorted(list(groups),key=key):
            try:out.append(next(groups[g]))
            except StopIteration:del groups[g]
            if len(out)==n:break
    return out

def natural():
    assert (DATA/'natural_preparation.json').exists()
    run_fastani('natural_fastani',DATA/'natural_fastani.txt',DATA/'natural_fastani.txt')
    if (DATA/'natural_selection.json').exists():return
    objects={r['object_id']:r['path'] for r in read_tsv(DATA/'natural_candidate_objects.tsv')}
    candidates=add_refs(read_tsv(DATA/'natural_candidate_pairs.tsv'),objects,read_fastani(['natural_fastani']))
    tsv(DATA/'natural_candidate_references.tsv',candidates)
    species=collections.defaultdict(list)
    for p in candidates:
        if p['pair_class']=='within_species':species[p['species_a']].append(p)
    selected={};reasons=collections.defaultdict(list)
    def add(rs,reason):
        for r in rs:selected[r['pair_id']]=r;reasons[r['pair_id']].append(reason)
    for rs in species.values():add(sorted(rs,key=lambda r:key(r['pair_id']))[:4],'species_coverage')
    available={}
    for label in ['90-95','95-99','99-100']:
        pool=[p for p in candidates if p['ani_bin']==label];available[label]=len(pool)
        add(balanced_take(pool,600),'ani_'+label)
    add([p for p in candidates if p['pair_class'].startswith('cross_')],'distant_control')
    final=[dict(r,selection_reason=';'.join(reasons[r['pair_id']])) for r in selected.values()];final.sort(key=lambda r:r['pair_id'])
    tsv(DATA/'natural_pairs.tsv',final);shutil.copy2(DATA/'natural_candidate_objects.tsv',DATA/'natural_objects.tsv')
    audit=[]
    for label in ['90-95','95-99','99-100']:
        audit += sorted([p for p in final if p['ani_bin']==label],key=lambda r:key('alignment:'+r['pair_id']))[:20]
    tsv(DATA/'natural_alignment_pairs.tsv',audit)
    atomic(DATA/'natural_selection.json',dict(pairs=len(final),species=len(species),candidate_bin_counts=available,final_bin_counts=dict(collections.Counter(p['ani_bin'] for p in final)),pair_classes=dict(collections.Counter(p['pair_class'] for p in final)),alignment_pairs=len(audit),selection_before_any_native_sketch=True,candidate_reference_sha256=sha(DATA/'natural_candidate_references.tsv'),pairs_sha256=sha(DATA/'natural_pairs.tsv')))
    print(json.dumps(json.loads((DATA/'natural_selection.json').read_text()),indent=2),flush=True)

def mixed():
    assert (DATA/'mixed_preparation.json').exists()
    run_fastani('mixed_forward_fastani',DATA/'mixed_queries.txt',DATA/'mixed_refs.txt')
    run_fastani('mixed_reverse_fastani',DATA/'mixed_refs.txt',DATA/'mixed_queries.txt')
    objects={r['object_id']:r['path'] for r in read_tsv(DATA/'mixed_objects.tsv')}
    refs=add_refs(read_tsv(DATA/'mixed_pairs.tsv'),objects,read_fastani(['mixed_forward_fastani','mixed_reverse_fastani']))
    tsv(DATA/'mixed_references.tsv',refs);tsv(DATA/'mixed_alignment_pairs.tsv',[p for p in refs if int(p['alignment_audit'])])
    atomic(DATA/'mixed_reference_summary.json',dict(pairs=len(refs),reference_available=sum(int(p['reference_available']) for p in refs),alignment_pairs=sum(int(p['alignment_audit']) for p in refs)))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('panel',choices=['natural','mixed']);a=ap.parse_args()
    globals()[a.panel]()
