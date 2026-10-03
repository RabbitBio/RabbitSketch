#!/usr/bin/env python3
"""Known edits, preserved gap counts, and replayable compressed event scripts."""
import numpy as np
from common import *

LETTERS=np.array(list(b'ACGT'),dtype=np.uint8)
LOOK=np.zeros(256,dtype=np.uint8);LOOK[LETTERS]=np.arange(4)

def event_lengths(total,rng):
    if not total:return np.array([],dtype=np.int32)
    vals=np.minimum(rng.geometric(.5,size=total),10).astype(np.int32)
    stop=np.searchsorted(np.cumsum(vals),total)+1;vals=vals[:stop]
    vals[-1]-=int(vals.sum())-total
    return vals

def mutate(seq,target,indel_share,seed):
    rng=np.random.default_rng(seed);a=np.frombuffer(seq,dtype=np.uint8);valid=np.isin(a,LETTERS)
    nvalid=int(valid.sum());edited=round((1-target)*nvalid/(1-(1-target)*indel_share/2))
    ni=round(edited*indel_share/2);nd=ni;ns=edited-ni-nd
    dl=event_lengths(nd,rng);occupied=~valid.copy();dp=[]
    # All deletion spans are valid original coordinates and are disjoint.
    for length in dl:
        while True:
            pos=int(rng.integers(0,len(a)-int(length)+1))
            if not occupied[pos:pos+length].any():break
        occupied[pos:pos+length]=True;dp.append(pos)
    dp=np.array(dp,dtype=np.int64);order=np.argsort(dp);dp=dp[order];dl=dl[order]
    il=event_lengths(ni,rng);ip=np.sort(rng.choice(np.flatnonzero(~occupied),len(il),replace=False))
    occupied[ip]=True
    sp=np.sort(rng.choice(np.flatnonzero(~occupied),ns,replace=False))
    sv=LETTERS[(LOOK[a[sp]]+rng.integers(1,4,len(sp)))%4]
    ib=LETTERS[rng.integers(0,4,ni)]
    edits=dict(sub_pos=sp,sub_value=sv,del_pos=dp,del_len=dl,ins_pos=ip,ins_len=il,ins_bases=ib)
    b=replay(seq,edits)
    matches=nvalid-nd-ns
    truth=dict(valid_reference_bases=nvalid,reference_bases=len(a),query_bases=len(b),matches=matches,substitutions=ns,insertions=ni,deletions=nd,indel_events=len(il)+len(dl),insertion_events=len(il),deletion_events=len(dl),column_identity=matches/(nvalid+ni),homologous_ungapped_identity=matches/(nvalid-nd),edit_burden=(ns+ni+nd)/nvalid,achieved_indel_base_share=(ni+nd)/edited if edited else 0)
    validate_edits(seq,b,edits,truth)
    return b,edits,truth

def replay(seq,e):
    a=np.frombuffer(seq,dtype=np.uint8).copy();a[e['sub_pos']]=e['sub_value']
    deleted=np.zeros(len(a),dtype=bool)
    for p,n in zip(e['del_pos'],e['del_len']):deleted[p:p+n]=True
    keep=np.flatnonzero(~deleted)
    lengths=np.zeros(len(a),dtype=np.int32);lengths[e['ins_pos']]=e['ins_len']
    dest=np.arange(len(keep))+np.cumsum(lengths[keep]);b=np.empty(len(keep)+len(e['ins_bases']),dtype=np.uint8);b[dest]=a[keep]
    if len(e['ins_pos']):
        starts=np.searchsorted(keep,e['ins_pos'])+np.cumsum(e['ins_len'])-e['ins_len']
        offsets=np.arange(len(e['ins_bases']))-np.repeat(np.cumsum(e['ins_len'])-e['ins_len'],e['ins_len'])
        b[np.repeat(starts,e['ins_len'])+offsets]=e['ins_bases']
    return b.tobytes()

def validate_edits(seq,b,e,t):
    a=np.frombuffer(seq,dtype=np.uint8);occupied=np.zeros(len(a),dtype=np.uint8)
    assert np.all(a[e['sub_pos']]!=e['sub_value'])
    assert np.all(np.isin(a[e['sub_pos']],LETTERS))
    for p,n in zip(e['del_pos'],e['del_len']):
        assert n>=1 and n<=10 and p+n<=len(a)
        assert np.all(np.isin(a[p:p+n],LETTERS))
        occupied[p:p+n]+=1
    occupied[e['ins_pos']]+=1;occupied[e['sub_pos']]+=1
    assert occupied.max(initial=0)<=1
    assert len(b)==len(a)+int(e['ins_len'].sum())-int(e['del_len'].sum())
    assert len(e['ins_bases'])==int(e['ins_len'].sum())
    assert b.count(b'N')==seq.count(b'N')
    assert t['matches']+t['substitutions']+t['deletions']==t['valid_reference_bases']
    assert abs(t['column_identity']-t['matches']/(t['matches']+t['substitutions']+t['insertions']+t['deletions']))<1e-12

def save_fasta(path,oid,seq):
    with path.open('wb') as f:
        f.write(b'>'+oid.encode()+b'\n')
        f.write(b'\n'.join(seq[i:i+100] for i in range(0,len(seq),100))+b'\n')

def main():
    if (DATA/'mixed_preparation.json').exists():print('Mixed preparation already complete');return
    old=ROOT/'data/regimes_20261001'
    refs=read_tsv(old/'controlled_provenance.tsv');oldobjects={r['object_id']:r['path'] for r in read_tsv(old/'ani_objects.tsv')}
    gen=DATA/'mixed_genomes';gen.mkdir(exist_ok=True);scripts=DATA/'mixed_edits';scripts.mkdir(exist_ok=True)
    objects=[];pairs=[];auditrefs=set(r['reference'] for r in sorted(refs,key=lambda r:key(r['reference']))[:5])
    for r in refs:
        rid=r['reference'];seq=fasta_records(Path(oldobjects[rid]).read_bytes())[0][1]
        oid='mixed_'+rid;path=gen/(oid+'.fna');save_fasta(path,oid,seq);objects.append(dict(object_id=oid,path=str(path)))
        for target in [.90,.95,.97,.99]:
            for q in [0,.1,.25,.5]:
                pid=f'{oid}_id{round(target*100)}_indel{round(q*100)}';seed=int(key(pid)[:16],16)
                output=gen/(pid+'.fna');ep=scripts/(pid+'.npz');meta=scripts/(pid+'.json')
                if meta.exists():
                    truth=json.loads(meta.read_text());assert truth['output_sha256']==sha(output) and truth['edits_sha256']==sha(ep)
                else:
                    b,edits,truth=mutate(seq,target,q,seed);save_fasta(output,pid,b);np.savez_compressed(ep,**edits)
                    with np.load(ep) as saved:assert replay(seq,saved)==b
                    truth.update(output_sha256=sha(output),edits_sha256=sha(ep),reference_sha256=sha(path));atomic(meta,truth)
                objects.append(dict(object_id=pid,path=str(output)))
                pairs.append(dict(pair_id=pid,a=oid,b=pid,reference=rid,genus=r['genus'],target_identity=target,indel_share=q,mutation_model=f'sub{round((1-q)*100)}_indel{round(q*100)}',data_seed=seed,alignment_audit=int(rid in auditrefs),**truth))
        print(rid,r['genus'],'16 mutations validated',flush=True)
    tsv(DATA/'mixed_objects.tsv',objects);tsv(DATA/'mixed_pairs.tsv',pairs)
    (DATA/'mixed_refs.txt').write_text(''.join(r['path']+'\n' for r in objects if '_id' not in r['object_id']))
    (DATA/'mixed_queries.txt').write_text(''.join(r['path']+'\n' for r in objects if '_id' in r['object_id']))
    atomic(DATA/'mixed_preparation.json',dict(objects=len(objects),pairs=len(pairs),references=len(refs),audit_pairs=sum(p['alignment_audit'] for p in pairs),original_provenance_sha256=sha(old/'controlled_provenance.tsv'),protocol_sha256=sha(HERE/'PROTOCOL.md'),generator_sha256=sha(__file__)))

if __name__=='__main__':main()
