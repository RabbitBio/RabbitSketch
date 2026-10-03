import numpy as np
import pytest
from prepare_mixed import mutate,replay
from common import fasta_records,fingerprint
from run_references import add_refs,balanced_take,ani_bin
from run_alignment import parse_report

@pytest.mark.parametrize('q',[0,.1,.25,.5])
@pytest.mark.parametrize('target',[.90,.95,.97,.99])
def test_mutation_identity_and_independent_replay(q,target):
    seq=b'ACGT'*2000+b'NNNN'+b'TGCA'*2000
    b,e,t=mutate(seq,target,q,301)
    # Independent, string-based event application; does not share vectorized replay.
    deleted={int(p)+i for p,n in zip(e['del_pos'],e['del_len']) for i in range(int(n))}
    subs={int(p):bytes([int(v)]) for p,v in zip(e['sub_pos'],e['sub_value'])}
    ins={};off=0
    for p,n in zip(e['ins_pos'],e['ins_len']):ins[int(p)]=e['ins_bases'][off:off+n].tobytes();off+=n
    expected=b''.join(ins.get(i,b'')+subs.get(i,seq[i:i+1]) for i in range(len(seq)) if i not in deleted)
    assert b==expected and replay(seq,e)==expected
    assert abs(t['column_identity']-target)<2/t['valid_reference_bases']
    assert abs(t['achieved_indel_base_share']-q)<3/(t['substitutions']+t['insertions']+t['deletions'])
    assert mutate(seq,target,q,301)[0]==b
    assert (t['homologous_ungapped_identity']>=t['column_identity'])

def test_fingerprint_normalization_and_order():
    a=fasta_records(b'>a\nACGAT\n>b\nNCT\n')
    b=fasta_records(b'>x\nAGr\n>y\nATCGT\n')
    assert fingerprint(a)==fingerprint(b)
    assert fingerprint(a)!=fingerprint(fasta_records(b'>c\nACGATNCT\n'))

def test_directional_reference_missing_not_zero():
    p=[dict(pair_id='p',a='a',b='b')];obj={'a':'A','b':'B'}
    assert add_refs(p,obj,{('A','B'):(.99,.8)})[0]['fastani']==''
    assert add_refs(p,obj,{('A','B'):(.99,.8),('B','A'):(.97,.1)})[0]['reference_available']==0
    r=add_refs(p,obj,{('A','B'):(.99,.8),('B','A'):(.97,.7)})[0]
    assert r['fastani']==pytest.approx(.98) and r['ani_bin']=='95-99'
    assert [ani_bin(x) for x in [.90,.95,.99,1.0,.899,None]]==['90-95','95-99','99-100','99-100','below90','missing']

def test_balance_does_not_duplicate_or_fill_from_other_bins():
    rows=[dict(pair_id=str(i),genus_a='a' if i<10 else 'b') for i in range(12)]
    chosen=balanced_take(rows,4)
    assert len({r['pair_id'] for r in chosen})==4
    assert sum(r['genus_a']=='a' for r in chosen)==2
    assert len(balanced_take(rows,100))==12

def test_alignment_one_to_one_coverage_excludes_short_matches():
    txt='TotalBases 1000 2000\n1-to-1 Alignments\nTotalLength 300 300\nAvgIdentity 99.0 98.0\nM-to-M Alignments\n'
    r=parse_report(txt)
    assert r['nucmer']==pytest.approx(.985) and r['qualified']==0
