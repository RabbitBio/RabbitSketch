from pathlib import Path
import json, subprocess, sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).parent))
from prepare_weighted import probability_jaccard
BIN=Path('/home/bigssd/zt/rabbitsketch2_benchmark/build_validated/rabbitsketch_pilot')

def run(tmp_path,algorithm,parameter=128,mode='accuracy',weighted=False):
    if weighted:
        (tmp_path/'a').write_text('10\t1\n20\t3\n30\t7\n');(tmp_path/'b').write_text('10\t10\n20\t30\n30\t70\n')
    else:
        (tmp_path/'a').write_text('>a\nACGTACGTACGTACGTACGTACGTACGT\n>b\nAACTGGTTACCGTTACGTGACGTTACGT\n');(tmp_path/'b').write_bytes((tmp_path/'a').read_bytes())
    (tmp_path/'objects').write_bytes(f'object_id\tpath\r\na\t{tmp_path}/a\r\nb\t{tmp_path}/b\r\n'.encode())
    (tmp_path/'pairs').write_text('pair_id\ta\tb\np\ta\tb\n')
    (tmp_path/'config').write_text(f'#regimes-v1\nalgorithm\t{algorithm}\nparameter\t{parameter}\nk\t20\nseed\t42\nmode\t{mode}\nobjects\t{tmp_path}/objects\npairs\t{tmp_path}/pairs\npasses\t2\n')
    r=subprocess.run([str(BIN),'--screening-manifest',str(tmp_path/'config'),'-t','3'],capture_output=True,text=True)
    assert r.returncode==0,r.stderr
    return [json.loads(s) for s in r.stdout.splitlines() if s.startswith('{')]

@pytest.mark.parametrize('algorithm,parameter',[('minhash',128),('bindash',128),('bindash_hll',384),('setsketch',8),('hll',8),('fracminhash',1)])
def test_identical_cache_and_distinct_allpairs(tmp_path,algorithm,parameter):
    rows=run(tmp_path,algorithm,parameter)
    q=next(r for r in rows if r['type']=='pair');assert q['jaccard']==pytest.approx(1,abs=1e-8)
    rows=run(tmp_path,algorithm,parameter,mode='allpairs')
    pairs=[r for r in rows if r['type']=='allpairs'];assert len(pairs)==2
    assert all(r['pairs']==1 and r['checksum']==pytest.approx(1,abs=1e-8) and r['invalid']==0 for r in pairs)

def test_exact_contig_boundaries(tmp_path):
    rows=run(tmp_path,'exact',1)
    p=next(r for r in rows if r['type']=='pair')
    seqs=['ACGTACGTACGTACGTACGTACGTACGT','AACTGGTTACCGTTACGTGACGTTACGT'];words=set()
    for s in seqs:
        for i in range(len(s)-19):
            w=s[i:i+20];rc=w.translate(str.maketrans('ACGT','TGCA'))[::-1];words.add(min(w,rc))
    assert p['n_a']==p['n_b']==p['intersection']==len(words)

def test_weighted_scale_invariance(tmp_path):
    rows=run(tmp_path,'probminhash',4096,mode='weighted',weighted=True)
    assert next(r for r in rows if r['type']=='pair')['jaccard']==1
    assert probability_jaccard([1,3,7],[10,30,70])==pytest.approx(1)
    assert sum(np.minimum([1,3,7],[10,30,70]))/sum(np.maximum([1,3,7],[10,30,70]))==pytest.approx(.1)

def test_probability_jaccard_independent_formula():
    rng=np.random.default_rng(67)
    for _ in range(100):
        a=rng.integers(0,8,17).astype(float);b=rng.integers(0,8,17).astype(float)
        slow=sum(1/sum(np.maximum(a/a[i],b/b[i])) for i in range(len(a)) if a[i]>0 and b[i]>0)
        assert probability_jaccard(a,b)==pytest.approx(slow,abs=1e-13)
        assert probability_jaccard(a,b)==pytest.approx(probability_jaccard(b,a),abs=1e-13)
