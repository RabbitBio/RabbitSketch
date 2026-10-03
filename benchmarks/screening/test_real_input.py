import json
import os
from pathlib import Path
import random
import subprocess

import pytest


def words(records, k):
    result=set()
    for record in records:
        for i in range(len(record)-k+1):
            w=record[i:i+k].upper()
            if set(w)<=set('ACGT'):
                reverse=w.translate(str.maketrans('ACGT','TGCA'))[::-1]
                result.add(min(w,reverse))
    return result


@pytest.mark.parametrize('k',[1,4,20,24,32])
def test_exact_fasta_against_independent_naive_words(tmp_path,k):
    binary=os.environ['SCREENING_MAIN_BINARY']
    rng=random.Random(742)
    a=[''.join(rng.choice('ACGT') for _ in range(70)), 'acgtNNACGT', 'ACGT'*15]
    b=[a[0][:45]+'R'+'TGCA'*10, 'A'*50, a[2][::-1]]
    for name,records in [('a',a),('b',b)]:
        (tmp_path/(name+'.fa')).write_text(''.join(f'>record {i}\n{s[:33]}\n{s[33:]}\n' for i,s in enumerate(records)))
    manifest=tmp_path/'jobs.tsv'
    manifest.write_text('\t'.join(map(str,['0'*24,'exact',1,k,1,1,0,0,42,1,6,16,
                                         tmp_path/'a.fa',tmp_path/'b.fa']))+'\n')
    p=subprocess.run([binary,'--screening-manifest',str(manifest),'-t','2'],capture_output=True,text=True,check=True)
    row=[json.loads(line) for line in p.stdout.splitlines() if line.startswith('{') and 'job_id' in line][0]
    wa,wb=words(a,k),words(b,k)
    assert (row['n_a'],row['n_b'],row['intersection'])==(len(wa),len(wb),len(wa&wb))
    assert row['ground_truth']==len(wa&wb)/len(wa|wb)
    assert row['measurement_profile']=='main_openmp_real_pair'


@pytest.mark.parametrize('algorithm,parameter',[('minhash',64),('fracminhash',1),('hll',8),
                         ('setsketch',8),('probminhash',64),('bindash',64)])
def test_native_backends_consume_the_same_fasta(tmp_path,algorithm,parameter):
    rng=random.Random(946)
    a=[''.join(rng.choice('ACGT') for _ in range(100)), 'ACGT'*20]
    b=[a[0], 'TGCA'*20]
    for name,records in [('a',a),('b',b)]:
        (tmp_path/(name+'.fa')).write_text(''.join('>contig\n'+s+'\n' for s in records))
    wa,wb=words(a,20),words(b,20)
    cmd=[os.environ['SCREENING_BINARY'],'--algorithm',algorithm,'--parameter',str(parameter),
         '--k','20','--n-a',str(len(wa)),'--n-b',str(len(wb)),'--intersection',str(len(wa&wb)),
         '--fasta-a',str(tmp_path/'a.fa'),'--fasta-b',str(tmp_path/'b.fa'),
         '--query-min-ms','1','--verify-data','1']
    row=json.loads(subprocess.check_output(cmd,text=True))
    assert row['status']=='ok'
    assert row['verified_intersection']==len(wa&wb)
    assert 0<=row['estimate']<=1
    if algorithm=='fracminhash':
        assert row['estimate']==pytest.approx(len(wa&wb)/len(wa|wb),abs=1e-14)
