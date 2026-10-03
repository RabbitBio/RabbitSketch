#!/usr/bin/env python3
"""Independent alignment audit; pairs selected without any sketch errors."""
import fcntl,re,subprocess,time
from common import *
TOOLS=ROOT/'tools/mummer-4.0.1-install/bin'

def parse_report(text):
    block=re.split(r'\n1-to-1\s',text)[1].split('\nM-to-M')[0]
    identity=re.search(r'AvgIdentity\s+([\d.]+)\s+([\d.]+)',block)
    lengths=re.search(r'TotalLength\s+(\d+)\s+(\d+)',block);totals=re.search(r'TotalBases\s+(\d+)\s+(\d+)',text)
    if not identity or not lengths or not totals:raise ValueError('missing one-to-one identity or lengths')
    x,y=[float(identity[i])/100 for i in [1,2]];a,b=[int(lengths[i])/int(totals[i]) for i in [1,2]]
    return dict(nucmer=(x+y)/2,nucmer_ref=x,nucmer_query=y,one_to_one_af_a=a,one_to_one_af_b=b,qualified=int(min(a,b)>=.2))

def main():
    out=RUN/'alignment';out.mkdir(exist_ok=True);objects={};plan=[]
    for panel in ['natural','mixed']:
        objects.update({r['object_id']:r['path'] for r in read_tsv(DATA/(panel+'_objects.tsv'))})
        plan += [dict(r,panel=panel) for r in read_tsv(DATA/(panel+'_alignment_pairs.tsv'))]
    tsv(out/'plan.tsv',plan)
    ident=dict(plan_sha256=sha(out/'plan.tsv'),binaries={n:sha(TOOLS/n) for n in ['nucmer','dnadiff']},source_sha256=sha(__file__),threads=64,options=['--maxmatch','-c','65','-l','20'])
    if (out/'manifest.json').exists():assert json.loads((out/'manifest.json').read_text())==ident
    else:atomic(out/'manifest.json',ident)
    lock=(ROOT/'runs/.active_benchmark.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
    results=[]
    for i,p in enumerate(plan):
        d=out/p['pair_id'];d.mkdir(exist_ok=True);result=d/'result.json';prefix=d/'align'
        if result.exists():results.append(json.loads(result.read_text()));continue
        print(f'{i+1}/{len(plan)} {p["panel"]} {p["pair_id"]}',flush=True);start=time.monotonic()
        cmds=[[str(TOOLS/'nucmer'),'--maxmatch','-c','65','-l','20','-t','64','-p',str(prefix),objects[p['a']],objects[p['b']]],[str(TOOLS/'dnadiff'),'-d',str(prefix)+'.delta','-p',str(prefix)]]
        try:
            for n,cmd in enumerate(cmds):
                with (d/f'command{n}.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=600)
            report=Path(str(prefix)+'.report');r=dict(pair_id=p['pair_id'],panel=p['panel'],status='ok',**parse_report(report.read_text()),report_sha256=sha(report))
        except (subprocess.SubprocessError,ValueError,IndexError) as e:r=dict(pair_id=p['pair_id'],panel=p['panel'],status='failed',error=str(e))
        r.update(commands=cmds,wall_seconds=time.monotonic()-start);atomic(result,r);results.append(r)
    atomic(out/'results.json',results);atomic(out/'completion.json',dict(total=len(results),failed=sum(r['status']!='ok' for r in results),results_sha256=sha(out/'results.json')))

if __name__=='__main__':main()
