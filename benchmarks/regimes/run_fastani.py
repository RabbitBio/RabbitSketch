#!/usr/bin/env python3
from pathlib import Path
import datetime, fcntl, hashlib, json, os, shutil, subprocess, time
ROOT=Path('/home/bigssd/zt/rabbitsketch2_benchmark')
out=ROOT/'runs/regimes_20261001/fastani';out.mkdir(parents=True,exist_ok=True)
inp=ROOT/'data/regimes_20261001/fastani_genomes.txt'
lock=(ROOT/'runs/.active_benchmark.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
cmd=[shutil.which('fastANI'),'--ql',str(inp),'--rl',str(inp),'-t','64','--minFraction','0.2','-o',str(out/'references.tsv')]
meta=dict(command=cmd,version=subprocess.run([cmd[0],'--version'],capture_output=True,text=True).stderr.strip(),input_list_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),threads=64)
if (out/'completion.json').exists():raise SystemExit('already completed')
(out/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
start=time.monotonic()
with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
    ret=subprocess.run(cmd,stdout=stdout,stderr=stderr,env=dict(os.environ,OMP_NUM_THREADS='64',OMP_DYNAMIC='FALSE')).returncode
meta.update(returncode=ret,wall_seconds=time.monotonic()-start,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
if ret:raise RuntimeError(meta)
meta.update(rows=sum(1 for _ in (out/'references.tsv').open()),sha256=hashlib.sha256((out/'references.tsv').read_bytes()).hexdigest())
(out/'completion.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta),flush=True)
