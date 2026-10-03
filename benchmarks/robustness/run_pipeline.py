#!/usr/bin/env python3
"""One reference/benchmark process at a time; resumable preparation prerequisites."""
import subprocess,sys,time
from common import *

def waitfor(name):
    while not (DATA/name).exists():time.sleep(10)
def command(script,*args):
    print('Starting',script,*args,flush=True)
    subprocess.run([sys.executable,str(HERE/script),*args],check=True)

if __name__=='__main__':
    waitfor('mixed_preparation.json');command('run_references.py','mixed')
    waitfor('natural_preparation.json');command('run_references.py','natural')
    command('run_native.py');command('run_alignment.py')
    atomic(RUN/'pipeline_completed.json',dict(status='raw_experiments_complete',native_progress=sha(RUN/'progress.json'),alignment_completion=sha(RUN/'alignment/completion.json')))
