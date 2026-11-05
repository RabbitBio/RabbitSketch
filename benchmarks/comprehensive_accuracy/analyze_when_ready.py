#!/usr/bin/env python3
"""Analyze newly completed panels while the single native runner continues."""
import os
import subprocess
import sys
import time
from common import *


def call(script, *args):
    env = os.environ.copy(); env['OPENBLAS_NUM_THREADS'] = '1'
    subprocess.run([sys.executable, str(HERE / script), *args], env=env, check=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    expected_real = len(json.loads((RUN / 'plan.json').read_text())['real_jobs'])
    last = (-1, -1)
    while True:
        current = (len(list((RUN / 'synthetic').glob('*/completion.json'))),
                   len(list((RUN / 'real').glob('*/completion.json'))))
        if current != last:
            if current[0] != last[0]: call('analyze_synthetic.py', '--completed-only')
            if current[1] != last[1] and current[1]: call('analyze_real.py', '--completed-only')
            last = current
            atomic(OUT / 'pipeline_status.json', dict(phase='waiting_for_native',
                synthetic_batches=current[0], real_jobs=current[1], pid=os.getpid()))
        if current == (12, expected_real): break
        time.sleep(30)
    call('analyze_synthetic.py')
    call('constraint_profiles.py')
    call('analyze_real.py')
    call('downstream.py')
    call('figures.py')
    atomic(OUT / 'pipeline_status.json', dict(phase='analysis_complete',
        synthetic_batches=12, real_jobs=expected_real, pid=os.getpid()))


if __name__ == '__main__':
    try: main()
    except Exception as exc:
        atomic(OUT / 'pipeline_error.json', dict(error=str(exc), pid=os.getpid()))
        raise
