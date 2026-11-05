#!/usr/bin/env python3
"""Continue the dependency-ordered experiment after calibration has completed."""
import subprocess
import sys
import time

from scenario_common import *


def call(script, *arguments, log_name):
    with (RUN / log_name).open('w') as log:
        subprocess.run([sys.executable, str(HERE / script), *arguments], stdout=log, stderr=subprocess.STDOUT, check=True)


def main():
    started = time.monotonic()
    for marker in [RUN / 'development_complete.json', OUT / 'development_scoring_complete.json']:
        while not marker.exists():
            if time.monotonic() - started > 12 * 3600: raise TimeoutError('development did not complete')
            time.sleep(5)
    atomic(RUN / 'pipeline_status.json', dict(stage='fit_model'))
    call('advisor.py', '--stage', 'fit', log_name='fit.log')
    atomic(RUN / 'pipeline_status.json', dict(stage='validation_pilot'))
    call('run.py', '--stage', 'validation-pilot', log_name='validation_pilot.log')
    call('measure.py', '--split', 'validation', '--stage', 'features', log_name='validation_features.log')
    call('advisor.py', '--stage', 'freeze-decisions', log_name='freeze_decisions.log')
    atomic(RUN / 'pipeline_status.json', dict(stage='validation_outcomes'))
    with (RUN / 'validation_native.log').open('w') as log:
        native = subprocess.Popen([sys.executable, str(HERE / 'run.py'), '--stage', 'validation-outcomes'],
                                  stdout=log, stderr=subprocess.STDOUT)
        try:
            while not (RUN / 'validation-outcomes_plan.json').exists():
                if native.poll() is not None: raise RuntimeError('validation runner failed before its plan')
                time.sleep(1)
            call('measure.py', '--split', 'validation', '--stage', 'scores', '--wait', log_name='validation_scoring.log')
            if native.wait() != 0: raise RuntimeError('validation native runner failed')
        except BaseException:
            if native.poll() is None: native.terminate(); native.wait()
            raise
    atomic(RUN / 'pipeline_status.json', dict(stage='evaluate'))
    call('evaluate.py', log_name='evaluation.log')
    atomic(RUN / 'pipeline_status.json', dict(stage='evaluation_complete'))


if __name__ == '__main__':
    try: main()
    except BaseException as exc:
        atomic(RUN / 'pipeline_status.json', dict(stage='failed', error=repr(exc)))
        raise
