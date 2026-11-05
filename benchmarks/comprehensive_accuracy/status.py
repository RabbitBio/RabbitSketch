#!/usr/bin/env python3
"""Read-only progress; completed archives are separate from partial active output."""
import json
from pathlib import Path

RUN = Path('/home/bigssd/zt/rabbitsketch2_benchmark/runs/comprehensive_accuracy_20261030')
plan = json.loads((RUN / 'plan.json').read_text())
synthetic = [json.loads(p.read_text()) for p in (RUN / 'synthetic').glob('*/completion.json')]
real = list((RUN / 'real').glob('*/completion.json'))
active = json.loads((RUN / 'active_job.json').read_text()) if (RUN / 'active_job.json').exists() else None
if active and any((RUN / stage / active['job'] / 'completion.json').exists() for stage in ['synthetic', 'real']):
    active = None
result = dict(synthetic_completed_batches=len(synthetic), synthetic_planned_batches=12,
    synthetic_completed_measurements=sum(r['synthetic_rows'] for r in synthetic),
    synthetic_planned_measurements=plan['synthetic_measurements'],
    real_completed_jobs=len(real), real_planned_jobs=len(plan['real_jobs']), active_job=active,
    publication_verified=(RUN / 'complete.json').exists())
print(json.dumps(result, indent=2))
