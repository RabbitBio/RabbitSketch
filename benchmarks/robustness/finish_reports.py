#!/usr/bin/env python3
"""Finish a native analysis once the separately running reference audit ends."""
import subprocess,sys,time
import pandas as pd
from common import *
from analyze import alignment,figures

if __name__=='__main__':
    status=json.loads((ANALYSIS/'analysis_status.json').read_text());assert status['native_statistics_complete']
    while not (RUN/'pipeline_completed.json').exists():time.sleep(10)
    done=json.loads((RUN/'alignment/completion.json').read_text());assert done['failed']==0 and done['total']==140
    z=pd.read_csv(ANALYSIS/'all_native_rows.csv.gz',low_memory=False);base=pd.read_csv(ANALYSIS/'exact_truth_and_models.csv',low_memory=False)
    alignment(z,base);figures(pd.read_csv(ANALYSIS/'budget_accuracy.csv'),base)
    status['alignment_complete']=True;atomic(ANALYSIS/'analysis_status.json',status)
    subprocess.run([sys.executable,str(HERE/'write_report.py')],check=True)
    print('Reference audit integrated; report ready for review',flush=True)
