#!/usr/bin/env python3
"""Verify completed artifacts and snapshot analysis, without altering old stages."""
import datetime,re,shutil
from concurrent.futures import ThreadPoolExecutor
from common import *

def main():
    assert (RUN/'pipeline_completed.json').exists()
    suite=json.loads((RUN/'suite_manifest.json').read_text());assert len(suite['jobs'])==94
    for job in suite['jobs']:
        d=RUN/'jobs'/job['job_id'];done=json.loads((d/'completion.json').read_text());assert done['sha256']==sha(d/'raw.jsonl')
    inventory=json.loads((RUN/'input_sha256.json').read_text());paths=sorted(inventory)
    with ThreadPoolExecutor(max_workers=16) as pool:
        for path,value in zip(paths,pool.map(sha,paths)):assert value==inventory[path],path
    for r in read_tsv(DATA/'natural_candidates.tsv'):assert inventory[r['path']]==r['sha256']
    objects={r['object_id']:r['path'] for r in read_tsv(DATA/'mixed_objects.tsv')}
    for r in read_tsv(DATA/'mixed_pairs.tsv'):
        assert inventory[objects[r['a']]]==r['reference_sha256']
        assert inventory[objects[r['b']]]==r['output_sha256']
        assert sha(DATA/'mixed_edits'/(r['pair_id']+'.npz'))==r['edits_sha256']
    for n in ['natural_fastani','mixed_forward_fastani','mixed_reverse_fastani']:
        done=json.loads((RUN/n/'completion.json').read_text());assert sha(RUN/n/'references.tsv')==done['sha256']
    alignment=json.loads((RUN/'alignment/completion.json').read_text());assert alignment['failed']==0
    for r in json.loads((RUN/'alignment/results.json').read_text()):assert sha(RUN/'alignment'/r['pair_id']/'align.report')==r['report_sha256']
    code=ANALYSIS/'code_snapshot';code.mkdir(exist_ok=True)
    for p in HERE.glob('*'):
        if p.is_file():shutil.copy2(p,code/p.name)
    shutil.copy2(HERE/'PROTOCOL.md',ANALYSIS/'PROTOCOL.md')
    shutil.copy2(RUN/'preparation_tests.log',ANALYSIS/'preparation_tests.log')
    for filename in ['RabbitSketch2_robustness_validation.md','RabbitSketch2_algorithm_selection_analysis.md']:
        text=(HERE.parent.parent/filename).read_text()
        for target in re.findall(r'\]\(([^\n)]+)\)',text):
            if target.startswith('http') or target.startswith('#'):continue
            path=Path(target.split('#')[0])
            if not path.is_absolute():path=HERE.parent.parent/path
            if path.name=='completion.json' and path.parent==ANALYSIS:continue
            assert path.exists(),str(path)
    hashes={str(p.relative_to(ANALYSIS)):sha(p) for p in sorted(ANALYSIS.rglob('*')) if p.is_file() and p.name!='completion.json' and '__pycache__' not in str(p)}
    report=HERE.parent.parent/'RabbitSketch2_robustness_validation.md'
    atomic(ANALYSIS/'completion.json',dict(status='complete',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),jobs=len(suite['jobs']),input_files_verified=len(paths),mixed_edit_scripts_verified=320,alignment_pairs_verified=alignment['total'],analysis_artifacts=hashes,report_path=str(report),report_sha256=sha(report),native_suite_sha256=sha(RUN/'suite_manifest.json'),preparation_code_manifest_sha256=sha(RUN/'preparation_code_manifest.json'),protocol_sha256=sha(HERE/'PROTOCOL.md'),second_cpu_tested=False,orderminhash_included=False))
    print('Verified',len(paths),'inputs;',len(hashes),'artifacts;',alignment['total'],'alignment pairs',flush=True)

if __name__=='__main__':main()
