#!/usr/bin/env python3
"""Check shipped main/test defaults against the measured binary and native oracle."""
from validate import REPO, RUN, OLD, NEW, PREVIOUS, environment, pairs, sha
import json
import subprocess


def main():
    directory = RUN/'default_validation'; directory.mkdir(exist_ok=True)
    rows = []
    panels = [('real128', PREVIOUS/'smoke.list'),
              ('ambiguous', PREVIOUS/'revisions/kssd_encoder_v2/validation/ambiguous.list'),
              ('formats', RUN/'validation/formats.list')]
    for panel, listing in panels:
        reference = None
        full_reference = None
        for entry, threads in [('measured',1), ('main_default',1), ('main_default',8),
                               ('main_explicit',1), ('test_default',1), ('test_default',8)]:
            tag = f'{panel}_{entry}_t{threads}'; output = directory/(tag+'.tsv')
            if entry == 'test_default':
                command = [REPO/'examples/exe_test_Kssd', listing, '.05', str(threads), output]
            else:
                binary = NEW if entry == 'measured' else REPO/'examples/rabbitsketch'
                command = [binary, '--kssd', *(['--index'] if entry in ['measured','main_explicit'] else []),
                           '-i', listing, '-l', '-t', str(threads), '-o', output]
            process = subprocess.run(list(map(str,command)), cwd=REPO/'examples',
                                     env=environment(threads), capture_output=True, text=True, check=True)
            (directory/(tag+'.log')).write_text(process.stdout+process.stderr)
            assert 'build CSR index:' in process.stderr and 'sketch time:' in process.stderr
            observed = pairs(output)
            full = sorted(output.read_text().splitlines())
            if reference is None: reference, full_reference = observed, full
            assert observed == reference and full == full_reference, tag
            rows.append(dict(case=tag, match=True, pairs=len(observed), command=list(map(str,command))))
        print(panel, 'main/default/test/measured MATCH', len(reference), flush=True)
    for k in [22,32]:
        reference = None
        for entry, binary, indexed in [('native_before',OLD,False),
                                       ('main_default',REPO/'examples/rabbitsketch',True),
                                       ('native_after',REPO/'examples/rabbitsketch',False)]:
            tag=f'k{k}_{entry}'; output=directory/(tag+'.tsv')
            command=[binary,'--kssd','-k',str(k),'-i',RUN/'validation/formats.list','-l','-t','2','-o',output]
            process=subprocess.run(list(map(str,command)),cwd=REPO/'examples',
                env=environment(2,indexed),capture_output=True,text=True,check=True)
            (directory/(tag+'.log')).write_text(process.stdout+process.stderr)
            observed=pairs(output)
            if reference is None: reference=observed
            assert observed==reference,tag
            rows.append(dict(case=tag,match=True,pairs=len(observed),command=list(map(str,command))))
        print('k',k,'default index / native MATCH',flush=True)
    result=dict(all_passed=True,cases=rows,
        main_sha256=sha(REPO/'examples/rabbitsketch'),test_sha256=sha(REPO/'examples/exe_test_Kssd'),
        measured_binary_sha256=sha(NEW))
    (RUN/'default_validation.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__': main()
