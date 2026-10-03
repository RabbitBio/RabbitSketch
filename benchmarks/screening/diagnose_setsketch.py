#!/usr/bin/env python3
"""Reproduce the pilot's sparse-register anomaly using the same main and input.

Adds diagnostic output to a COPY of the frozen worker; native source and the
registered run are read-only. Reuses the native library and CLI objects from
the original build. Timings from this instrumented executable are not results.
"""
import datetime as dt
import json
import os
from pathlib import Path
import subprocess

import numpy as np

from run import DEFAULT_ROOT, HERE, atomic_json, file_sha
from run_main import write_manifest


def main():
    root = DEFAULT_ROOT
    run = root/'runs/pilot'
    source = root/'source/RabbitSketch-main'
    build = root/'build_main'
    output = HERE/'results/pilot_20261001/setsketch_diagnostic'
    output.mkdir(parents=True, exist_ok=True)
    worker = (run/'code_snapshot/synthetic_pilot.cpp').read_text()
    anchor = '                return r;\n            };\n            if (algorithm == "hll") {'
    extra = '''                if constexpr (std::is_same<T, Sketch::SetSketch>::value) {
                    const auto& a = *sketches[0]; const auto& b = *sketches[1];
                    const auto& ca = a.getCore(); const auto& cb = b.getCore();
                    size_t both_zero = 0, both_max = 0;
                    for (size_t i = 0; i < ca.size(); ++i) {
                        both_zero += ca[i] == 0 && cb[i] == 0;
                        both_max += ca[i] == a.getQ()+1 && cb[i] == a.getQ()+1;
                    }
                    out.put("diagnostic_both_zero_fraction", double(both_zero)/ca.size());
                    out.integer("diagnostic_range_fallback", (both_zero+both_max)*10 > ca.size()*3);
                    out.put("diagnostic_inclexcl", a.jaccard_index_inclexcl(b));
                    out.put("diagnostic_static_fallback", Sketch::SetSketch::jaccardFromCores(
                        ca.data(), cb.data(), a.getM(), a.getBaseInvPow(), a.getFactor(),
                        a.cardinality(), b.cardinality()));
                    out.put("diagnostic_scope", "estimates only; instrumented timings are invalid");
                }
'''
    if worker.count(anchor) != 1:
        raise ValueError('frozen worker changed; instrumentation anchor is ambiguous')
    instrumented = '#include <type_traits>\n' + worker.replace(anchor, extra+anchor)
    cpp = output/'instrumented_worker.cpp'
    cpp.write_text(instrumented)
    obj, binary = output/'instrumented_worker.o', output/'rabbitsketch_setsketch_diagnostic'
    compile_cmd = ['/usr/bin/c++','-O3','-DNDEBUG','-fopenmp','-std=gnu++17',
                   '-DRABBITSKETCH_SCREENING_EMBEDDED','-I'+str(source/'examples'),
                   '-I'+str(source/'src'),'-c',str(cpp),'-o',str(obj)]
    original_objects = [build/'CMakeFiles/rabbitsketch_pilot.dir'/
                        (str(source/'examples'/name).lstrip('/')+'.o') for name in ('main.cpp','sub_command.cpp')]
    original_objects.append(build/'CMakeFiles/rabbitsketch_pilot.dir/main_batch.cpp.o')
    library = build/'rabbitsketch/librabbitsketch_static.a'
    link_cmd = ['/usr/bin/c++','-O3','-fopenmp',*[str(p) for p in original_objects],str(obj),
                '-o',str(binary),str(library),'-Wl,-rpath,/usr/local/lib','-ldeflate','-lz','-lpthread']
    subprocess.run(compile_cmd, check=True)
    subprocess.run(link_cmd, check=True)
    plan = json.loads((run/'plan.json').read_text())
    records, summaries, commands = [], [], []
    env = dict(os.environ, OMP_DYNAMIC='FALSE', OMP_NUM_THREADS='64', OMP_THREAD_LIMIT='64')
    for precision in (12,14,16):
        jobs = [j for j in plan['jobs'] if j['algorithm']=='setsketch' and j['parameter']==precision and
                j['k']==20 and j['n_a']==10000 and j['n_b']==100000 and j['target_jaccard']==.001]
        assert len(jobs)==20
        manifest = output/f'p{precision}.tsv'
        write_manifest(manifest, jobs)
        command = [str(binary),'--screening-manifest',str(manifest),'-t','64']
        commands.append(command)
        proc = subprocess.run(command, check=True, text=True, capture_output=True, env=env, timeout=120)
        (output/f'p{precision}.stdout').write_text(proc.stdout)
        (output/f'p{precision}.stderr').write_text(proc.stderr)
        rows = [json.loads(s) for s in proc.stdout.splitlines() if s.startswith('{')]
        event = rows.pop(0)
        assert event['threads']==64 and len(rows)==20
        for row in rows:
            original = json.loads((run/'rows'/(row['job_id']+'.json')).read_text())
            assert all(row[key]==original[key] for key in ('estimate','ground_truth','input_digest_a','input_digest_b','sketch_payload_bytes'))
            row['diagnostic_precision'] = precision
        records.extend(rows)
        truth = rows[0]['ground_truth']
        summary = dict(precision=precision, k=20, n_a=10000, n_b=100000,
                       ground_truth=truth, replicates=len(rows),
                       default_estimate_mean=float(np.mean([r['estimate'] for r in rows])),
                       default_mae=float(np.mean([abs(r['estimate']-truth) for r in rows])),
                       inclexcl_mae=float(np.mean([abs(r['diagnostic_inclexcl']-truth) for r in rows])),
                       static_fallback_mae=float(np.mean([abs(r['diagnostic_static_fallback']-truth) for r in rows])),
                       both_zero_fraction_mean=float(np.mean([r['diagnostic_both_zero_fraction'] for r in rows])),
                       range_fallback_count=sum(r['diagnostic_range_fallback'] for r in rows),
                       exact_match_to_original=True)
        summaries.append(summary)
        print(json.dumps(summary),flush=True)
    with (output/'diagnostic_rows.jsonl').open('w') as stream:
        for row in records:
            stream.write(json.dumps(row,sort_keys=True)+'\n')
    atomic_json(output/'summary.json',summaries)
    atomic_json(output/'manifest.json',dict(
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), source_run=str(run),
        diagnostic_script_sha256=file_sha(__file__), instrumented_cpp_sha256=file_sha(cpp),
        native_cpp_sha256=file_sha(source/'src/SetSketch.cpp'),
        library_sha256=file_sha(library), original_objects_sha256={str(p):file_sha(p) for p in original_objects},
        binary_sha256=file_sha(binary), compile_command=compile_cmd, link_command=link_cmd,
        run_commands=commands, verified_default_estimates_and_inputs_match_original_pairs=60,
        timing_scope='not for performance; diagnostic queries run before normal timing'))


if __name__ == '__main__':
    main()
