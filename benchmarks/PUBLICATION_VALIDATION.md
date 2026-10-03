# Publication validation — 2026-10-03

This check packages the existing workflows and frozen evidence against source
revision `8db783ea75f597068b72705272ae207c94194e83`. It does not rerun the genome
experiments or change the recorded statistical conclusions.

## Checks performed

| Check | Result |
| --- | --- |
| Screening, main adapter, real-input and metric workflow tests | 62 passed |
| Independent-genome and mixed-mutation workflow tests | 20 passed |
| Recorded native partitioned SetSketch test binary | 450 model cases, 48 API cases and 20 sparse replicates passed |
| Native main adapter patch and SetSketch correction | Both pass `git apply --check` against the recorded baseline |
| Python and shell syntax | 47 Python files and 4 shell files passed |
| Repository-relative Markdown links in the publication | No missing tracked targets |
| Added artifacts | No compiled executables, object files or detected private keys/tokens |

The adapter tests used the recorded `build_validated/synthetic_pilot` and
`build_validated/rabbitsketch_pilot` binaries. The independent estimator check
used `build_validated/test_partitioned_setsketch`; its sparse MAE was
`0.000466766`. These are checks against the existing validated native build,
not a claim of a fresh full-library build or a second CPU evaluation.

## Advisor portability

The advisor can now read the bundled tables, an explicit `--data-dir`, or
`RABBITSKETCH_ADVISOR_DATA`. Its configuration selection logic is unchanged.
The bundled gzip payload decompresses to the original CSV byte for byte.

For `family=substitution`, `ANI=99–100`, `pair budget=8192 bytes` and
`target MAE=0.001`, all three data-source paths produced identical complete
JSON responses:

- Incident query reuse 1: MinHash, K=512, `validated_in_this_suite`.
- Incident query reuse 10,000: BinDash, bins=512, `validated_in_this_suite`.

The status refers to the original calibration/holdout suite. The independent
validation failures remain documented in the
[robustness report](../RabbitSketch2_robustness_validation.md); publication
does not turn those failures into transferable accuracy guarantees.

## Evidence integrity

[publication_manifest.json](results/publication_manifest.json) records the
published subset and omitted large tables. CSV and plot evidence is unchanged;
Markdown links are repository-relative. Whitespace checks exclude stored
`.patch` assets because a single space on a blank context line is part of the
unified diff format. Both source patches passed their applicability checks.
