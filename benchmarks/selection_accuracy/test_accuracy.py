"""Guard scientific failure modes: cohorts, budgets and dependent resampling."""
import copy
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

spec = importlib.util.spec_from_file_location("selection_accuracy_analysis", Path(__file__).with_name("analyze.py"))
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def cohort():
    jobs, rows = [], []
    for algorithm in ["bindash", "setsketch"]:
        for replicate in range(3):
            job = dict(job_id=f"{algorithm}-{replicate}", algorithm=algorithm, parameter=1,
                       regime_id="r", config_id=algorithm, k=24, n_a=100, n_b=100,
                       intersection=20, ground_truth=20 / 180, replicate=replicate,
                       data_seed=replicate, hash_seed=42, bbits=16, half_subk=6,
                       size_ratio=1, target_jaccard=20 / 180)
            jobs.append(job)
            rows.append(dict(job, status="ok", estimate=.1 + replicate * .01,
                             payload_bytes_a=4, payload_bytes_b=4, sketch_payload_bytes=8,
                             empty_sample_a=0, empty_sample_b=0,
                             input_digest_a=2**64 - 1 - replicate,
                             input_digest_b=2**63 + replicate, query_seconds=999))
    plan = dict(jobs=jobs, regimes=["r"], config=dict(replicates=3, analysis=dict(
        pair_budgets=[8], bootstrap_samples=100, bootstrap_seed=4,
        family_alpha=.05, minimum_relative_improvement=.1)))
    return plan, rows


def test_input_digests_remain_exact_and_timings_are_excluded():
    plan, rows = cohort()
    d = analysis.validate_rows(plan, rows, strict=True)
    assert int(d.input_digest_a.iloc[0]) == 2**64 - 1
    assert "query_seconds" not in d
    analysis.check_paired_inputs(d)


@pytest.mark.parametrize("corruption", ["missing", "duplicate", "digest", "truth", "payload", "negative_payload", "failure", "nonfinite"])
def test_corrupt_or_partial_run_is_rejected(corruption):
    plan, rows = cohort()
    if corruption == "missing":
        rows.pop()
    elif corruption == "duplicate":
        rows.append(copy.deepcopy(rows[0]))
    elif corruption == "digest":
        rows[-1]["input_digest_b"] += 1
    elif corruption == "truth":
        rows[0]["measured_ground_truth"] = .9
    elif corruption == "payload":
        rows[0]["sketch_payload_bytes"] = 9
    elif corruption == "negative_payload":
        rows[0].update(payload_bytes_a=-4, payload_bytes_b=12)
    elif corruption == "failure":
        rows[0]["status"] = "error"
    else:
        rows[0]["estimate"] = float("nan")
    with pytest.raises(ValueError):
        analysis.validate_rows(plan, rows, strict=True)


def test_budget_selection_uses_maximum_bytes_and_never_accuracy():
    data = pd.DataFrame([
        dict(parameter=1, config_id="oversized", mean_payload=90, max_payload=120, complete=True, mae=.001),
        dict(parameter=2, config_id="chosen", mean_payload=80, max_payload=80, complete=True, mae=.5),
        dict(parameter=3, config_id="less_payload", mean_payload=70, max_payload=70, complete=True, mae=0),
        dict(parameter=4, config_id="incomplete", mean_payload=95, max_payload=95, complete=False, mae=0),
    ]).assign(regime_id="r", algorithm="fracminhash")
    chosen = analysis.select_by_budget(data, [100])
    assert chosen.config_id.tolist() == ["chosen"]


def test_cross_run_pairing_is_validated():
    plan, rows = cohort()
    d = analysis.validate_rows(plan, rows, strict=True)
    d.loc[d.index[-1], "data_seed"] = 99
    with pytest.raises(ValueError, match="mismatch"):
        analysis.check_paired_inputs(d)


def test_pairing_across_budgets_preserves_dependence():
    differences = np.array([[.2, -.1, .3, .1, .4], [.2, -.1, .3, .1, .4]])
    bands, critical = analysis.paired_bands([differences], 2000, 7, .05)
    for values in bands[0].values():
        assert values[0] == values[1]
    _, single = analysis.paired_bands([differences[:1]], 2000, 7, .05)
    assert critical == single
    # Adding an identical, dependent budget cannot inflate the family.


def test_multiple_independent_regimes_do_not_narrow_family_band():
    x = np.array([[.2, -.1, .3, .1, .4]])
    one, c1 = analysis.paired_bands([x], 2000, 7, .05)
    many, c2 = analysis.paired_bands([x, x, x], 2000, 7, .05)
    assert c2 >= c1
    assert many[0]["simultaneous_low"][0] <= one[0]["simultaneous_low"][0]


@pytest.mark.parametrize("args,expected", [
    ((.1, .08, -.04, -.001, .01, .1), "setsketch"),
    ((.08, .1, .001, .04, .01, .1), "bindash"),
    ((.1, .08, -.04, .001, .01, .1), "unresolved"),
    ((.1, .095, -.01, -.001, .01, .1), "unresolved"),
    ((.1, .08, -.02, -.02, 0, .1), "unresolved"),
    ((0, 0, 0, 0, 0, .1), "both_exact"),
])
def test_preference_requires_evidence_and_effect_size(args, expected):
    assert analysis.preference(*args) == expected


def test_matched_comparison_rejects_inadequate_budget():
    plan, rows = cohort()
    d = analysis.validate_rows(plan, rows, strict=True)
    plan["config"]["analysis"]["pair_budgets"] = [7]
    with pytest.raises(ValueError, match="budget"):
        analysis.new_comparisons(d, plan)


def test_duplicate_replicate_cannot_be_counted_as_independent():
    plan, rows = cohort()
    d = analysis.validate_rows(plan, rows, strict=True)
    d = pd.concat([d, d.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        analysis.configuration_summary(d)
