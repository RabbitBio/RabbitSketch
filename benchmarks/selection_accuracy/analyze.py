#!/usr/bin/env python3
"""Accuracy-only conditional evidence; no timing fields or trained selector."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = Path("/home/bigssd/zt/rabbitsketch2_benchmark")
RESULTS = HERE.parent / "results/selection_accuracy_20261012"
NAMES = {"bindash": "BinDash", "setsketch": "SetSketch-fixed",
         "minhash": "FastKMV", "fracminhash": "FracMinHash", "hll": "HLL",
         "kssd": "KSSD", "probminhash": "ProbMinHash-unweighted",
         "bindash_hll": "BinDash+HLL"}
GROUP = ["regime_id", "config_id"]
NATIVE = ["estimate", "payload_bytes_a", "payload_bytes_b", "sketch_payload_bytes",
          "empty_sample_a", "empty_sample_b", "input_digest_a", "input_digest_b"]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def validate_rows(plan, rows, strict):
    """Check the complete registered input cohort before constructing a DataFrame."""
    expected = {r["job_id"]: r for r in plan["jobs"]}
    if len(expected) != len(plan["jobs"]):
        raise ValueError("duplicate registered job")
    seen, result = set(), []
    digests = {}
    for row in rows:
        job_id = row.get("job_id")
        if job_id not in expected or job_id in seen:
            raise ValueError("unknown or duplicate native job")
        seen.add(job_id)
        job = expected[job_id]
        if any(row.get(k) != v for k, v in job.items()):
            raise ValueError(f"registered metadata differs: {job_id}")
        valid = row.get("status") == "ok"
        if valid:
            value = row.get("estimate")
            valid = (value is not None and math.isfinite(value) and 0 <= value <= 1)
        if strict and not valid:
            raise ValueError(f"failed/nonfinite/invalid result: {job_id}")
        if valid:
            for key in ["payload_bytes_a", "payload_bytes_b", "sketch_payload_bytes"]:
                if not isinstance(row.get(key), int) or row[key] < 0:
                    raise ValueError("payload must be a nonnegative integer")
            truth = job["intersection"] / (job["n_a"] + job["n_b"] - job["intersection"])
            if not math.isclose(truth, row["ground_truth"], rel_tol=0, abs_tol=1e-15):
                raise ValueError("truth differs from exact set counts")
            measured = row.get("measured_ground_truth", row["ground_truth"])
            if not math.isclose(truth, measured, rel_tol=0, abs_tol=1e-15):
                raise ValueError("worker truth differs from plan")
            if row["payload_bytes_a"] + row["payload_bytes_b"] != row["sketch_payload_bytes"]:
                raise ValueError("inconsistent pair payload")
            key = (job["regime_id"], job["replicate"])
            signature = (row["input_digest_a"], row["input_digest_b"], job["data_seed"])
            if digests.setdefault(key, signature) != signature:
                raise ValueError("paired inputs differ between configurations")
        result.append(dict(job, **{k: row.get(k) for k in NATIVE}, valid=valid,
                           status=row.get("status"), expected=plan["config"]["replicates"]))
    if seen != set(expected):
        raise ValueError(f"incomplete run: missing {len(set(expected) - seen)} rows")
    return pd.DataFrame(result)


def read_run(directory, strict=False):
    plan = json.loads((directory / "plan.json").read_text())
    with (directory / "raw.jsonl").open() as stream:
        frame = validate_rows(plan, (json.loads(line) for line in stream), strict)
    return plan, frame


def check_paired_inputs(frame):
    # uint64 digests remain exact; do not coerce them through floating point.
    columns = ["input_digest_a", "input_digest_b", "data_seed", "ground_truth"]
    if (frame[frame.valid].groupby(["regime_id", "replicate"])[columns].nunique() > 1).any().any():
        raise ValueError("input/truth mismatch across frozen runs")


def configuration_summary(frame):
    d = frame.copy()
    d["error"] = np.where(d.valid, abs(d.estimate - d.ground_truth), np.nan)
    d["signed_error"] = np.where(d.valid, d.estimate - d.ground_truth, np.nan)
    d["squared_error"] = d.error ** 2
    d["empty"] = (d.empty_sample_a != 0) | (d.empty_sample_b != 0)
    d["zero"] = d.valid & (d.estimate == 0)
    meta = ["algorithm", "parameter", "k", "n_a", "n_b", "size_ratio",
            "target_jaccard", "ground_truth", "expected", "bbits", "half_subk"]
    if d.duplicated(GROUP + ["replicate"]).any():
        raise ValueError("duplicate configuration replicate")
    agg = d.groupby(GROUP).agg(
        observed=("replicate", "nunique"), valid=("valid", "sum"),
        empty_count=("empty", "sum"), mean_payload=("sketch_payload_bytes", "mean"),
        max_payload=("sketch_payload_bytes", "max"), mae=("error", "mean"),
        mse=("squared_error", "mean"), bias=("signed_error", "mean"),
        p95=("error", lambda x: x.quantile(.95)), zero_rate=("zero", "mean")).reset_index()
    agg = agg.merge(d[GROUP + meta].drop_duplicates(), on=GROUP, validate="one_to_one")
    agg["rmse"] = np.sqrt(agg.pop("mse"))
    agg["complete"] = (agg.observed == agg.expected) & (agg.valid == agg.expected)
    agg["all_nonempty"] = agg.empty_count == 0
    agg["relative_mae"] = agg.mae / agg.ground_truth.replace(0, np.nan)
    return agg


def select_by_budget(summary, budgets):
    """Largest tested query payload; never select using realized error."""
    rows = []
    for budget in budgets:
        eligible = summary[summary.complete & (summary.max_payload <= budget)]
        chosen = eligible.sort_values(["mean_payload", "parameter", "config_id"],
                                      ascending=[False, True, True]).drop_duplicates(
                                          ["regime_id", "algorithm"])
        rows.append(chosen.assign(budget=budget))
    return pd.concat(rows, ignore_index=True)


def paired_bands(groups, samples, seed, alpha):
    """groups: independent regimes, each a [budget, independent input] array.

    Budgets share resampling indices. The maximum standardized centered
    bootstrap deviation controls the complete reported family approximately;
    fixed empirical standard errors avoid unstable nested studentization.
    """
    rng = np.random.default_rng(seed)
    maximum = np.zeros(samples)
    output = []
    for values in groups:
        values = np.asarray(values, dtype=float)
        if values.ndim != 2 or values.shape[1] < 2 or not np.isfinite(values).all():
            raise ValueError("paired differences must be finite with independent replicates")
        n = values.shape[1]
        mean = values.mean(axis=1)
        se = values.std(axis=1, ddof=1) / np.sqrt(n)
        indices = rng.integers(0, n, size=(samples, n))
        means = values[:, indices].mean(axis=2)
        nondegenerate = se > 0
        if nondegenerate.any():
            z = abs((means[nondegenerate] - mean[nondegenerate, None]) /
                    se[nondegenerate, None])
            maximum = np.maximum(maximum, z.max(axis=0))
        lo, hi = np.quantile(means, [alpha / 2, 1 - alpha / 2], axis=1)
        output.append(dict(mean=mean, se=se, pointwise_low=lo, pointwise_high=hi))
    critical = float(np.quantile(maximum, 1 - alpha))
    for row in output:
        row["simultaneous_low"] = row["mean"] - critical * row["se"]
        row["simultaneous_high"] = row["mean"] + critical * row["se"]
    return output, critical


def preference(mae_bd, mae_ss, low, high, se, margin):
    if mae_bd == 0 and mae_ss == 0:
        return "both_exact"
    if se <= 0:
        return "unresolved"
    if high < 0 and mae_ss <= (1 - margin) * mae_bd:
        return "setsketch"
    if low > 0 and mae_bd <= (1 - margin) * mae_ss:
        return "bindash"
    return "unresolved"


def new_comparisons(frame, plan):
    config = plan["config"]
    policy = config["analysis"]
    summary = configuration_summary(frame)
    selected = select_by_budget(summary, policy["pair_budgets"])
    required = len(plan["regimes"]) * len(policy["pair_budgets"]) * 2
    if len(selected) != required or not selected.all_nonempty.all():
        raise ValueError("missing/empty/budget-infeasible matched configuration")
    if not (selected.max_payload == selected.budget).all():
        raise ValueError("new configurations do not exactly match the planned budget")
    errors = frame.assign(error=abs(frame.estimate - frame.ground_truth))
    matrices, records = [], []
    for regime_id, cohort in selected.groupby("regime_id", sort=True):
        differences, group_records = [], []
        for budget, cell in cohort.groupby("budget", sort=True):
            paired = {}
            for algorithm in ["bindash", "setsketch"]:
                row = cell[cell.algorithm == algorithm].iloc[0]
                f = errors[(errors.regime_id == regime_id) & (errors.config_id == row.config_id)]
                f = f.sort_values("replicate")
                if f.replicate.tolist() != list(range(config["replicates"])):
                    raise ValueError("nonmatching replicate cohort")
                paired[algorithm] = f.error.to_numpy()
            differences.append(paired["setsketch"] - paired["bindash"])
            base = cell.iloc[0]
            group_records.append(dict(regime_id=regime_id, k=int(base.k), n_b=int(base.n_b),
                                      size_ratio=int(base.size_ratio), target_jaccard=base.target_jaccard,
                                      ground_truth=base.ground_truth, budget=int(budget),
                                      containment_smaller=base.ground_truth * (1 + base.size_ratio) /
                                          (1 + base.ground_truth),
                                      replicates=config["replicates"],
                                      mae_bindash=paired["bindash"].mean(),
                                      mae_setsketch=paired["setsketch"].mean()))
        matrices.append(np.asarray(differences))
        records.append(group_records)
    bands, critical = paired_bands(matrices, policy["bootstrap_samples"],
                                  policy["bootstrap_seed"], policy["family_alpha"])
    rows = []
    for group_records, band in zip(records, bands):
        for i, row in enumerate(group_records):
            row.update({name: float(values[i]) for name, values in band.items()})
            row["ss_over_bd"] = (row["mae_setsketch"] / row["mae_bindash"]
                                 if row["mae_bindash"] > 0 else np.nan)
            row["preference"] = preference(row["mae_bindash"], row["mae_setsketch"],
                row["simultaneous_low"], row["simultaneous_high"], row["se"],
                policy["minimum_relative_improvement"])
            rows.append(row)
    return summary, pd.DataFrame(rows), critical


def frozen_screen(root, output, provenance):
    frames = []
    for name in ["pilot", "setsketch_validated_20261001", "frac_budget_extension_20261001"]:
        directory = root / "runs" / name
        _, data = read_run(directory)
        for filename in ["raw.jsonl", "plan.json", "manifest.json"]:
            provenance[str(directory / filename)] = sha(directory / filename)
        if name == "pilot":
            data = data[data.algorithm != "setsketch"]
        frames.append(data)
    data = pd.concat(frames, ignore_index=True)
    check_paired_inputs(data)
    summary = configuration_summary(data)
    summary.to_csv(output / "frozen_seven_configurations.csv", index=False)
    selected = select_by_budget(summary, [2048, 8192, 32768])
    # Availability is part of the decision, never convert an empty sketch to a win.
    usable = selected[selected.all_nonempty].copy()
    best = usable.groupby(["regime_id", "budget"]).mae.min().rename("best_mae")
    usable = usable.join(best, on=["regime_id", "budget"])
    usable["eligible_candidates"] = usable.groupby(["regime_id", "budget"]).algorithm.transform("nunique")
    usable["relative_to_best"] = usable.mae / usable.best_mae.replace(0, np.nan)
    usable["descriptive_best"] = np.isclose(usable.mae, usable.best_mae, rtol=0, atol=1e-15)
    usable.to_csv(output / "frozen_seven_budget_accuracy.csv", index=False)
    cells = summary.regime_id.nunique()
    coverage = selected.groupby(["budget", "algorithm"]).agg(
        budget_feasible_regimes=("regime_id", "nunique"),
        all_nonempty_regimes=("all_nonempty", "sum")).reset_index()
    coverage["total_regimes"] = cells
    wins = usable.groupby(["budget", "algorithm"]).descriptive_best.sum().rename("descriptive_best_cells")
    coverage = coverage.join(wins, on=["budget", "algorithm"])
    coverage.to_csv(output / "frozen_seven_coverage.csv", index=False)
    return dict(rows=len(data), regimes=cells, candidates=sorted(data.algorithm.unique()))


def read_evidence(path, provenance):
    provenance[str(path)] = sha(path)
    return pd.read_csv(path, low_memory=False)


def frozen_tasks(root, output, provenance):
    analysis = root / "analysis"
    for name in ["budget_accuracy", "paired_cluster_intervals", "budget_coverage"]:
        d = read_evidence(analysis / "robustness_20261001" / (name + ".csv"), provenance)
        if "metric" in d:
            d = d[d.metric == "j_error"]
        d.to_csv(output / ("real_" + name + ".csv"), index=False)
    for name in ["containment_budget_accuracy", "cardinality_accuracy", "weighted_accuracy"]:
        d = read_evidence(analysis / "regimes_20261001" / (name + ".csv"), provenance)
        if name.startswith("containment"):
            d = d[d.split == "validation"]
            d["public_metric_supported"] = d.algorithm != "kssd"
            d["all_outputs_available"] = d.availability == 1
            d["point_best_available_supported"] = False
            # No recommendation from the MAE of a selectively available subset.
            ok = d.public_metric_supported & d.all_outputs_available
            best = d[ok].groupby(["k", "budget", "nominal_ratio", "nominal_containment"]).mae.transform("min")
            d.loc[ok, "point_best_available_supported"] = np.isclose(d.loc[ok, "mae"], best, atol=1e-14, rtol=0)
        d.to_csv(output / ("heldout_" + name + ".csv" if name.startswith("containment") else name + ".csv"), index=False)

    # Directly examine retained query hashes, rather than equating ratio with support.
    selected = read_evidence(analysis / "regimes_20261001/containment_selected_rows.csv.gz", provenance)
    selected = selected[(selected.split == "validation") & (selected.algorithm == "fracminhash")]
    keys = ["pair_id", "k", "algorithm", "parameter"]
    choices = selected[keys + ["budget", "nominal_ratio", "nominal_containment"]].drop_duplicates()
    truth = read_evidence(analysis / "regimes_20261001/containment_exact_truth.csv", provenance)
    directory = root / "runs/containment_holdout_20261001"
    manifest_path = directory / "suite_manifest.json"
    provenance[str(manifest_path)] = sha(manifest_path)
    jobs = json.loads(manifest_path.read_text())["jobs"]
    rows = []
    for job in jobs:
        if job.get("algorithm") != "fracminhash":
            continue
        path = directory / "jobs" / job["job_id"] / "raw.jsonl"
        completion = json.loads((path.parent / "completion.json").read_text())
        actual = sha(path)
        if actual != completion["sha256"]:
            raise ValueError("frozen containment raw digest changed")
        provenance[str(path)] = actual
        for line in path.open():
            r = json.loads(line)
            if r["type"] == "pair":
                keep = keys + ["seed", "payload_bytes_a", "payload_bytes_b", "empty_a", "empty_b",
                               "containment_derived_ab", "containment_derived_ba"]
                rows.append({key: r[key] for key in keep})
    raw = pd.DataFrame(rows)
    if raw.duplicated(keys + ["seed"]).any():
        raise ValueError("duplicate frozen containment observation")
    d = raw.merge(choices, on=keys, validate="many_to_many").merge(
        truth[["pair_id", "k", "truth_ab"]], on=["pair_id", "k"], validate="many_to_one")
    counts = d.groupby(keys + ["budget"]).seed.nunique()
    if len(d) != len(selected) or not (counts == 3).all():
        raise ValueError("incomplete selected FracMinHash cohort")
    if (d.payload_bytes_a % 8 != 0).any() or (d.payload_bytes_a + d.payload_bytes_b > d.budget).any():
        raise ValueError("incorrect retained-hash payload")
    d["retained_query_hashes"] = d.payload_bytes_a // 8
    d["available"] = (d.empty_a == 0) & (d.empty_b == 0) & np.isfinite(d.containment_derived_ab)
    cap = np.maximum(1, np.maximum(d.containment_derived_ab, d.containment_derived_ba))
    d["error"] = np.where(d.available, abs(np.maximum(0, d.containment_derived_ab / cap) - d.truth_ab), np.nan)
    # Reproduce the frozen aggregate exactly before extending it with support counts.
    gkeys = ["k", "budget", "nominal_ratio", "nominal_containment"]
    old = selected.groupby(gkeys).agg(mae=("error", "mean"), availability=("available", "mean"))
    check = d.groupby(gkeys).agg(mae=("error", "mean"), availability=("available", "mean"))
    if not np.allclose(old.sort_index(), check.sort_index(), equal_nan=True, atol=1e-14, rtol=0):
        raise ValueError("retained-sample analysis does not reproduce frozen containment results")
    support = d.groupby(gkeys).agg(
        observations=("pair_id", "size"), pairs=("pair_id", "nunique"),
        mean_query_hashes=("retained_query_hashes", "mean"),
        minimum_query_hashes=("retained_query_hashes", "min"),
        maximum_query_hashes=("retained_query_hashes", "max"),
        availability=("available", "mean"), mae_available=("error", "mean"),
        p95_available=("error", lambda x: x.quantile(.95))).reset_index()
    support.to_csv(output / "fracminhash_query_support.csv", index=False)
    d.to_csv(output / "fracminhash_query_support_rows.csv.gz", index=False,
             compression={"method": "gzip", "mtime": 0})


def figures(output, comparisons):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    plt.rcParams.update({"pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
    states = ["bindash", "setsketch", "unresolved", "both_exact"]
    colors = ["#D55E00", "#009E73", "#D4D7DC", "#F2D16B"]
    cm = ListedColormap(colors)
    cm.set_bad("white")
    norm = BoundaryNorm(np.arange(-.5, 4.5), cm.N)
    js = sorted(comparisons.target_jaccard.unique())
    fig, axes = plt.subplots(2, 3, figsize=(15, 6), sharex=True, sharey=True)
    for row, cardinality in enumerate(sorted(comparisons.n_b.unique())):
        for col, budget in enumerate(sorted(comparisons.budget.unique())):
            ax = axes[row, col]
            sub = comparisons[(comparisons.n_b == cardinality) & (comparisons.budget == budget)]
            matrix = sub.pivot(index="size_ratio", columns="target_jaccard", values="preference")
            matrix = matrix.reindex(index=[1, 2, 4], columns=js)
            mapping = dict(zip(states, range(4)))
            matrix = matrix.apply(lambda column: column.map(mapping))
            ax.imshow(matrix.to_numpy(dtype=float), cmap=cm, norm=norm, aspect="auto")
            ax.set_title(f"Larger set = {cardinality:,}; pair payload = {budget // 1024} KiB", fontsize=10)
            ax.set_yticks(range(3), [1, 2, 4])
            ax.set_xticks(range(len(js)), [f"{j:g}" for j in js], rotation=60)
            ax.set_xticks(np.arange(-.5, len(js), 1), minor=True)
            ax.set_yticks(np.arange(-.5, 3, 1), minor=True)
            ax.grid(which="minor", color="white", linewidth=1)
            ax.tick_params(which="minor", bottom=False, left=False)
            if col == 0:
                ax.set_ylabel("Larger / smaller set")
            if row == 1:
                ax.set_xlabel("Exact Jaccard (explanatory variable)")
    labels = ["BinDash supported", "SetSketch-fixed supported", "Unresolved", "Both zero error", "Infeasible J > 1/ratio"]
    handles = [Patch(facecolor=c, label=l, edgecolor="#888888") for c, l in zip(colors + ["white"], labels)]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=9)
    fig.suptitle("Accuracy-only preference: paired simultaneous 95% bands + >=10% MAE reduction\nk=24; 30 independent inputs per cell; fixed hash seed", fontsize=12)
    fig.tight_layout(rect=[0, .065, 1, .92])
    for ext in ["png", "pdf"]:
        fig.savefig(output / ("overlap_selection_map." + ext), dpi=180)
    plt.close(fig)

    d = pd.read_csv(output / "fracminhash_query_support.csv")
    d = d[(d.k == 24) & (d.nominal_containment == .5)]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7))
    for budget, g in d.groupby("budget"):
        g = g.sort_values("nominal_ratio")
        label = f"{budget // 1024} KiB/pair"
        for ax, value in zip(axes, ["mean_query_hashes", "mae_available", "availability"]):
            ax.plot(g.nominal_ratio, g[value], "o-", label=label)
            ax.set_xscale("log")
            ax.set_xlabel("Reference / query length")
            ax.grid(alpha=.2)
    axes[0].set_yscale("log")
    for ax, label in zip(axes, ["Mean retained query hashes", "Containment MAE (available only)", "Available output fraction"]):
        ax.set_ylabel(label)
    axes[2].set_ylim(-.03, 1.03)
    axes[0].legend(fontsize=8)
    fig.suptitle("FracMinHash support: held-out references, k=24, nominal containment=0.5", fontsize=11)
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(output / ("containment_sample_support." + ext), dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--output", type=Path, default=RESULTS)
    args = parser.parse_args()
    run = args.run or args.root / "runs/selection_accuracy_20261012"
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    provenance = {}
    plan, data = read_run(run, strict=True)
    if plan["config"] != json.loads((HERE / "config.json").read_text()):
        raise ValueError("run differs from registered config")
    identity = json.loads((run / "manifest.json").read_text())["identity"]
    frozen_identity_path = args.root / "runs/robustness_20261001/suite_manifest.json"
    frozen_identity = json.loads(frozen_identity_path.read_text())
    if identity["binary_sha256"] != frozen_identity["binary_sha256"]:
        raise ValueError("new run did not use the frozen external-validation binary")
    source_path = args.root / "runs/setsketch_validated_20261001/manifest.json"
    corrected_source = json.loads(source_path.read_text())["identity"]["source"]
    if identity["source"] != corrected_source:
        raise ValueError("SetSketch correction/source differs from the validated screen")
    provenance[str(frozen_identity_path)] = sha(frozen_identity_path)
    provenance[str(source_path)] = sha(source_path)
    check_paired_inputs(data)
    summary, comparisons, critical = new_comparisons(data, plan)
    data.drop(columns=["query_min_ms"]).to_csv(output / "new_measurements.csv.gz", index=False,
                                               compression={"method": "gzip", "mtime": 0})
    summary.to_csv(output / "new_configuration_accuracy.csv", index=False)
    comparisons.to_csv(output / "new_paired_comparisons.csv", index=False)
    for name in ["plan.json", "raw.jsonl", "manifest.json"]:
        provenance[str(run / name)] = sha(run / name)
    (output / "native_manifest.json").write_bytes((run / "manifest.json").read_bytes())
    print(json.dumps(dict(new_rows=len(data), comparisons=len(comparisons),
                          preferences=comparisons.preference.value_counts().to_dict(),
                          simultaneous_critical=critical)), flush=True)
    frozen = frozen_screen(args.root, output, provenance)
    frozen_tasks(args.root, output, provenance)
    figures(output, comparisons)
    save_json(output / "analysis_manifest.json", dict(
        new_measurements=len(data), new_regimes=len(plan["regimes"]),
        infeasible_regimes=len(plan["skipped_regimes"]),
        paired_contrasts=len(comparisons), simultaneous_critical=critical,
        frozen_screen=frozen, timing_used=False, inputs_sha256=provenance,
        frozen_external_binary_matches=True, validated_corrected_source_matches=True,
        analysis_sha256=sha(__file__), config_sha256=sha(HERE / "config.json"),
        protocol_sha256=sha(HERE / "PROTOCOL.md"),
        numpy=np.__version__, pandas=pd.__version__))
    print(f"Accuracy evidence written to {output}", flush=True)


if __name__ == "__main__":
    main()
