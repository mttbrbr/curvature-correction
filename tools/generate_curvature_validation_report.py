#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Validate curvature V&V data and generate Markdown, HTML and JSON reports.

The input is a versioned JSON document.  Profile/scalar datasets may contain
inline arrays or point to a CSV file, making the same tool suitable for small
committed snapshots and freshly generated OpenFOAM post-processing results.
Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from html import escape
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable

SCHEMA_VERSION = 2
EPS = 1.0e-30
VALID_CLASSIFICATIONS = {"improves", "neutral", "worsens", "non-conclusive"}
VALID_OPERATORS = {"le", "ge"}


class DataError(ValueError):
    """Raised when validation input is incomplete or malformed."""


@dataclass(frozen=True)
class Metrics:
    count: int
    rmse: float
    normalized_rmse: float
    mae: float
    mean_bias: float
    maximum_absolute_error: float
    mean_relative_error: float | None
    maximum_relative_error: float | None


def finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise DataError(f"{label} must be a number")
    result = float(value)
    if not math.isfinite(result):
        raise DataError(f"{label} contains NaN or infinity")
    return result


def evaluate_criterion(criterion: dict[str, Any]) -> str:
    """Compute a criterion status from its declared observation and threshold."""
    if "status" in criterion:
        raise DataError("acceptance criterion status is computed and must not be supplied")
    if bool(criterion.get("not_applicable", False)):
        return "not-applicable"
    if criterion.get("observed") is None:
        return "not-verified"
    operator = str(required(criterion, "operator", "acceptance criterion"))
    if operator not in VALID_OPERATORS:
        raise DataError(f"invalid acceptance operator {operator}")
    observed = finite_number(criterion["observed"], "acceptance criterion.observed")
    threshold = finite_number(
        required(criterion, "threshold", "acceptance criterion"),
        "acceptance criterion.threshold",
    )
    passed = observed <= threshold if operator == "le" else observed >= threshold
    return "pass" if passed else "fail"


def weighted_metrics(model: Iterable[float], reference: Iterable[float], weights: Iterable[float]) -> Metrics:
    values = [float(item) for item in model]
    refs = [float(item) for item in reference]
    wts = [float(item) for item in weights]
    if not values or len(values) != len(refs) or len(values) != len(wts):
        raise DataError("model, reference and weight arrays must have equal non-zero length")
    if not all(math.isfinite(item) for item in values + refs + wts):
        raise DataError("metric arrays contain NaN or infinity")
    if any(item <= 0 for item in wts):
        raise DataError("all metric weights must be positive")

    total_weight = sum(wts)
    errors = [value - ref for value, ref in zip(values, refs)]
    rmse = math.sqrt(sum(weight * error * error for weight, error in zip(wts, errors)) / total_weight)
    reference_rms = math.sqrt(sum(weight * ref * ref for weight, ref in zip(wts, refs)) / total_weight)
    relative = [error / ref for error, ref in zip(errors, refs) if abs(ref) > EPS]
    return Metrics(
        count=len(values),
        rmse=rmse,
        normalized_rmse=rmse / max(reference_rms, EPS),
        mae=sum(weight * abs(error) for weight, error in zip(wts, errors)) / total_weight,
        mean_bias=sum(weight * error for weight, error in zip(wts, errors)) / total_weight,
        maximum_absolute_error=max(abs(error) for error in errors),
        mean_relative_error=(sum(relative) / len(relative)) if relative else None,
        maximum_relative_error=max((abs(item) for item in relative), default=None),
    )


def required(mapping: dict[str, Any], key: str, context: str) -> Any:
    if key not in mapping:
        raise DataError(f"missing {context}.{key}")
    return mapping[key]


def load_series(dataset: dict[str, Any], input_dir: Path) -> dict[str, list[float]]:
    columns = ("x", "baseline", "corrected", "reference", "weight")
    if "csv" in dataset:
        path = (input_dir / str(dataset["csv"])).resolve()
        try:
            with path.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
        except OSError as error:
            raise DataError(f"cannot read dataset CSV {path}: {error}") from error
        if not rows:
            raise DataError(f"dataset CSV {path} is empty")
        result: dict[str, list[float]] = {}
        for column in columns:
            if column not in rows[0] or all(row.get(column, "") == "" for row in rows):
                continue
            result[column] = [finite_number(float(row[column]), f"{path}:{column}") for row in rows]
    else:
        result = {}
        for column in columns:
            if column in dataset:
                if not isinstance(dataset[column], list):
                    raise DataError(f"dataset {dataset.get('id', '?')} {column} must be an array")
                result[column] = [finite_number(item, f"dataset {dataset.get('id', '?')}:{column}") for item in dataset[column]]

    for key in ("baseline", "corrected"):
        if key not in result:
            raise DataError(f"dataset {dataset.get('id', '?')} lacks {key} values")
    size = len(result["baseline"])
    if size == 0 or any(len(items) != size for items in result.values()):
        raise DataError(f"dataset {dataset.get('id', '?')} columns have inconsistent lengths")
    result.setdefault("x", [float(index) for index in range(size)])
    result.setdefault("weight", [1.0] * size)
    if any(weight <= 0 for weight in result["weight"]):
        raise DataError(f"dataset {dataset.get('id', '?')} has non-positive weights")
    return result


def classify(
    baseline: Metrics | None,
    corrected: Metrics | None,
    *,
    converged: bool,
    grid_independent: bool | None,
    neutral_tolerance: float,
) -> str:
    if baseline is None or corrected is None or not converged or grid_independent is False:
        return "non-conclusive"
    scale = max(baseline.rmse, corrected.rmse, EPS)
    change = (baseline.rmse - corrected.rmse) / scale
    if change > neutral_tolerance:
        return "improves"
    if change < -neutral_tolerance:
        return "worsens"
    return "neutral"


def classify_neutral(
    baseline: list[float],
    corrected: list[float],
    *,
    converged: bool,
    grid_independent: bool | None,
    neutral_tolerance: float,
) -> str:
    if not converged or grid_independent is False:
        return "non-conclusive"
    maximum_change = max(
        abs(after - before) / max(abs(before), EPS)
        for before, after in zip(baseline, corrected)
    )
    return "neutral" if maximum_change <= neutral_tolerance else "worsens"


def svg_plot(path: Path, title: str, x: list[float], series: list[tuple[str, list[float], str]]) -> None:
    width, height = 920, 520
    left, right, top, bottom = 78, 24, 48, 62
    all_y = [value for _, values, _ in series for value in values]
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(all_y), max(all_y)
    if xmax == xmin:
        xmin, xmax = xmin - 0.5, xmax + 0.5
    if ymax == ymin:
        ymin, ymax = ymin - 0.5, ymax + 0.5
    pad = 0.05 * (ymax - ymin)
    ymin, ymax = ymin - pad, ymax + pad

    def sx(value: float) -> float:
        return left + (value - xmin) * (width - left - right) / (xmax - xmin)

    def sy(value: float) -> float:
        return height - bottom - (value - ymin) * (height - top - bottom) / (ymax - ymin)

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="17">{escape(title)}</text>',
    ]
    for tick in range(6):
        value = ymin + (ymax - ymin) * tick / 5
        y_value = sy(value)
        lines.append(f'<line x1="{left}" y1="{y_value:.2f}" x2="{width-right}" y2="{y_value:.2f}" stroke="#e5e5e5"/>')
        lines.append(f'<text x="{left-8}" y="{y_value+4:.2f}" text-anchor="end" font-family="sans-serif" font-size="10">{value:.4g}</text>')
    lines.append(f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#222"/>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#222"/>')
    for index, (name, values, colour) in enumerate(series):
        points = " ".join(f"{sx(xv):.2f},{sy(yv):.2f}" for xv, yv in zip(x, values))
        lines.append(f'<polyline points="{points}" fill="none" stroke="{colour}" stroke-width="2"/>')
        legend_x = 100 + index * 230
        lines.append(f'<line x1="{legend_x}" y1="{height-22}" x2="{legend_x+24}" y2="{height-22}" stroke="{colour}" stroke-width="3"/>')
        lines.append(f'<text x="{legend_x+31}" y="{height-18}" font-family="sans-serif" font-size="11">{escape(name)}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n")


def analyse(data: dict[str, Any], input_path: Path, output_dir: Path) -> dict[str, Any]:
    if required(data, "schema_version", "root") != SCHEMA_VERSION:
        raise DataError(f"schema_version must be {SCHEMA_VERSION}")
    metadata = required(data, "metadata", "root")
    implementation = required(data, "implementation", "root")
    benchmarks = required(data, "benchmarks", "root")
    criteria = required(data, "acceptance_criteria", "root")
    if not isinstance(metadata, dict) or not isinstance(implementation, dict):
        raise DataError("metadata and implementation must be objects")
    if not isinstance(benchmarks, list) or not isinstance(criteria, list):
        raise DataError("benchmarks and acceptance_criteria must be arrays")
    for key in (
        "campaign_id",
        "commit",
        "openfoam_version",
        "compiler",
        "archive_manifest_sha256",
    ):
        required(metadata, key, "metadata")

    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    neutral_tolerance = finite_number(metadata.get("neutral_tolerance", 0.02), "metadata.neutral_tolerance")
    analysed_benchmarks: list[dict[str, Any]] = []
    missing_reference: list[str] = []
    non_converged: list[str] = []

    for benchmark in benchmarks:
        if not isinstance(benchmark, dict):
            raise DataError("each benchmark must be an object")
        benchmark_id = str(required(benchmark, "id", "benchmark"))
        required(benchmark, "name", f"benchmark {benchmark_id}")
        required(benchmark, "category", f"benchmark {benchmark_id}")
        required(benchmark, "case_revision", f"benchmark {benchmark_id}")
        digest = str(required(benchmark, "input_sha256", f"benchmark {benchmark_id}"))
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise DataError(f"benchmark {benchmark_id}.input_sha256 must be SHA-256")
        meshes = required(benchmark, "meshes", f"benchmark {benchmark_id}")
        if not isinstance(meshes, list) or len(meshes) < 3:
            raise DataError(f"benchmark {benchmark_id}.meshes must contain at least three meshes")
        converged = bool(required(benchmark, "converged", f"benchmark {benchmark_id}"))
        if not converged:
            non_converged.append(benchmark_id)
        datasets = benchmark.get("datasets", [])
        if not isinstance(datasets, list):
            raise DataError(f"benchmark {benchmark_id}.datasets must be an array")
        analysed_sets = []
        for dataset in datasets:
            if not isinstance(dataset, dict):
                raise DataError(f"benchmark {benchmark_id} dataset must be an object")
            dataset_id = str(required(dataset, "id", f"benchmark {benchmark_id} dataset"))
            series = load_series(dataset, input_path.parent)
            has_reference = "reference" in series
            if not has_reference:
                missing_reference.append(f"{benchmark_id}/{dataset_id}")
            publication = str(dataset.get("reference_publication", "aggregate-only"))
            if publication not in {"open", "aggregate-only"}:
                raise DataError(
                    f"dataset {dataset_id}.reference_publication must be open or aggregate-only"
                )
            if has_reference:
                required(dataset, "reference_citation", f"dataset {dataset_id}")
            base_metrics = weighted_metrics(series["baseline"], series["reference"], series["weight"]) if has_reference else None
            corrected_metrics = weighted_metrics(series["corrected"], series["reference"], series["weight"]) if has_reference else None
            if dataset.get("expected_behavior") == "neutral":
                classification = classify_neutral(
                    series["baseline"],
                    series["corrected"],
                    converged=converged,
                    grid_independent=benchmark.get("grid_independent"),
                    neutral_tolerance=neutral_tolerance,
                )
            else:
                classification = classify(
                    base_metrics,
                    corrected_metrics,
                    converged=converged,
                    grid_independent=benchmark.get("grid_independent"),
                    neutral_tolerance=neutral_tolerance,
                )
            if classification not in VALID_CLASSIFICATIONS:
                raise AssertionError(classification)
            plot_name = f"{benchmark_id}-{dataset_id}.svg".replace("/", "-")
            plot_series = [
                (str(dataset.get("baseline_label", "baseline")), series["baseline"], "#2864b7"),
                (str(dataset.get("corrected_label", "corrected")), series["corrected"], "#c53d36"),
            ]
            if has_reference and publication == "open":
                plot_series.append((str(dataset.get("reference_label", "reference")), series["reference"], "#222222"))
            svg_plot(plots_dir / plot_name, str(dataset.get("quantity", dataset_id)), series["x"], plot_series)
            analysed_sets.append(
                {
                    "id": dataset_id,
                    "quantity": str(dataset.get("quantity", dataset_id)),
                    "classification": classification,
                    "baseline_metrics": asdict(base_metrics) if base_metrics else None,
                    "corrected_metrics": asdict(corrected_metrics) if corrected_metrics else None,
                    "baseline_to_corrected_percent": [
                        100.0 * (corrected - baseline) / max(abs(baseline), EPS)
                        for baseline, corrected in zip(series["baseline"], series["corrected"])
                    ],
                    "plot": f"plots/{plot_name}",
                    "weighting": str(dataset.get("weighting", "uniform point weights")),
                    "reference_citation": dataset.get("reference_citation"),
                    "reference_publication": publication,
                }
            )
        item = dict(benchmark)
        item["datasets"] = analysed_sets
        analysed_benchmarks.append(item)

    critical_failures = []
    critical_incomplete = []
    analysed_criteria = []
    for criterion in criteria:
        if not isinstance(criterion, dict):
            raise DataError("each acceptance criterion must be an object")
        status = evaluate_criterion(criterion)
        analysed = dict(criterion)
        analysed["status"] = status
        analysed_criteria.append(analysed)
        if bool(criterion.get("critical", False)):
            identifier = str(
                criterion.get("id", criterion.get("criterion", "unknown"))
            )
            if status == "fail":
                critical_failures.append(identifier)
            if status != "pass":
                critical_incomplete.append(identifier)

    return {
        "schema_version": SCHEMA_VERSION,
        "metadata": metadata,
        "implementation": implementation,
        "benchmarks": analysed_benchmarks,
        "acceptance_criteria": analysed_criteria,
        "missing_benchmarks": data.get("missing_benchmarks", []),
        "merge_blockers": data.get("merge_blockers", []),
        "validation_status": data.get("validation_status", []),
        "merge_conditions": data.get("merge_conditions", []),
        "known_issues": data.get("known_issues", []),
        "commands": data.get("commands", []),
        "non_converged": non_converged,
        "missing_reference": missing_reference,
        "critical_failures": critical_failures,
        "critical_incomplete": critical_incomplete,
    }


def number(value: Any) -> str:
    if value is None:
        return "n/a"
    return f"{value:.6g}" if isinstance(value, float) else str(value)


def markdown_report(result: dict[str, Any]) -> str:
    meta = result["metadata"]
    impl = result["implementation"]
    lines = [
        f"# {meta.get('title', 'Curvature-correction validation report')}",
        "",
        f"Generated from schema v{result['schema_version']}. Commit: `{meta.get('commit', 'not recorded')}`.",
        "",
        "## Executive summary",
        "",
        str(meta.get("executive_summary", "Not provided.")),
        "",
        f"**Final recommendation: {meta.get('recommendation', 'not validated')}**",
        "",
        "## Implementation audit",
        "",
        str(impl.get("summary", "Not provided.")),
        "",
        "### Effective formulas",
        "",
        str(impl.get("formulas", "Not provided.")),
        "",
        "### Constants",
        "",
        "| Model | Constant | Default | Units | Role |",
        "| --- | --- | ---: | --- | --- |",
    ]
    for item in impl.get("constants", []):
        lines.append(f"| {item['model']} | `{item['name']}` | {item['default']} | {item.get('units', 'dimensionless')} | {item['role']} |")
    lines.extend(["", "### Code map", ""])
    for item in impl.get("code_map", []):
        lines.append(f"- `{item['path']}` — {item['role']}")
    lines.extend(["", "### Numerical and modelling risks", ""])
    for item in impl.get("risks", []):
        recommendation = str(item.get("recommendation", "not provided")).rstrip(".")
        lines.append(f"- **{item.get('severity', 'unspecified')}** — {item['description']} Recommendation: {recommendation}.")
    lines.extend(["", "### Literature comparison", "", str(impl.get("literature_comparison", "Not provided.")), ""])
    lines.extend(["## Merge blockers", ""])
    for item in result["merge_blockers"]:
        lines.append(f"- {item}")
    lines.extend(["", "## Current validation status", ""])
    for item in result["validation_status"]:
        lines.append(f"- {item}")

    lines.extend([
        "## Verification",
        "",
        str(meta.get("verification_summary", "See acceptance criteria and benchmark matrix.")),
        "",
        "## Validation matrix",
        "",
        "| Benchmark | Category | Models | Grid/cells | Reference | Converged | Grid-independent | Result |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ])
    for benchmark in result["benchmarks"]:
        classifications = sorted({item["classification"] for item in benchmark["datasets"]}) or ["non-conclusive"]
        lines.append(
            f"| {benchmark['name']} | {benchmark['category']} | {benchmark.get('models', 'n/a')} | "
            f"{benchmark.get('grid', 'n/a')} / {benchmark.get('cells', 'n/a')} | {benchmark.get('reference', 'none')} | "
            f"{benchmark['converged']} | {benchmark.get('grid_independent', 'not verified')} | {', '.join(classifications)} |"
        )

    lines.extend(["", "## Quantitative results", ""])
    for benchmark in result["benchmarks"]:
        lines.extend([f"### {benchmark['name']}", "", str(benchmark.get("discussion", "No discussion supplied.")), ""])
        for dataset in benchmark["datasets"]:
            lines.extend([f"#### {dataset['quantity']} — {dataset['classification']}", "", f"![{dataset['quantity']}]({dataset['plot']})", ""])
            lines.append("| Variant | RMSE | nRMSE | MAE | Bias | Max abs | Mean rel. | Max rel. |")
            lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
            for label, metrics in (("baseline", dataset["baseline_metrics"]), ("corrected", dataset["corrected_metrics"])):
                if metrics is None:
                    lines.append(f"| {label} | n/a | n/a | n/a | n/a | n/a | n/a | n/a |")
                else:
                    lines.append(
                        f"| {label} | {number(metrics['rmse'])} | {number(metrics['normalized_rmse'])} | "
                        f"{number(metrics['mae'])} | {number(metrics['mean_bias'])} | {number(metrics['maximum_absolute_error'])} | "
                        f"{number(metrics['mean_relative_error'])} | {number(metrics['maximum_relative_error'])} |"
                    )
            lines.extend(["", f"Weights: {dataset['weighting']}.", ""])

    lines.extend(["## Grid convergence", "", str(meta.get("grid_convergence_summary", "Not verified.")), ""])
    lines.extend(["## Sensitivity and robustness", "", str(meta.get("sensitivity_summary", "Not verified.")), ""])
    lines.extend(["## Acceptance criteria", "", "| Criterion | Critical | Observed | Rule | Status | Evidence |", "| --- | --- | ---: | --- | --- | --- |"])
    for criterion in result["acceptance_criteria"]:
        rule = "n/a" if criterion.get("not_applicable") else f"{criterion.get('operator', '?')} {criterion.get('threshold', '?')} {criterion.get('unit', '')}".strip()
        lines.append(f"| {criterion['criterion']} | {criterion.get('critical', False)} | {criterion.get('observed', 'n/a')} | {rule} | {criterion['status']} | {criterion.get('evidence', '')} |")
    lines.extend(["", "## Missing or non-conclusive work", ""])
    for item in result["missing_benchmarks"]:
        lines.append(f"- {item}")
    if result["non_converged"]:
        lines.append(f"- Non-converged cases: {', '.join(result['non_converged'])}")
    if result["missing_reference"]:
        lines.append(f"- Datasets without reference: {', '.join(result['missing_reference'])}")
    lines.extend(["", "## Known issues and limitations", ""])
    for item in result["known_issues"]:
        lines.append(f"- {item}")
    lines.extend(["", "## Conditions before definitive merge", ""])
    for index, item in enumerate(result["merge_conditions"], start=1):
        lines.append(f"{index}. {item}")
    lines.extend(["", "## Reproduction", ""])
    for command in result["commands"]:
        lines.extend(["```bash", command, "```", ""])
    lines.extend(["## Conclusion", "", str(meta.get("conclusion", meta.get("recommendation", "not validated"))), ""])
    return "\n".join(lines)


def html_report(markdown: str, title: str) -> str:
    """Render a self-contained, deliberately small Markdown subset."""
    rendered = []
    in_code = False
    for line in markdown.splitlines():
        if line.startswith("```"):
            rendered.append("</code></pre>" if in_code else "<pre><code>")
            in_code = not in_code
        elif in_code:
            rendered.append(escape(line))
        elif line.startswith("#### "):
            rendered.append(f"<h4>{escape(line[5:])}</h4>")
        elif line.startswith("### "):
            rendered.append(f"<h3>{escape(line[4:])}</h3>")
        elif line.startswith("## "):
            rendered.append(f"<h2>{escape(line[3:])}</h2>")
        elif line.startswith("# "):
            rendered.append(f"<h1>{escape(line[2:])}</h1>")
        elif line.startswith("![") and "](" in line and line.endswith(")"):
            alt, target = line[2:].split("](", 1)
            rendered.append(f'<img src="{escape(target[:-1])}" alt="{escape(alt)}">')
        elif line.startswith("- "):
            rendered.append(f"<p>• {escape(line[2:])}</p>")
        elif line.startswith("|"):
            rendered.append(f"<pre>{escape(line)}</pre>")
        elif line:
            rendered.append(f"<p>{escape(line)}</p>")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>{escape(title)}</title>
<style>body{{font:15px/1.5 system-ui,sans-serif;max-width:1120px;margin:2rem auto;padding:0 1rem;color:#222}}h1,h2,h3{{color:#183b63}}pre{{background:#f5f7fa;padding:.55rem;overflow:auto}}img{{max-width:100%;border:1px solid #ddd}}code{{font-family:ui-monospace,monospace}}</style></head>
<body>{''.join(rendered)}</body></html>\n"""


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="versioned validation JSON")
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="campaign-specific output directory",
    )
    parser.add_argument(
        "--strict-acceptance",
        action="store_true",
        help="return 1 unless every critical acceptance criterion passes",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        data = json.loads(args.input.read_text())
        if not isinstance(data, dict):
            raise DataError("input root must be an object")
        result = analyse(data, args.input.resolve(), args.output_dir)
        markdown = markdown_report(result)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "curvature-validation-report.md").write_text(markdown + "\n")
        (args.output_dir / "curvature-validation-report.html").write_text(
            html_report(markdown, str(result["metadata"].get("title", "Curvature validation")))
        )
        (args.output_dir / "aggregated-metrics.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    except (OSError, json.JSONDecodeError, DataError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    print(f"Markdown: {args.output_dir / 'curvature-validation-report.md'}")
    print(f"HTML: {args.output_dir / 'curvature-validation-report.html'}")
    print(f"Metrics: {args.output_dir / 'aggregated-metrics.json'}")
    blocking = (
        result["critical_incomplete"]
        if args.strict_acceptance
        else result["critical_failures"]
    )
    if blocking:
        label = "incomplete" if args.strict_acceptance else "failures"
        print(
            f"ERROR: critical validation {label}: " + ", ".join(blocking),
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
