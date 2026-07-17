#!/usr/bin/env python3
"""Analyse skin-friction sensitivity across the nested TMBWG grids."""

from __future__ import annotations

import bisect
import csv
import math
import sys
from pathlib import Path

from analyse import (
    extract_cf,
    near_zero_reduction,
    read_experiment,
    read_reference,
    write_svg,
)

LEVELS = ("level3-129x49", "level2-257x97", "level1-513x193")
MODELS = ("sst", "sstcc-on", "sa", "sacc-on")


def interpolated_value(curve: list[tuple[float, float]], x_value: float) -> float:
    coordinates = [x for x, _ in curve]
    index = bisect.bisect_left(coordinates, x_value)
    if index == 0:
        return curve[0][1]
    if index == len(curve):
        return curve[-1][1]
    x0, y0 = curve[index - 1]
    x1, y1 = curve[index]
    fraction = (x_value - x0) / (x1 - x0)
    return y0 + fraction*(y1 - y0)


def relative_rms(
    sampled: list[tuple[float, float]],
    reference: list[tuple[float, float]],
) -> float:
    pairs = [
        (value, interpolated_value(reference, x_value))
        for x_value, value in sampled
    ]
    numerator = sum((first - second)**2 for first, second in pairs)
    denominator = max(sum(second**2 for _, second in pairs), 1e-30)
    return math.sqrt(numerator/denominator)


def main() -> int:
    tutorial = Path(sys.argv[1]).resolve()
    results = tutorial / "runs" / "grid-study"
    curves: dict[str, dict[str, list[tuple[float, float]]]] = {}

    for level in LEVELS:
        curves[level] = {}
        for model in MODELS:
            curve, _ = extract_cf(results / level / model)
            curves[level][model] = curve
            with (results / level / model / "cf.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("x_m", "Cf_Uref"))
                writer.writerows(curve)

    metrics: dict[str, float] = {}
    for level in LEVELS:
        metrics[f"{level}_sst_reduction"] = near_zero_reduction(
            curves[level]["sst"], curves[level]["sstcc-on"]
        )
        metrics[f"{level}_sa_reduction"] = near_zero_reduction(
            curves[level]["sa"], curves[level]["sacc-on"]
        )

    for model in MODELS:
        coarse_medium = relative_rms(
            curves[LEVELS[0]][model], curves[LEVELS[1]][model]
        )
        medium_fine = relative_rms(
            curves[LEVELS[1]][model], curves[LEVELS[2]][model]
        )
        metrics[f"{model}_129_to_257_relative_rms"] = coarse_medium
        metrics[f"{model}_257_to_513_relative_rms"] = medium_fine
        metrics[f"{model}_convergence_ratio"] = medium_fine/max(coarse_medium, 1e-30)

    with (results / "summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("metric", "value"))
        writer.writerows((name, f"{value:.12g}") for name, value in metrics.items())

    reference = tutorial / "resources" / "reference"
    experiment = read_experiment(
        reference / "Smits_curvature_expdata_reworked2012.dat"
    )
    write_svg(
        results / "sst-grid-convergence.svg",
        [
            ("SST 129x49", curves[LEVELS[0]]["sst"], "#9ecae1", False),
            ("SST 257x97", curves[LEVELS[1]]["sst"], "#4292c6", False),
            ("SST 513x193", curves[LEVELS[2]]["sst"], "#08519c", False),
            ("SST-RC 129x49", curves[LEVELS[0]]["sstcc-on"], "#fcae91", False),
            ("SST-RC 257x97", curves[LEVELS[1]]["sstcc-on"], "#fb6a4a", False),
            ("SST-RC 513x193", curves[LEVELS[2]]["sstcc-on"], "#cb181d", False),
            (
                "TMBWG SST-RC",
                read_reference(reference / "smits_cf_cfl3d_sstrc.dat"),
                "#7f0000",
                False,
            ),
            ("Experiment", experiment, "#111111", True),
        ],
    )
    write_svg(
        results / "sa-grid-convergence.svg",
        [
            ("SA 129x49", curves[LEVELS[0]]["sa"], "#a1d99b", False),
            ("SA 257x97", curves[LEVELS[1]]["sa"], "#41ab5d", False),
            ("SA 513x193", curves[LEVELS[2]]["sa"], "#006d2c", False),
            ("SA-RC 129x49", curves[LEVELS[0]]["sacc-on"], "#dadaeb", False),
            ("SA-RC 257x97", curves[LEVELS[1]]["sacc-on"], "#9e9ac8", False),
            ("SA-RC 513x193", curves[LEVELS[2]]["sacc-on"], "#6a51a3", False),
            (
                "TMBWG SA-RC",
                read_reference(reference / "smits_cf_cfl3d_sarc.dat"),
                "#3f007d",
                False,
            ),
            ("Experiment", experiment, "#111111", True),
        ],
    )

    failures = []
    for level in LEVELS:
        for family in ("sst", "sa"):
            value = metrics[f"{level}_{family}_reduction"]
            if not math.isfinite(value) or value <= 0.01:
                failures.append(
                    f"{family.upper()} correction response failed on {level}"
                )
    for model in MODELS:
        ratio = metrics[f"{model}_convergence_ratio"]
        if not math.isfinite(ratio) or ratio >= 1.0:
            failures.append(
                f"{model} grid differences do not contract ({ratio:.3g})"
            )
        value = metrics[f"{model}_257_to_513_relative_rms"]
        if not math.isfinite(value) or value > 0.05:
            failures.append(
                f"{model} medium-to-fine relative RMS exceeds 5% ({value:.3g})"
            )

    for name, value in metrics.items():
        print(f"{name}={value:.6g}")
    if failures:
        for failure in failures:
            print(f"ERROR: {failure}", file=sys.stderr)
        return 1
    print("All multi-grid validation checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
