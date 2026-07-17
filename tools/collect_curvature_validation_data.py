#!/usr/bin/env python3
"""Collect reproducible, lightweight validation data from completed tutorials."""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
from pathlib import Path
import sys


def interpolate(curve: list[tuple[float, float]], x_value: float) -> float:
    x = [item[0] for item in curve]
    index = bisect.bisect_left(x, x_value)
    if index == 0:
        return curve[0][1]
    if index == len(curve):
        return curve[-1][1]
    x0, y0 = curve[index - 1]
    x1, y1 = curve[index]
    return y0 + (x_value - x0) * (y1 - y0) / (x1 - x0)


def at_bend(curve: list[tuple[float, float]]) -> float:
    return min(curve, key=lambda item: abs(item[0]))[1]


def write_profile(
    path: Path,
    baseline: list[tuple[float, float]],
    corrected: list[tuple[float, float]],
    reference: list[tuple[float, float]],
) -> None:
    xmin, xmax = reference[0][0], reference[-1][0]
    rows = [
        (x_value, base, correction)
        for (x_value, base), (_, correction) in zip(baseline, corrected)
        if xmin <= x_value <= xmax
    ]
    if len(rows) < 2:
        raise RuntimeError(f"fewer than two model points overlap reference range [{xmin}, {xmax}]")
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("x", "baseline", "corrected", "reference", "weight"))
        for index, (x_value, base, correction) in enumerate(rows):
            left = rows[max(0, index - 1)][0]
            right = rows[min(len(rows) - 1, index + 1)][0]
            weight = 0.5 * (right - left)
            writer.writerow((f"{x_value:.12g}", f"{base:.12g}", f"{correction:.12g}", f"{interpolate(reference, x_value):.12g}", f"{weight:.12g}"))


def correction_from_indicator(indicator: float, cr2: float, *, sst: bool) -> float:
    raw = 2.0 * (1.0 - math.atan(cr2 * indicator)) - 1.0
    return min(max(raw, 0.0), 1.25) if sst else raw


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tutorial",
        type=Path,
        default=Path("tutorials/incompressible/simpleFoam/smitsCurve"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/curvature-validation"),
    )
    args = parser.parse_args()
    tutorial = args.tutorial.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(tutorial / "tools"))
    from analyse import extract_cf, read_experiment  # type: ignore

    results = tutorial / "runs/grid-study"
    summary_path = results / "summary.csv"
    if not summary_path.is_file():
        print(f"ERROR: run Allrun-gridStudy first; missing {summary_path}", file=sys.stderr)
        return 2

    levels = ("level3-129x49", "level2-257x97", "level1-513x193")
    curves: dict[str, dict[str, list[tuple[float, float]]]] = {}
    for level in levels:
        curves[level] = {}
        for model in ("sst", "sstcc-on", "sa", "sacc-on"):
            curves[level][model], _ = extract_cf(results / level / model)

    experiment = read_experiment(
        tutorial / "resources/reference/Smits_curvature_expdata_reworked2012.dat"
    )
    fine = curves[levels[-1]]
    write_profile(output / "smits-sst-fine.csv", fine["sst"], fine["sstcc-on"], experiment)
    write_profile(output / "smits-sa-fine.csv", fine["sa"], fine["sacc-on"], experiment)

    with (output / "smits-grid-bend.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("level", "cells", "sst", "sst_rc", "sa", "sa_rc"))
        for level, cells in zip(levels, (6144, 24576, 98304)):
            writer.writerow(
                (
                    level,
                    cells,
                    at_bend(curves[level]["sst"]),
                    at_bend(curves[level]["sstcc-on"]),
                    at_bend(curves[level]["sa"]),
                    at_bend(curves[level]["sacc-on"]),
                )
            )

    summary = {}
    with summary_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            summary[row["metric"]] = float(row["value"])
    (output / "smits-grid-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )

    indicators = [value / 20 for value in range(-40, 41)]
    for family, default, perturbed in (("sst", 2.0, 1.8), ("sa", 12.0, 10.8)):
        with (output / f"{family}-cr2-sensitivity.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("x", "baseline", "corrected", "weight"))
            for indicator in indicators:
                writer.writerow(
                    (
                        indicator,
                        correction_from_indicator(indicator, perturbed, sst=family == "sst"),
                        correction_from_indicator(indicator, default, sst=family == "sst"),
                        1,
                    )
                )

    print(f"Collected validation data in {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
