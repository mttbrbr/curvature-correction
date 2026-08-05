#!/usr/bin/env python3
"""Analyse TMR NACA 0012 SA/SA-RC forces against CC0 CFL3D data."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import re
import sys


ANGLES = (0, 10, 15)
MODELS = ("sa", "sarc")
ACCURACY_LIMITS = {
    "quick": {0: (None, 1.10), 10: (0.12, 0.45)},
    "medium": {0: (None, 0.25), 10: (0.05, 0.10)},
}


def numeric_rows(path: Path, columns: int) -> list[tuple[float, ...]]:
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) != columns:
            continue
        try:
            rows.append(tuple(map(float, fields)))
        except ValueError:
            continue
    if not rows:
        raise RuntimeError(f"no numeric rows in {path}")
    return rows


def force_coefficients(case: Path) -> tuple[float, float, int, float, float]:
    files = sorted(
        case.glob("postProcessing/forceCoeffs/*/coefficient.dat"),
        key=lambda path: float(path.parent.name),
    )
    if not files:
        raise RuntimeError(f"forceCoeffs output absent in {case}")
    rows = numeric_rows(files[-1], 13)
    tail = rows[-min(5, len(rows)):]
    cd = sum(row[1] for row in tail) / len(tail)
    cl = sum(row[4] for row in tail) / len(tail)
    cl_span = max(row[4] for row in tail) - min(row[4] for row in tail)
    cd_span = max(row[1] for row in tail) - min(row[1] for row in tail)
    log = (case / "log.simpleFoam").read_text()
    iterations = len(re.findall(r"^Time = ", log, flags=re.MULTILINE))
    if not all(math.isfinite(value) for value in (cl, cd)):
        raise RuntimeError(f"non-finite force coefficient in {case}")
    return cl, cd, iterations, cl_span, cd_span


def references(directory: Path) -> dict[str, dict[int, tuple[float, float]]]:
    result = {}
    for model in MODELS:
        rows = numeric_rows(directory / f"n0012clcd-cfl3d-{model}.dat", 3)
        result[model] = {
            int(round(alpha)): (cl, cd)
            for alpha, cl, cd in rows
            if int(round(alpha)) in ANGLES
        }
    return result


def relative(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1e-12)


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: analyse.py TUTORIAL GRID_LEVEL", file=sys.stderr)
        return 2
    tutorial = Path(sys.argv[1]).resolve()
    level = sys.argv[2]
    results = tutorial / "runs" / level
    reference = references(tutorial / "resources" / "reference")
    rows = []
    failures = []
    angles = tuple(
        angle
        for angle in ANGLES
        if all((results / f"alpha{angle}-{model}").is_dir() for model in MODELS)
    )
    if not angles:
        raise RuntimeError(f"no complete SA/SA-RC angle pairs in {results}")

    for angle in angles:
        local = {}
        for model in MODELS:
            case = results / f"alpha{angle}-{model}"
            cl, cd, iterations, cl_span, cd_span = force_coefficients(case)
            ref_cl, ref_cd = reference[model][angle]
            local[model] = (cl, cd)
            cl_error = 0.0 if angle == 0 else relative(cl, ref_cl)
            cd_error = relative(cd, ref_cd)
            rows.append(
                {
                    "grid": level,
                    "angle": angle,
                    "model": model,
                    "cl": cl,
                    "cd": cd,
                    "reference_cl": ref_cl,
                    "reference_cd": ref_cd,
                    "cl_relative_error": cl_error,
                    "cd_relative_error": cd_error,
                    "iterations": iterations,
                    "cl_tail_span": cl_span,
                    "cd_tail_span": cd_span,
                    "acceptance_role": "gate" if angle in (0, 10) else "diagnostic",
                }
            )
            if angle in (0, 10) and cd_span / max(abs(cd), 1e-12) > 0.08:
                failures.append(
                    f"alpha={angle} {model}: final CD tail span exceeds 8%"
                )
            if angle == 0 and cl_span > 0.01:
                failures.append(f"alpha=0 {model}: final CL tail span exceeds 0.01")
            if angle == 10 and cl_span / max(abs(cl), 1e-12) > 0.01:
                failures.append(
                    f"alpha=10 {model}: final CL tail span exceeds 1%"
                )
            limits = ACCURACY_LIMITS.get(level, {}).get(angle)
            if not (-3.0 < cl < 3.0 and 0.0 < cd < 1.0):
                failures.append(
                    f"alpha={angle} {model}: coefficients outside finite sanity bounds"
                )
            if angle == 0 and abs(cl) > 0.02:
                failures.append(f"alpha=0 {model}: absolute CL exceeds 0.02")
            if limits is not None:
                cl_limit, cd_limit = limits
                if cl_limit is not None and cl_error > cl_limit:
                    failures.append(
                        f"alpha={angle} {model}: CL relative error exceeds {cl_limit:.0%}"
                    )
                if cd_error > cd_limit:
                    failures.append(
                        f"alpha={angle} {model}: CD relative error exceeds {cd_limit:.0%}"
                    )

        if angle == 10:
            if local["sarc"][0] < local["sa"][0]:
                failures.append(f"alpha={angle}: SA-RC does not increase CL")
            if local["sarc"][1] > local["sa"][1]:
                failures.append(f"alpha={angle}: SA-RC does not reduce CD")

    with (results / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{row['grid']} alpha={row['angle']:>2} {row['model']}: "
            f"CL={row['cl']:.6g} (ref {row['reference_cl']:.6g}), "
            f"CD={row['cd']:.6g} (ref {row['reference_cd']:.6g})"
        )
    if failures:
        for failure in failures:
            print(f"ERROR: {failure}", file=sys.stderr)
        return 1
    print(f"TMR NACA 0012 {level} benchmark passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
