#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare reconstructed OpenFOAM volume fields after restart/MPI runs."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def latest_time(case: Path) -> Path:
    times = []
    for path in case.iterdir():
        if not path.is_dir():
            continue
        try:
            value = float(path.name)
        except ValueError:
            continue
        if value > 0:
            times.append((value, path))
    if not times:
        raise RuntimeError(f"no result time found in {case}")
    return max(times)[1]


def internal_values(path: Path) -> list[float]:
    text = path.read_text()
    uniform_vector = re.search(
        rf"internalField\s+uniform\s+\(\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\)\s*;",
        text,
    )
    if uniform_vector:
        return [float(value) for value in uniform_vector.groups()]

    uniform_scalar = re.search(rf"internalField\s+uniform\s+({FLOAT})\s*;", text)
    if uniform_scalar:
        return [float(uniform_scalar.group(1))]

    vector_values = re.search(
        r"internalField\s+nonuniform\s+List<vector>\s+\d+\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if vector_values:
        return [
            float(value)
            for vector in re.findall(
                rf"\(\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\)",
                vector_values.group(1),
            )
            for value in vector
        ]

    scalar_values = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+\d+\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if scalar_values:
        return [float(value) for value in re.findall(FLOAT, scalar_values.group(1))]

    raise RuntimeError(f"cannot read internalField from {path}")


def compare(reference: Path, candidate: Path) -> tuple[float, float]:
    first = internal_values(reference)
    second = internal_values(candidate)
    if len(first) != len(second):
        raise RuntimeError(
            f"field sizes differ: {reference} has {len(first)}, "
            f"{candidate} has {len(second)}"
        )
    if not all(math.isfinite(value) for value in first + second):
        raise RuntimeError(f"non-finite value in {reference} or {candidate}")
    maximum = max((abs(a - b) for a, b in zip(first, second)), default=0.0)
    scale = max((abs(value) for value in first), default=0.0)
    return maximum, maximum/max(scale, 1e-30)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--candidate", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tolerance", type=float, default=5e-5)
    parser.add_argument("fields", nargs="+")
    args = parser.parse_args()

    reference_time = latest_time(args.reference)
    rows = []
    failed = False
    for case in args.candidate:
        candidate_time = latest_time(case)
        for field in args.fields:
            absolute, relative = compare(
                reference_time / field, candidate_time / field
            )
            passed = relative <= args.tolerance
            failed |= not passed
            rows.append((case.name, field, absolute, relative, passed))
            print(
                f"{case.name}/{field}: abs={absolute:.6g}, "
                f"rel={relative:.6g}, pass={passed}"
            )

    with args.output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("candidate", "field", "max_abs", "relative_linf", "pass"))
        writer.writerows(
            (case, field, f"{absolute:.12g}", f"{relative:.12g}", str(passed).lower())
            for case, field, absolute, relative, passed in rows
        )
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
