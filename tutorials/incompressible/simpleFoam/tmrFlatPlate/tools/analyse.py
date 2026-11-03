#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Analyse the incompressible TMR 2DZP flat-plate benchmark."""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
NU = 2.0e-7
CASE_NAMES = ("sst", "sstcc-off", "sstcc-on", "sa", "sacc-off", "sacc-on")


def latest_time(case: Path) -> Path:
    times = []
    for path in case.iterdir():
        try:
            value = float(path.name)
        except (ValueError, TypeError):
            continue
        if path.is_dir() and value > 0:
            times.append((value, path))
    if not times:
        raise RuntimeError(f"no result time found in {case}")
    return max(times)[1]


def list_block(text: str, kind: str, field: Path) -> str:
    match = re.search(
        rf"internalField\s+nonuniform\s+List<{kind}>\s+\d+\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if not match:
        raise RuntimeError(f"cannot read {kind} internalField from {field}")
    return match.group(1)


def internal_vectors(field: Path) -> list[tuple[float, float, float]]:
    block = list_block(field.read_text(), "vector", field)
    return [
        tuple(map(float, values))
        for values in re.findall(rf"\(\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\)", block)
    ]


def internal_scalars(field: Path) -> list[float]:
    text = field.read_text()
    uniform = re.search(rf"internalField\s+uniform\s+({FLOAT})\s*;", text)
    if uniform:
        return [float(uniform.group(1))]
    return [float(value) for value in re.findall(FLOAT, list_block(text, "scalar", field))]


def mesh_points(path: Path) -> list[tuple[float, float, float]]:
    return [
        tuple(map(float, values))
        for values in re.findall(
            rf"\(\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\)", path.read_text()
        )
    ]


def mesh_faces(path: Path) -> list[list[int]]:
    return [
        list(map(int, values.split()))
        for values in re.findall(r"\d+\(\s*([0-9\s]+?)\s*\)", path.read_text())
    ]


def patch_block(path: Path, patch: str) -> str:
    match = re.search(
        rf"\b{re.escape(patch)}\s*\{{(.*?)\n\s*\}}", path.read_text(), re.S
    )
    if not match:
        raise RuntimeError(f"patch {patch} not found in {path}")
    return match.group(1)


def patch_range(boundary: Path, patch: str) -> tuple[int, int]:
    block = patch_block(boundary, patch)
    count = int(re.search(r"\bnFaces\s+(\d+)\s*;", block).group(1))
    start = int(re.search(r"\bstartFace\s+(\d+)\s*;", block).group(1))
    return start, count


def patch_vectors(field: Path, patch: str) -> list[tuple[float, float, float]]:
    return [
        tuple(map(float, values))
        for values in re.findall(
            rf"\(\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\)",
            patch_block(field, patch),
        )
    ]


def read_plot3d(path: Path) -> tuple[int, int, list[float], list[float]]:
    tokens = path.read_text().split()
    if int(tokens[0]) != 1:
        raise RuntimeError("the flat-plate analyser expects one PLOT3D block")
    ni, nj = map(int, tokens[1:3])
    values = list(map(float, tokens[3:]))
    count = ni * nj
    if len(values) != 2 * count:
        raise RuntimeError(f"unexpected PLOT3D value count in {path}")
    return ni, nj, values[:count], values[count:]


def cell_centres(
    ni: int, nj: int, x: list[float], y: list[float]
) -> tuple[list[float], list[float]]:
    cx, cy = [], []
    for j in range(nj - 1):
        for i in range(ni - 1):
            ids = (j * ni + i, j * ni + i + 1, (j + 1) * ni + i, (j + 1) * ni + i + 1)
            cx.append(sum(x[index] for index in ids) / 4.0)
            cy.append(sum(y[index] for index in ids) / 4.0)
    return cx, cy


def extract_cf(case: Path) -> list[tuple[float, float]]:
    mesh = case / "constant" / "polyMesh"
    result = latest_time(case)
    points = mesh_points(mesh / "points")
    faces = mesh_faces(mesh / "faces")
    start, count = patch_range(mesh / "boundary", "plate")
    shear = patch_vectors(result / "wallShearStress", "plate")
    if len(shear) != count:
        raise RuntimeError(f"expected {count} wall values in {case}, found {len(shear)}")
    curve = []
    for face, tau in zip(faces[start : start + count], shear):
        x_coord = sum(points[index][0] for index in face) / len(face)
        # rho=1 and U_inf=1, so Cf=2*|tau_w|.
        curve.append((x_coord, 2.0 * math.hypot(tau[0], tau[1])))
    return sorted(curve)


def profiles(
    case: Path, grid: Path, curve: list[tuple[float, float]]
) -> list[dict[str, object]]:
    ni, nj, x, y = read_plot3d(grid)
    cx, cy = cell_centres(ni, nj, x, y)
    velocity = internal_vectors(latest_time(case) / "U")
    nci, ncj = ni - 1, nj - 1
    if len(velocity) != nci * ncj:
        raise RuntimeError(f"field/grid size mismatch in {case}")
    result = []
    for wall_x, cf in curve:
        i = min(range(nci), key=lambda index: abs(cx[index] - wall_x))
        samples = [(0.0, 0.0)]
        for j in range(ncj):
            index = j * nci + i
            samples.append((cy[index], velocity[index][0]))
        samples.sort()
        edge_speed = max(speed for _, speed in samples)
        normalized = [
            (height, min(max(speed / edge_speed, 0.0), 1.0))
            for height, speed in samples
        ]
        edge_index = next(
            (
                index
                for index, (_, speed) in enumerate(normalized)
                if index > 0 and speed >= 0.999
            ),
            len(normalized) - 1,
        )
        integration_samples = normalized[: edge_index + 1]
        theta = sum(
            0.5
            * ((ua * (1.0 - ua)) + (ub * (1.0 - ub)))
            * (yb - ya)
            for (ya, ua), (yb, ub) in zip(
                integration_samples, integration_samples[1:]
            )
        )
        utau = math.sqrt(max(cf, 0.0) / 2.0)
        result.append(
            {
                "x": wall_x,
                "cf": cf,
                "retheta": theta / NU,
                "upyp": [
                    (height * utau / NU, speed / utau)
                    for height, speed in samples
                    if height > 0 and utau > 0
                ],
            }
        )
    return result


def numeric_rows(path: Path) -> list[list[float]]:
    rows = []
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("#") or "VARIABLE" in line.upper() or "ZONE" in line.upper():
            continue
        values = re.findall(FLOAT, line)
        if len(values) >= 2:
            rows.append([float(value) for value in values])
    return rows


def interpolate(data: list[tuple[float, float]], x: float) -> float:
    data = sorted(data)
    if not data[0][0] <= x <= data[-1][0]:
        raise ValueError(x)
    for (xa, ya), (xb, yb) in zip(data, data[1:]):
        if xa <= x <= xb:
            fraction = (x - xa) / max(xb - xa, 1e-300)
            return ya + fraction * (yb - ya)
    return data[-1][1]


def errors(
    predicted: list[tuple[float, float]], reference: list[tuple[float, float]]
) -> tuple[float, float, float]:
    pairs = []
    for x, value in predicted:
        try:
            expected = interpolate(reference, x)
        except ValueError:
            continue
        pairs.append((value, expected))
    if not pairs:
        raise RuntimeError("no overlap with reference data")
    rmse = math.sqrt(sum((value - expected) ** 2 for value, expected in pairs) / len(pairs))
    scale = math.sqrt(sum(expected**2 for _, expected in pairs) / len(pairs))
    bias = sum(value - expected for value, expected in pairs) / len(pairs)
    return rmse, rmse / max(scale, 1e-300), bias


def curve_difference(
    baseline: list[tuple[float, float]], candidate: list[tuple[float, float]]
) -> tuple[float, float]:
    if len(baseline) != len(candidate):
        raise RuntimeError("curves use different meshes")
    relative = [
        abs(b - a) / max(abs(a), 1e-300)
        for (_, a), (_, b) in zip(baseline, candidate)
    ]
    return max(relative), math.sqrt(sum(value * value for value in relative) / len(relative))


def analyse_case(
    case: Path,
    grid: Path,
    ks: list[tuple[float, float]],
    cfd_reference: list[tuple[float, float]],
    coles: list[tuple[float, float]],
) -> tuple[list[tuple[float, float]], dict[str, float]]:
    curve = extract_cf(case)
    sampled = profiles(case, grid, curve)
    cf_retheta = [
        (float(row["retheta"]), float(row["cf"]))
        for row in sampled
        if 4000.0 <= float(row["retheta"]) <= 13000.0
    ]
    ks_rmse, ks_nrmse, ks_bias = errors(cf_retheta, ks)
    cfd_rmse, cfd_nrmse, _ = errors(cf_retheta, cfd_reference)
    station = min(sampled, key=lambda row: abs(float(row["retheta"]) - 10000.0))
    upyp = [(math.log10(x), y) for x, y in station["upyp"] if x > 0]
    coles_log = [(math.log10(x), y) for x, y in coles if x > 0]
    up_rmse, up_nrmse, _ = errors(upyp, coles_log)
    metrics = {
        "ks_cf_rmse": ks_rmse,
        "ks_cf_nrmse": ks_nrmse,
        "ks_cf_bias": ks_bias,
        "cfl3d_cf_rmse": cfd_rmse,
        "cfl3d_cf_nrmse": cfd_nrmse,
        "coles_up_rmse": up_rmse,
        "coles_up_nrmse": up_nrmse,
        "profile_retheta": float(station["retheta"]),
        "profile_x": float(station["x"]),
    }
    with (case / "cf-retheta.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("Re_theta", "Cf"))
        writer.writerows(cf_retheta)
    with (case / "upyp.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("y_plus", "u_plus"))
        writer.writerows(station["upyp"])
    return curve, metrics


def main() -> int:
    tutorial = Path(sys.argv[1]).resolve()
    grid = Path(sys.argv[2]).resolve()
    runs = tutorial / "runs"
    reference = tutorial / "resources" / "reference"
    ks = [(row[0], row[1]) for row in numeric_rows(reference / "cf_K-S.dat")]
    coles = [(row[0], row[1]) for row in numeric_rows(reference / "coles-upyp.dat")]
    cfd = {
        "sst": [(row[1], row[0]) for row in numeric_rows(reference / "sst-cf-cfl3d.dat")],
        "sa": [(row[1], row[0]) for row in numeric_rows(reference / "sa-cf-cfl3d.dat")],
    }
    curves, metrics = {}, {}
    for name in CASE_NAMES:
        family = "sst" if name.startswith("sst") else "sa"
        curves[name], metrics[name] = analyse_case(
            runs / name, grid, ks, cfd[family], coles
        )
        fr1 = latest_time(runs / name) / "fr1"
        if fr1.exists():
            values = internal_scalars(fr1)
            if not all(math.isfinite(value) for value in values):
                raise RuntimeError(f"non-finite fr1 in {name}")
            metrics[name]["fr1_min"] = min(values)
            metrics[name]["fr1_max"] = max(values)
            metrics[name]["fr1_mean_abs_departure"] = sum(
                abs(value - 1.0) for value in values
            ) / len(values)

    comparisons = {}
    for label, first, second in (
        ("sst_off", "sst", "sstcc-off"),
        ("sa_off", "sa", "sacc-off"),
        ("sst_on", "sst", "sstcc-on"),
        ("sa_on", "sa", "sacc-on"),
    ):
        comparisons[f"{label}_cf_max_relative"] , comparisons[
            f"{label}_cf_rms_relative"
        ] = curve_difference(curves[first], curves[second])

    with (runs / "summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("case", "metric", "value"))
        for name in CASE_NAMES:
            for metric, value in sorted(metrics[name].items()):
                writer.writerow((name, metric, f"{value:.12g}"))
        for metric, value in sorted(comparisons.items()):
            writer.writerow(("comparison", metric, f"{value:.12g}"))

    failures = []
    for family in ("sst", "sa"):
        if comparisons[f"{family}_off_cf_max_relative"] > 1e-10:
            failures.append(f"{family.upper()} correction-off does not recover native Cf")
        if comparisons[f"{family}_on_cf_max_relative"] > 0.01:
            failures.append(f"{family.upper()} RC changes neutral flat-plate Cf by more than 1%")
    for name in ("sstcc-on", "sacc-on"):
        if metrics[name]["ks_cf_nrmse"] > 0.15:
            failures.append(f"{name} differs from Karman-Schoenherr by more than 15% nRMSE")

    for metric, value in sorted(comparisons.items()):
        print(f"{metric}={value:.6g}")
    for name in CASE_NAMES:
        print(
            f"{name}: KS nRMSE={metrics[name]['ks_cf_nrmse']:.4g}, "
            f"CFL3D nRMSE={metrics[name]['cfl3d_cf_nrmse']:.4g}, "
            f"Coles nRMSE={metrics[name]['coles_up_nrmse']:.4g}"
        )
        if "fr1_min" in metrics[name]:
            print(
                f"{name}: fr1=[{metrics[name]['fr1_min']:.6g}, "
                f"{metrics[name]['fr1_max']:.6g}]"
            )
    if failures:
        for failure in failures:
            print(f"ERROR: {failure}", file=sys.stderr)
        return 1
    print("TMR 2DZP neutral benchmark passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
