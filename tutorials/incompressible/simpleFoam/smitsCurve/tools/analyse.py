#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Extract convex-wall Cf, validate ON/OFF behaviour, and draw an SVG."""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
UREF = 31.9
CASE_NAMES = ("sst", "sstcc-off", "sstcc-on", "sa", "sacc-off", "sacc-on")


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


def read_points(path: Path) -> list[tuple[float, float, float]]:
    text = path.read_text()
    return [
        tuple(map(float, values))
        for values in re.findall(rf"\(\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\)", text)
    ]


def read_faces(path: Path) -> list[list[int]]:
    text = path.read_text()
    return [
        list(map(int, values.split()))
        for values in re.findall(r"\d+\(\s*([0-9\s]+?)\s*\)", text)
    ]


def patch_range(boundary: Path, patch: str) -> tuple[int, int]:
    text = boundary.read_text()
    match = re.search(rf"\b{re.escape(patch)}\s*\{{(.*?)\n\s*\}}", text, re.S)
    if not match:
        raise RuntimeError(f"patch {patch} not found in {boundary}")
    block = match.group(1)
    count = int(re.search(r"\bnFaces\s+(\d+)\s*;", block).group(1))
    start = int(re.search(r"\bstartFace\s+(\d+)\s*;", block).group(1))
    return start, count


def patch_vectors(field: Path, patch: str) -> list[tuple[float, float, float]]:
    text = field.read_text()
    match = re.search(rf"\b{re.escape(patch)}\s*\{{(.*?)\n\s*\}}", text, re.S)
    if not match:
        raise RuntimeError(f"patch {patch} not found in {field}")
    block = match.group(1)
    return [
        tuple(map(float, values))
        for values in re.findall(rf"\(\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\)", block)
    ]


def internal_scalars(field: Path) -> list[float]:
    text = field.read_text()
    uniform = re.search(rf"internalField\s+uniform\s+({FLOAT})\s*;", text)
    if uniform:
        return [float(uniform.group(1))]
    values = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+\d+\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if not values:
        raise RuntimeError(f"cannot read scalar internalField from {field}")
    return [float(value) for value in re.findall(FLOAT, values.group(1))]


def extract_cf(case: Path) -> tuple[list[tuple[float, float]], tuple[float, float] | None]:
    mesh = case / "constant" / "polyMesh"
    time = latest_time(case)
    points = read_points(mesh / "points")
    faces = read_faces(mesh / "faces")
    start, count = patch_range(mesh / "boundary", "convexWall")
    shear = patch_vectors(time / "wallShearStress", "convexWall")
    if len(shear) != count:
        raise RuntimeError(f"expected {count} wall-shear values in {case}, found {len(shear)}")

    curve = []
    for face, tau in zip(faces[start : start + count], shear):
        x = sum(points[index][0] for index in face) / len(face)
        magnitude = math.sqrt(sum(component * component for component in tau))
        curve.append((x, 2.0 * magnitude / (UREF * UREF)))
    curve.sort()

    fr_field = time / "fr1"
    fr_range = None
    if fr_field.exists():
        values = internal_scalars(fr_field)
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError(f"non-finite fr1 values in {case}")
        fr_range = (min(values), max(values))
    return curve, fr_range


def read_reference(path: Path) -> list[tuple[float, float]]:
    result = []
    for line in path.read_text().splitlines():
        values = re.findall(FLOAT, line)
        if len(values) == 2 and not line.lstrip().startswith("#"):
            result.append((float(values[0]), float(values[1])))
    return result


def read_experiment(path: Path) -> list[tuple[float, float]]:
    text = path.read_text()
    marker = 'variables="s, m","cf_ue","cf_uref (approx)"'
    section = text.split(marker, 1)[1].split("----------", 1)[0]
    result = []
    for line in section.splitlines():
        values = re.findall(FLOAT, line)
        if len(values) != 3:
            continue
        s_coord, _, cf_ref = map(float, values)
        x_coord = (
            -0.0635 + (s_coord + 0.0665) * math.cos(math.radians(30.0))
            if s_coord < -0.066497
            else s_coord
        )
        result.append((x_coord, cf_ref))
    return result


def relative_curve_error(first: list[tuple[float, float]], second: list[tuple[float, float]]) -> float:
    if len(first) != len(second):
        raise RuntimeError("curves use different meshes")
    denominator = max(max(abs(value) for _, value in first), 1e-30)
    return max(abs(a[1] - b[1]) for a, b in zip(first, second)) / denominator


def near_zero_reduction(base: list[tuple[float, float]], corrected: list[tuple[float, float]]) -> float:
    index = min(range(len(base)), key=lambda item: abs(base[item][0]))
    base_value = base[index][1]
    return (base_value - corrected[index][1]) / max(abs(base_value), 1e-30)


def write_svg(path: Path, series: list[tuple[str, list[tuple[float, float]], str, bool]]) -> None:
    width, height = 1000, 620
    left, right, top, bottom = 90, 30, 40, 75
    xmin, xmax = -0.25, 1.3
    visible = [y for _, data, _, _ in series for x, y in data if xmin <= x <= xmax]
    ymax = max(visible) * 1.08

    def sx(value: float) -> float:
        return left + (value - xmin) * (width - left - right) / (xmax - xmin)

    def sy(value: float) -> float:
        return height - bottom - value * (height - top - bottom) / ymax

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="24" text-anchor="middle" font-family="sans-serif" font-size="18">Smits convex wall: skin-friction coefficient</text>',
    ]
    for tick in range(7):
        value = ymax * tick / 6
        y = sy(value)
        lines.append(f'<line x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}" stroke="#ddd"/>')
        lines.append(f'<text x="{left-8}" y="{y+4:.1f}" text-anchor="end" font-family="sans-serif" font-size="11">{value:.4f}</text>')
    for tick in range(7):
        value = xmin + (xmax - xmin) * tick / 6
        x = sx(value)
        lines.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{height-bottom}" stroke="#eee"/>')
        lines.append(f'<text x="{x:.1f}" y="{height-bottom+20}" text-anchor="middle" font-family="sans-serif" font-size="11">{value:.2f}</text>')
    lines.extend([
        f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black"/>',
        f'<text x="{(left+width-right)/2}" y="{height-25}" text-anchor="middle" font-family="sans-serif">x [m]</text>',
        f'<text transform="translate(24 {(top+height-bottom)/2}) rotate(-90)" text-anchor="middle" font-family="sans-serif">Cf based on Uref</text>',
    ])

    for index, (name, data, colour, points_only) in enumerate(series):
        visible_data = [(x, y) for x, y in data if xmin <= x <= xmax]
        if points_only:
            for x, y in visible_data:
                lines.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="4" fill="{colour}"/>')
        else:
            coords = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in visible_data)
            dash = ' stroke-dasharray="7,4"' if "TMBWG" in name else ""
            lines.append(f'<polyline points="{coords}" fill="none" stroke="{colour}" stroke-width="2"{dash}/>')
        legend_x = 610 + (index % 2) * 190
        legend_y = 52 + (index // 2) * 17
        lines.append(f'<rect x="{legend_x}" y="{legend_y-9}" width="18" height="3" fill="{colour}"/>')
        lines.append(f'<text x="{legend_x+24}" y="{legend_y-4}" font-family="sans-serif" font-size="11">{name}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    tutorial = Path(sys.argv[1]).resolve()
    runs = tutorial / "runs"
    curves = {}
    fr_ranges = {}
    for name in CASE_NAMES:
        curves[name], fr_ranges[name] = extract_cf(runs / name)
        with (runs / name / "cf.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("x_m", "Cf_Uref"))
            writer.writerows(curves[name])

    checks = {
        "sst_off_relative_error": relative_curve_error(curves["sst"], curves["sstcc-off"]),
        "sa_off_relative_error": relative_curve_error(curves["sa"], curves["sacc-off"]),
        "sst_near_bend_reduction": near_zero_reduction(curves["sst"], curves["sstcc-on"]),
        "sa_near_bend_reduction": near_zero_reduction(curves["sa"], curves["sacc-on"]),
    }
    with (runs / "summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("metric", "value"))
        for key, value in checks.items():
            writer.writerow((key, f"{value:.12g}"))
        for name in CASE_NAMES:
            if fr_ranges[name] is not None:
                writer.writerow((f"{name}_fr1_min", f"{fr_ranges[name][0]:.12g}"))
                writer.writerow((f"{name}_fr1_max", f"{fr_ranges[name][1]:.12g}"))

    reference = tutorial / "resources" / "reference"
    series = [
        ("OpenFOAM SST", curves["sst"], "#1f77b4", False),
        ("OpenFOAM SST-RC", curves["sstcc-on"], "#d62728", False),
        ("OpenFOAM SA", curves["sa"], "#2ca02c", False),
        ("OpenFOAM SA-RC", curves["sacc-on"], "#9467bd", False),
        ("TMBWG SST", read_reference(reference / "smits_cf_cfl3d_sst.dat"), "#6baed6", False),
        ("TMBWG SST-RC", read_reference(reference / "smits_cf_cfl3d_sstrc.dat"), "#fb6a4a", False),
        ("TMBWG SA", read_reference(reference / "smits_cf_cfl3d_sa.dat"), "#74c476", False),
        ("TMBWG SA-RC", read_reference(reference / "smits_cf_cfl3d_sarc.dat"), "#bcbddc", False),
    ]
    write_svg(runs / "cf-comparison.svg", series)

    failures = []
    if checks["sst_off_relative_error"] > 1e-5:
        failures.append("SSTCC-off differs from native SST")
    if checks["sa_off_relative_error"] > 1e-5:
        failures.append("SACC-off differs from native SA")
    if checks["sst_near_bend_reduction"] <= 0.01:
        failures.append("SST-RC does not reduce convex-wall Cf near the bend")
    if checks["sa_near_bend_reduction"] <= 0.01:
        failures.append("SA-RC does not reduce convex-wall Cf near the bend")

    sst_off = fr_ranges["sstcc-off"]
    sst_on = fr_ranges["sstcc-on"]
    sa_off = fr_ranges["sacc-off"]
    sa_on = fr_ranges["sacc-on"]
    if sst_off is None or max(abs(sst_off[0] - 1), abs(sst_off[1] - 1)) > 1e-12:
        failures.append("SSTCC-off fr1 is not uniformly one")
    if sa_off is None or max(abs(sa_off[0] - 1), abs(sa_off[1] - 1)) > 1e-12:
        failures.append("SACC-off fr1 is not uniformly one")
    if sst_on is None or sst_on[0] < -1e-12 or sst_on[1] > 1.25 + 1e-12 or sst_on[1] - sst_on[0] < 1e-3:
        failures.append("SST-RC fr1 is not finite, bounded, and non-uniform")
    if sa_on is None or sa_on[1] - sa_on[0] < 1e-3:
        failures.append("SA-RC fr1 is not finite and non-uniform")

    for key, value in checks.items():
        print(f"{key}={value:.6g}")
    for name in CASE_NAMES:
        if fr_ranges[name] is not None:
            print(f"{name}: fr1=[{fr_ranges[name][0]:.6g}, {fr_ranges[name][1]:.6g}]")
    if failures:
        for failure in failures:
            print(f"ERROR: {failure}", file=sys.stderr)
        return 1
    print("All curvature-correction validation checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
