#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Analyse the TMR rotating-pipe benchmark and its frame equivalence."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re
from typing import Iterable

UREF = 10.0
RADIUS = 0.03
FRAME_RMS_LIMIT = 0.01
FRAME_MAX_LIMIT = 0.02
SCHEME_RMS_LIMIT = 0.02
SCHEME_MAX_LIMIT = 0.04
OFF_LIMIT = 1.0e-10
CONTINUITY_LIMIT = 1.0e-5
LOCAL_CONTINUITY_LIMIT = 1.0e-4
RESIDUAL_LIMIT = 1.0e-5
NORMALIZATION_FLOOR = 1.0e-3
SAMPLE_FIELD_COMPONENTS = {
    "curvatureIndicator": 1,
    "fr1": 1,
    "fr1Clipped": 1,
    "k": 1,
    "nut": 1,
    "productionBase": 1,
    "productionCorrected": 1,
    "productionRatio": 1,
    "OmegaRot": 3,
    "U": 3,
    "turbulenceProperties:R": 6,
}


def numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        try:
            row = [float(field) for field in fields]
        except ValueError:
            continue
        if not all(math.isfinite(value) for value in row):
            raise RuntimeError(f"non-finite value in {path}")
        rows.append(row)
    if not rows:
        raise RuntimeError(f"no numeric rows in {path}")
    return rows


def parse_tecplot_zones(path: Path) -> dict[str, list[tuple[float, ...]]]:
    zones: dict[str, list[tuple[float, ...]]] = {}
    current: str | None = None
    for line in path.read_text().splitlines():
        match = re.match(r'ZONE\s+T="([^"]+)"', line.strip())
        if match:
            current = match.group(1)
            zones[current] = []
            continue
        if current is None or not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            values = tuple(float(field) for field in line.split())
        except ValueError:
            continue
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError(f"non-finite reference value in {path}")
        zones[current].append(values)
    if not zones or any(not values for values in zones.values()):
        raise RuntimeError(f"invalid Tecplot zones in {path}")
    return zones


def bulk_ratio_from_reference(rows: list[tuple[float, ...]]) -> float:
    """Area-integrate measured u/Uref and use an anchored 1/7 wall tail."""
    points = [(row[0], row[2]) for row in rows]
    measured = sum(
        (rb - ra) * (2 * ra * ua + 2 * rb * ub) / 2
        for (ra, ua), (rb, ub) in zip(points, points[1:])
    )
    radius, velocity = points[-1]
    wall_distance = 1 - radius
    coefficient = velocity / wall_distance ** (1 / 7)
    tail = 2 * coefficient * (
        wall_distance ** (8 / 7) / (8 / 7)
        - wall_distance ** (15 / 7) / (15 / 7)
    )
    return measured + tail


def latest_directory(root: Path) -> Path:
    candidates = []
    for path in root.iterdir():
        if path.is_dir():
            try:
                candidates.append((float(path.name), path))
            except ValueError:
                pass
    if not candidates:
        raise RuntimeError(f"no time directories below {root}")
    return max(candidates)[1]


def sample_file(case: Path, station: str, field: str) -> Path:
    time_dir = latest_directory(case / "postProcessing/sampleProfiles")
    matches = [
        path for path in time_dir.iterdir()
        if path.is_file() and path.name.startswith(station + "_") and field in path.name
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one {station}/{field} sample in {time_dir}, found {matches}"
        )
    return matches[0]


def sample_components(
    case: Path, station: str, field: str, components: int
) -> tuple[list[float], list[list[float]]]:
    path = sample_file(case, station, field)
    rows = numeric_rows(path)
    stem = path.name.removesuffix(path.suffix)
    prefix = station + "_"
    if not stem.startswith(prefix):
        raise RuntimeError(f"invalid sample filename {path.name}")
    sampled_fields = stem[len(prefix):].split("_")
    try:
        widths = [SAMPLE_FIELD_COMPONENTS[name] for name in sampled_fields]
    except KeyError as error:
        raise RuntimeError(f"unknown sampled field {error.args[0]} in {path}") from error
    if field not in sampled_fields:
        raise RuntimeError(f"field {field} is absent from {path}")
    field_index = sampled_fields.index(field)
    if widths[field_index] != components:
        raise RuntimeError(
            f"expected {components} components for {field}, found {widths[field_index]}"
        )
    coordinate_columns = len(rows[0]) - sum(widths)
    if coordinate_columns not in (1, 3):
        raise RuntimeError(
            f"unsupported sampled shape for {field}: {len(rows[0])} columns"
        )
    start = coordinate_columns + sum(widths[:field_index])
    end = start + components
    radii = []
    values = []
    for row in rows:
        if len(row) != coordinate_columns + sum(widths):
            raise RuntimeError(f"inconsistent sampled shape for {field}")
        if coordinate_columns == 1:
            radius = row[0] + 1.0e-6
        else:
            radius = math.hypot(row[1], row[2])
        radii.append(radius / RADIUS)
        values.append(row[start:end])
    return radii, values


def case_profile(case: Path) -> dict[str, list[float]]:
    rotation = case.name.split("-", 1)[0]
    station = "nonRotatingStation" if rotation == "n0" else "rotatingStation"
    radii, velocity = sample_components(case, station, "U", 3)
    r_r, stresses = sample_components(case, station, "turbulenceProperties:R", 6)
    if len(radii) != len(r_r):
        raise RuntimeError(f"U/R sample sizes differ in {case}")
    result = {
        "r": radii,
        "ux": [value[0] / UREF for value in velocity],
        "utheta": [value[2] / UREF for value in velocity],
        # symmTensor order: xx, xy, xz, yy, yz, zz.  At z=0 radial is +y.
        "rxr": [value[1] / UREF**2 for value in stresses],
    }
    family = case.name.split("-")[1]
    if family == "sst":
        r_k, values = sample_components(case, station, "k", 1)
        if len(r_k) != len(radii):
            raise RuntimeError(f"U/k sample sizes differ in {case}")
        result["k"] = [value[0] / UREF**2 for value in values]
    for field in (
        "nut", "fr1", "OmegaRot", "curvatureIndicator", "productionRatio",
        "fr1Clipped",
    ):
        try:
            components = SAMPLE_FIELD_COMPONENTS[field]
            r_f, values = sample_components(case, station, field, components)
        except RuntimeError:
            continue
        if len(r_f) == len(radii):
            result[field] = [
                value[0] if components == 1 else math.sqrt(sum(v * v for v in value))
                for value in values
            ]
    return result


def interpolate(curve: list[tuple[float, float]], x_value: float) -> float:
    if not curve[0][0] <= x_value <= curve[-1][0]:
        raise ValueError(x_value)
    for (xa, ya), (xb, yb) in zip(curve, curve[1:]):
        if xa <= x_value <= xb:
            fraction = (x_value - xa) / max(xb - xa, 1.0e-30)
            return ya + fraction * (yb - ya)
    return curve[-1][1]


def relative_errors(
    reference: Iterable[float], candidate: Iterable[float],
    characteristic_scale: float = 1.0,
) -> tuple[float, float]:
    first = list(reference)
    second = list(candidate)
    if not first or len(first) != len(second):
        raise RuntimeError("cannot compare empty or differently sized profiles")
    errors = [b - a for a, b in zip(first, second)]
    rms_scale = math.sqrt(sum(value * value for value in first) / len(first))
    max_scale = max(abs(value) for value in first)
    scale_floor = NORMALIZATION_FLOOR * abs(characteristic_scale)
    rms = math.sqrt(sum(value * value for value in errors) / len(errors))
    return (
        rms / max(rms_scale, scale_floor, 1.0e-30),
        max(map(abs, errors)) / max(max_scale, scale_floor, 1.0e-30),
    )


def weighted_metrics(
    model: list[float], reference: list[float], weights: list[float]
) -> dict[str, float]:
    total = sum(weights)
    errors = [value - target for value, target in zip(model, reference)]
    rmse = math.sqrt(sum(w * e * e for w, e in zip(weights, errors)) / total)
    ref_rms = math.sqrt(sum(w * r * r for w, r in zip(weights, reference)) / total)
    return {
        "rmse": rmse,
        "normalized_rmse": rmse / max(ref_rms, 1.0e-30),
        "bias": sum(w * e for w, e in zip(weights, errors)) / total,
        "maximum_absolute_error": max(map(abs, errors)),
    }


def reference_curves(directory: Path) -> dict[str, dict[str, list[tuple[float, float]]]]:
    specification = {
        "ux": ("zaets_exp_u.dat", 2),
        "utheta": ("zaets_exp_w.dat", 2),
        "rxr": ("zaets_exp_uv.dat", 2),
        "k": ("zaets_exp_k.dat", 2),
    }
    curves: dict[str, dict[str, list[tuple[float, float]]]] = {"n0": {}, "n06": {}}
    for quantity, (name, column) in specification.items():
        for zone, rows in parse_tecplot_zones(directory / name).items():
            rotation = "n06" if "N=0.6" in zone else "n0"
            curves[rotation][quantity] = [(row[0], row[column]) for row in rows]
    return curves


def annular_weights(radii: list[float]) -> list[float]:
    weights = []
    for index, radius in enumerate(radii):
        left = 0.0 if index == 0 else 0.5 * (radii[index - 1] + radius)
        right = 1.0 if index == len(radii) - 1 else 0.5 * (radius + radii[index + 1])
        weights.append(max(right * right - left * left, 1.0e-12))
    return weights


def write_svg(path: Path, title: str, rows: list[dict[str, float]]) -> None:
    width, height = 720, 460
    left, right, top, bottom = 70, 20, 42, 55
    xs = [row["r"] for row in rows]
    series = [("native", "baseline", "#2864b7"), ("RC", "corrected", "#c53d36"), ("reference", "reference", "#222222")]
    ys = [row[key] for _, key, _ in series for row in rows]
    ymin, ymax = min(ys), max(ys)
    padding = 0.05 * max(ymax - ymin, 1.0e-12)
    ymin -= padding
    ymax += padding

    def sx(value: float) -> float:
        return left + value * (width - left - right)

    def sy(value: float) -> float:
        return height - bottom - (value - ymin) * (height - top - bottom) / (ymax - ymin)

    content = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="25" text-anchor="middle" font-family="sans-serif">{title}</text>',
        f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#222"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#222"/>',
    ]
    for index, (label, key, colour) in enumerate(series):
        points = " ".join(f"{sx(x):.2f},{sy(row[key]):.2f}" for x, row in zip(xs, rows))
        content.append(f'<polyline points="{points}" fill="none" stroke="{colour}" stroke-width="2"/>')
        content.append(f'<text x="{100 + index*170}" y="{height-18}" fill="{colour}" font-family="sans-serif">{label}</text>')
    content.append("</svg>")
    path.write_text("\n".join(content) + "\n")


def convergence(
    case: Path, *, require_residual: bool = True
) -> tuple[bool, float | None, float | None, float | None]:
    text = (case / "log.simpleFoam").read_text()
    finite = re.search(r"(^|[^A-Za-z])(nan|inf)([^A-Za-z]|$)", text, re.I) is None
    values = re.findall(
        r"time step continuity errors\s*:\s*sum local = ([^,]+),\s*global = ([^,]+)", text
    )
    local_continuity = abs(float(values[-1][0])) if values else None
    continuity = abs(float(values[-1][1])) if values else None
    residual_matches = re.findall(
        r"Solving for (?:Ux|Uy|Uz|k|omega|nuTilda), Initial residual = ([^,]+)",
        text,
    )
    residual = max((abs(float(value)) for value in residual_matches[-12:]), default=None)
    passed = (
        "\nEnd\n" in text
        and finite
        and continuity is not None
        and continuity <= CONTINUITY_LIMIT
        and local_continuity is not None
        and local_continuity <= LOCAL_CONTINUITY_LIMIT
        and residual is not None
        and (not require_residual or residual <= RESIDUAL_LIMIT)
    )
    return passed, continuity, local_continuity, residual


def classify_performance(
    baseline: float, corrected: float, uncertainty: float
) -> str:
    """Classify a measured error change without turning it into a tuning gate."""
    delta = corrected - baseline
    if delta < -abs(uncertainty):
        return "improved"
    if delta > abs(uncertainty):
        return "degraded"
    return "neutral"


def qualify_frame_equivalence(
    corrected_rows: list[dict[str, object]],
    native_rows: list[dict[str, object]],
) -> tuple[str, str]:
    if not native_rows:
        return "not-verified", "native inertial/MRF control is absent"
    if not all(bool(row["pass"]) for row in native_rows):
        return "not-verified", "native inertial/MRF control fails the declared thresholds"
    if not corrected_rows:
        return "not-verified", "corrected inertial/MRF pair is absent"
    if all(bool(row["pass"]) for row in corrected_rows):
        return "pass", "native and corrected pairs satisfy the declared thresholds"
    return "fail", "corrected pair fails although the native frame control passes"


def analyse(results: Path, references: Path, mode: str) -> int:
    cases = sorted(path for path in results.iterdir() if (path / "log.simpleFoam").is_file())
    if not cases:
        raise RuntimeError(f"no completed cases in {results}")
    profiles = {case.name: case_profile(case) for case in cases}
    summary_rows = []
    solver_pass = True
    for case in cases:
        passed, continuity, local_continuity, residual = convergence(
            case, require_residual=mode in ("quick", "medium")
        )
        solver_pass &= passed
        summary_rows.append({
            "case": case.name,
            "completed": str(passed).lower(),
            "global_continuity": continuity,
            "local_continuity": local_continuity,
            "maximum_final_residual": residual,
        })
    with (results / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=(
            "case", "completed", "global_continuity", "local_continuity",
            "maximum_final_residual"
        ))
        writer.writeheader()
        writer.writerows(summary_rows)

    frame_rows = []
    for name, inertial in profiles.items():
        if not name.endswith("-on-inertial"):
            continue
        mrf_name = name.removesuffix("-inertial") + "-mrf"
        if mrf_name not in profiles:
            continue
        for quantity in sorted(set(inertial) & set(profiles[mrf_name]) - {"r", "OmegaRot"}):
            rms, maximum = relative_errors(inertial[quantity], profiles[mrf_name][quantity])
            frame_rows.append({"pair": f"{name}:{mrf_name}", "quantity": quantity, "relative_rms": rms, "relative_linf": maximum, "pass": rms <= FRAME_RMS_LIMIT and maximum <= FRAME_MAX_LIMIT})
    with (results / "frame-equivalence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("pair", "quantity", "relative_rms", "relative_linf", "pass"))
        writer.writeheader()
        writer.writerows(frame_rows)

    native_frame_rows = []
    for name, inertial in profiles.items():
        if not name.endswith("-native-inertial"):
            continue
        mrf_name = name.removesuffix("-inertial") + "-mrf"
        if mrf_name not in profiles:
            continue
        for quantity in sorted(
            set(inertial) & set(profiles[mrf_name]) - {"r", "OmegaRot"}
        ):
            rms, maximum = relative_errors(
                inertial[quantity], profiles[mrf_name][quantity]
            )
            native_frame_rows.append(
                {
                    "pair": f"{name}:{mrf_name}",
                    "quantity": quantity,
                    "relative_rms": rms,
                    "relative_linf": maximum,
                    "pass": rms <= FRAME_RMS_LIMIT and maximum <= FRAME_MAX_LIMIT,
                }
            )
    with (results / "native-frame-control.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "pair", "quantity", "relative_rms", "relative_linf", "pass"
            ),
        )
        writer.writeheader()
        writer.writerows(native_frame_rows)
    frame_status, frame_reason = qualify_frame_equivalence(
        frame_rows, native_frame_rows
    )

    off_rows = []
    for family in ("sst", "sa"):
        native_name = f"n06-{family}-native-mrf"
        off_name = f"n06-{family}-off-mrf"
        if native_name not in profiles or off_name not in profiles:
            continue
        for quantity in sorted(set(profiles[native_name]) & set(profiles[off_name]) - {"r"}):
            _, maximum = relative_errors(profiles[native_name][quantity], profiles[off_name][quantity])
            off_rows.append({"pair": f"{native_name}:{off_name}", "quantity": quantity, "relative_linf": maximum, "pass": maximum <= OFF_LIMIT})
    with (results / "off-equivalence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("pair", "quantity", "relative_linf", "pass"))
        writer.writeheader()
        writer.writerows(off_rows)
    off_status = "pass" if off_rows and all(row["pass"] for row in off_rows) else "not-verified"
    if mode == "quick" and off_rows and not all(row["pass"] for row in off_rows):
        off_status = "fail"

    scheme_rows = []
    for name, primary in profiles.items():
        if not name.endswith("-mrf"):
            continue
        alternate_name = name + "-upwind"
        if alternate_name not in profiles:
            continue
        for quantity in sorted(set(primary) & set(profiles[alternate_name]) - {"r", "OmegaRot"}):
            rms, maximum = relative_errors(primary[quantity], profiles[alternate_name][quantity])
            scheme_rows.append({
                "pair": f"{name}:{alternate_name}",
                "quantity": quantity,
                "relative_rms": rms,
                "relative_linf": maximum,
                "pass": rms <= SCHEME_RMS_LIMIT and maximum <= SCHEME_MAX_LIMIT,
            })
    with (results / "scheme-sensitivity.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=(
            "pair", "quantity", "relative_rms", "relative_linf", "pass"
        ))
        writer.writeheader()
        writer.writerows(scheme_rows)
    scheme_status = "not-verified"
    if scheme_rows:
        scheme_status = "pass" if all(row["pass"] for row in scheme_rows) else "fail"

    clipping_rows = []
    for name, profile in profiles.items():
        if "fr1" not in profile or "fr1Clipped" not in profile:
            continue
        count = len(profile["fr1"])
        clipped = sum(value > 0.5 for value in profile["fr1Clipped"])
        clipping_rows.append({
            "case": name,
            "sample_count": count,
            "clipped_count": clipped,
            "clipped_fraction": clipped / max(count, 1),
        })
    with (results / "clipping.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=(
            "case", "sample_count", "clipped_count", "clipped_fraction"
        ))
        writer.writeheader()
        writer.writerows(clipping_rows)

    curves = reference_curves(references)
    reference_metrics: dict[str, object] = {}
    for rotation in ("n0", "n06"):
        for family in ("sst", "sa"):
            native_name = f"{rotation}-{family}-native-inertial"
            corrected_name = f"{rotation}-{family}-on-inertial"
            if native_name not in profiles or corrected_name not in profiles:
                continue
            for quantity, curve in curves[rotation].items():
                if quantity not in profiles[native_name] or quantity not in profiles[corrected_name]:
                    continue
                radii = [r for r in profiles[native_name]["r"] if curve[0][0] <= r <= curve[-1][0]]
                indices = [profiles[native_name]["r"].index(r) for r in radii]
                baseline = [profiles[native_name][quantity][i] for i in indices]
                corrected = [profiles[corrected_name][quantity][i] for i in indices]
                target = [interpolate(curve, r) for r in radii]
                weights = annular_weights(radii)
                key = f"{rotation}-{family}-{quantity}"
                reference_metrics[key] = {
                    "baseline": weighted_metrics(baseline, target, weights),
                    "corrected": weighted_metrics(corrected, target, weights),
                }
                base_error = reference_metrics[key]["baseline"]["normalized_rmse"]
                corrected_error = reference_metrics[key]["corrected"]["normalized_rmse"]
                uncertainty = max(0.01 * base_error, NORMALIZATION_FLOOR)
                reference_metrics[key]["classification"] = classify_performance(
                    base_error, corrected_error, uncertainty
                )
                reference_metrics[key]["classification_uncertainty"] = uncertainty
                rows = [
                    {"r": r, "baseline": b, "corrected": c, "reference": t, "weight": w}
                    for r, b, c, t, w in zip(radii, baseline, corrected, target, weights)
                ]
                with (results / f"profile-{key}.csv").open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=("r", "baseline", "corrected", "reference", "weight"))
                    writer.writeheader()
                    writer.writerows(rows)
                write_svg(results / f"profile-{key}.svg", key, rows)

    criteria = [
        {"id": "solver-and-finite", "critical": True, "status": "pass" if solver_pass else "fail"},
        {"id": "frame-equivalence", "critical": True, "status": frame_status},
        {"id": "correction-off", "critical": True, "status": off_status},
        {"id": "physical-performance-classified", "critical": False, "status": "pass" if reference_metrics else "not-verified"},
        {"id": "grid-independence", "critical": True, "status": "not-verified"},
        {"id": "scheme-sensitivity", "critical": True, "status": scheme_status},
        {"id": "second-nonzero-reference", "critical": True, "status": "not-verified"},
    ]
    payload = {
        "schema_version": 2,
        "benchmark": "tmr-zaets-rotating-pipe",
        "mode": mode,
        "constants": {
            "Uref": UREF,
            "radius": RADIUS,
            "frame_rms_limit": FRAME_RMS_LIMIT,
            "frame_linf_limit": FRAME_MAX_LIMIT,
            "scheme_rms_limit": SCHEME_RMS_LIMIT,
            "scheme_linf_limit": SCHEME_MAX_LIMIT,
            "normalization_floor": NORMALIZATION_FLOOR,
        },
        "reference_metrics": reference_metrics,
        "frame_equivalence": frame_rows,
        "native_frame_control": native_frame_rows,
        "frame_qualification": {"status": frame_status, "reason": frame_reason},
        "off_equivalence": off_rows,
        "scheme_sensitivity": scheme_rows,
        "clipping": clipping_rows,
        "acceptance_criteria": criteria,
    }
    (results / "metrics.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for criterion in criteria:
        print(f"{criterion['id']}: {criterion['status']}")
    return int(any(item["status"] == "fail" for item in criteria))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--mode", choices=("smoke", "quick", "medium"), required=True)
    args = parser.parse_args()
    return analyse(args.results.resolve(), args.references.resolve(), args.mode)


if __name__ == "__main__":
    raise SystemExit(main())
