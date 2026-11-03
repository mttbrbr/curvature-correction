#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Analyse the fully developed Xiao et al. rotating-pipe validation."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re

RADIUS = 0.5
DIAMETER = 1.0
UBULK = 1.0
LEVELS = ("smoke", "quick", "medium")
CELL_COUNTS = {"smoke": 7488, "quick": 59904, "medium": 479232}
ROTATIONS = {"n0": 0.0, "n05": 0.5, "n1": 1.0}
OFF_LIMIT = 1.0e-10
SCHEME_LIMIT = 0.04
GRID_GCI_LIMIT = 0.02
CONTINUITY_LIMIT = 1.0e-5
LOCAL_CONTINUITY_LIMIT = 1.0e-4
RESIDUAL_LIMIT = 1.0e-5
COMPONENTS = {
    "U": 3, "turbulenceProperties:R": 6, "nut": 1, "k": 1,
    "fr1": 1, "fr1Clipped": 1, "curvatureIndicator": 1,
    "productionRatio": 1,
}


def numeric_rows(path: Path) -> list[list[float]]:
    output = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            row = [float(value) for value in line.split()]
        except ValueError:
            continue
        if not all(math.isfinite(value) for value in row):
            raise RuntimeError(f"non-finite sample in {path}")
        output.append(row)
    if not output:
        raise RuntimeError(f"no numeric samples in {path}")
    return output


def latest_directory(root: Path) -> Path:
    candidates = [(float(path.name), path) for path in root.iterdir() if path.is_dir() and _is_float(path.name)]
    if not candidates:
        raise RuntimeError(f"no time directory below {root}")
    return max(candidates)[1]


def _is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def sample_components(case: Path, field: str) -> tuple[list[float], list[list[float]]]:
    directory = latest_directory(case / "postProcessing/sampleProfiles")
    matches = [path for path in directory.iterdir() if path.name.startswith("radialProfile_") and field in path.name]
    if len(matches) != 1:
        raise RuntimeError(f"expected one {field} sample in {directory}, found {matches}")
    path = matches[0]
    names = path.name.removesuffix(path.suffix).removeprefix("radialProfile_").split("_")
    widths = [COMPONENTS[name] for name in names]
    index = names.index(field)
    rows = numeric_rows(path)
    coordinates = len(rows[0]) - sum(widths)
    start = coordinates + sum(widths[:index])
    radii, values = [], []
    for row in rows:
        radius = row[0] + 1.0e-6 if coordinates == 1 else math.hypot(row[1], row[2])
        radii.append(radius / RADIUS)
        values.append(row[start:start + widths[index]])
    return radii, values


def profile(case: Path) -> dict[str, list[float]]:
    radii, velocity = sample_components(case, "U")
    result = {
        "r": radii,
        "ux": [value[0] / UBULK for value in velocity],
        "utheta": [value[2] / UBULK for value in velocity],
    }
    for field in ("nut", "k", "fr1", "fr1Clipped", "curvatureIndicator", "productionRatio"):
        try:
            _, values = sample_components(case, field)
        except (RuntimeError, KeyError, ValueError):
            continue
        result[field] = [value[0] for value in values]
    return result


def relative_error(reference: list[float], candidate: list[float], floor: float = 1.0e-3) -> tuple[float, float]:
    if not reference or len(reference) != len(candidate):
        raise RuntimeError("incompatible profiles")
    errors = [b - a for a, b in zip(reference, candidate)]
    scale = max(math.sqrt(sum(value * value for value in reference) / len(reference)), floor)
    return math.sqrt(sum(value * value for value in errors) / len(errors)) / scale, max(map(abs, errors)) / scale


def reference_data(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())["references"]["xiao-2024-rotating-pipe"]


def fit_values(rotation: float, radii: list[float]) -> tuple[list[float], list[float]]:
    phi = 0.90 + 0.25 * math.log(rotation)
    alpha = 2.0 - 0.071 * rotation ** 1.2
    psi = 2.48 - 0.45 * rotation ** 0.57
    return ([phi / rotation * radius ** alpha for radius in radii], [radius ** psi for radius in radii])


def convergence(case: Path, residual_limit: float = RESIDUAL_LIMIT) -> tuple[bool, float | None, float | None, float | None, float | None]:
    text = (case / "log.simpleFoam").read_text()
    finite = re.search(r"(^|[^A-Za-z])(nan|inf)([^A-Za-z]|$)", text, re.I) is None
    continuity_values = re.findall(r"time step continuity errors\s*:\s*sum local = ([^,]+),\s*global = ([^,]+)", text)
    residual_values = re.findall(r"Solving for (?:Ux|Uy|Uz|k|omega|nuTilda), Initial residual = ([^,]+)", text)
    gradient_values = re.findall(r"pressure gradient\s*=\s*([-+0-9.eE]+)", text, re.I)
    local_continuity = abs(float(continuity_values[-1][0])) if continuity_values else None
    continuity = abs(float(continuity_values[-1][1])) if continuity_values else None
    residual = max((abs(float(value)) for value in residual_values[-12:]), default=None)
    gradient = abs(float(gradient_values[-1])) if gradient_values else None
    passed = bool(
        "\nEnd\n" in text and finite and continuity is not None
        and continuity <= CONTINUITY_LIMIT and local_continuity is not None
        and local_continuity <= LOCAL_CONTINUITY_LIMIT and residual is not None
        and residual <= residual_limit and gradient is not None
    )
    return passed, continuity, local_continuity, residual, gradient


def classify(baseline: float, corrected: float, uncertainty: float) -> str:
    if corrected < baseline - abs(uncertainty):
        return "improved"
    if corrected > baseline + abs(uncertainty):
        return "degraded"
    return "neutral"


def analyse_level(results: Path, reference_path: Path, level: str, mode: str) -> dict[str, object]:
    references = reference_data(reference_path)
    globals_by_n = {float(row["N"]): row for row in references["global_values"]}
    cases = sorted(path for path in results.iterdir() if (path / "log.simpleFoam").is_file())
    profiles = {case.name: profile(case) for case in cases}
    summary, friction = [], {}
    solver_pass = True
    for case in cases:
        residual_limit = 1.0e-3 if mode == "smoke" else RESIDUAL_LIMIT
        passed, continuity, local_continuity, residual, gradient = convergence(case, residual_limit)
        solver_pass &= passed
        factor = 2.0 * DIAMETER * gradient / UBULK ** 2 if gradient is not None else None
        friction[case.name] = factor
        summary.append({"case": case.name, "completed": passed, "global_continuity": continuity, "local_continuity": local_continuity, "maximum_final_residual": residual, "friction_factor": factor})
    with (results / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(summary[0]))
        writer.writeheader(); writer.writerows(summary)

    comparisons = {}
    for name, values in profiles.items():
        tokens = name.split("-")
        if len(tokens) != 3 or tokens[2] == "off":
            continue
        rotation, family, variant = tokens
        n_value = ROTATIONS[rotation]
        target = globals_by_n[n_value]
        friction_error = abs(friction[name] - target["friction_factor"]) / target["friction_factor"]
        re_tau = 44000.0 * math.sqrt(friction[name] / 32.0)
        record: dict[str, object] = {
            "friction_factor_relative_error": friction_error,
            "re_tau": re_tau,
            "re_tau_relative_error": abs(re_tau - target["re_tau"]) / target["re_tau"],
        }
        if n_value > 0:
            nonrotating = f"n0-{family}-{variant}"
            if nonrotating in friction:
                drag_reduction = 1.0 - friction[name] / friction[nonrotating]
                record.update({
                    "drag_reduction": drag_reduction,
                    "drag_reduction_absolute_error": abs(drag_reduction - target["drag_reduction"]),
                })
            axial_indices = [index for index, radius in enumerate(values["r"]) if radius <= 0.9]
            theta_indices = [index for index, radius in enumerate(values["r"]) if 0.2 <= radius <= 0.9]
            ucl = values["ux"][0]
            axial = [(ucl - values["ux"][index]) / n_value for index in axial_indices]
            theta = [values["utheta"][index] / n_value for index in theta_indices]
            axial_target, _ = fit_values(n_value, [values["r"][index] for index in axial_indices])
            _, theta_target = fit_values(n_value, [values["r"][index] for index in theta_indices])
            axial_error = relative_error(axial_target, axial)[0]
            theta_error = relative_error(theta_target, theta)[0]
            record.update({"axial_fit_normalized_rmse": axial_error, "azimuthal_fit_normalized_rmse": theta_error})
        comparisons[name] = record

    classifications = []
    for rotation in ROTATIONS:
        for family in ("sst", "sa"):
            native, corrected = f"{rotation}-{family}-native", f"{rotation}-{family}-on"
            if native not in comparisons or corrected not in comparisons:
                continue
            def score(name: str) -> float:
                errors = [
                    float(value) for key, value in comparisons[name].items()
                    if key.endswith("error") or key.endswith("rmse")
                ]
                return sum(errors) / len(errors)
            uncertainty = max(float(globals_by_n[ROTATIONS[rotation]]["friction_factor_relative_uncertainty"]), 0.005)
            classifications.append({"rotation": rotation, "family": family, "baseline_score": score(native), "corrected_score": score(corrected), "uncertainty": uncertainty, "classification": classify(score(native), score(corrected), uncertainty)})

    off_rows = []
    for family in ("sst", "sa"):
        native, off = f"n1-{family}-native", f"n1-{family}-off"
        if native not in profiles or off not in profiles:
            continue
        rms, maximum = relative_error(profiles[native]["ux"], profiles[off]["ux"])
        friction_delta = abs(friction[native] - friction[off]) / max(abs(friction[native]), 1.0e-30)
        off_rows.append({"family": family, "profile_rms": rms, "profile_linf": maximum, "friction_relative_error": friction_delta, "pass": max(rms, maximum, friction_delta) <= OFF_LIMIT})

    scheme_rows = []
    for family in ("sst", "sa"):
        for variant in ("native", "on"):
            primary, alternate = f"n1-{family}-{variant}", f"n1-{family}-{variant}-upwind"
            if primary not in profiles or alternate not in profiles:
                continue
            rms, maximum = relative_error(profiles[primary]["ux"], profiles[alternate]["ux"])
            friction_delta = abs(friction[primary] - friction[alternate]) / max(abs(friction[primary]), 1.0e-30)
            scheme_rows.append({"family": family, "variant": variant, "profile_rms": rms, "profile_linf": maximum, "friction_relative_error": friction_delta, "pass": max(rms, maximum, friction_delta) <= SCHEME_LIMIT})

    clipping = []
    for name, values in profiles.items():
        if "fr1Clipped" not in values:
            continue
        count = len(values["fr1Clipped"])
        clipped = sum(value > 0.5 for value in values["fr1Clipped"])
        clipping.append({"case": name, "sample_count": count, "clipped_count": clipped, "clipped_fraction": clipped / max(count, 1)})

    criteria = [
        {"id": "solver-and-finite", "critical": True, "status": "pass" if solver_pass else "fail"},
        {"id": "reference-coverage", "critical": True, "status": "pass" if len(classifications) == 6 else "not-verified"},
        {"id": "correction-off", "critical": True, "status": "pass" if off_rows and all(row["pass"] for row in off_rows) else "not-verified"},
        {"id": "scheme-sensitivity", "critical": True, "status": "pass" if scheme_rows and all(row["pass"] for row in scheme_rows) else "not-verified"},
        {"id": "physical-performance-classified", "critical": False, "status": "pass" if classifications else "not-verified"},
    ]
    payload = {"schema_version": 2, "benchmark": "xiao-2024-rotating-pipe", "mode": mode, "level": level, "cell_count": CELL_COUNTS[level], "friction_factors": friction, "reference_comparisons": comparisons, "physical_classification": classifications, "off_equivalence": off_rows, "scheme_sensitivity": scheme_rows, "clipping": clipping, "acceptance_criteria": criteria}
    (results / "metrics.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def grid_convergence(values: list[float]) -> dict[str, object]:
    coarse_difference = abs(values[0] - values[1])
    fine_difference = abs(values[1] - values[2])
    contracts = fine_difference < coarse_difference or coarse_difference <= 1.0e-14
    if coarse_difference <= 1.0e-14 or fine_difference <= 1.0e-14:
        order, gci = None, 0.0
    elif not contracts or values[2] == 0:
        order, gci = None, None
    else:
        order = math.log(coarse_difference / fine_difference) / math.log(2.0)
        gci = 1.25 * fine_difference / abs(values[2]) / (2.0 ** order - 1.0)
    return {"values": dict(zip(LEVELS, values)), "coarse_medium_difference": coarse_difference, "medium_fine_difference": fine_difference, "contracts": contracts, "observed_order": order, "fine_gci": gci, "gci_limit": GRID_GCI_LIMIT, "pass": contracts and gci is not None and gci <= GRID_GCI_LIMIT}


def analyse_validation(root: Path) -> dict[str, object]:
    payloads = {level: json.loads((root / level / "metrics.json").read_text()) for level in LEVELS if (root / level / "metrics.json").is_file()}
    grid = {}
    if len(payloads) == 3:
        common = set(payloads["smoke"]["friction_factors"])
        common &= set(payloads["quick"]["friction_factors"])
        common &= set(payloads["medium"]["friction_factors"])
        common = {case for case in common if len(case.split("-")) == 3 and not case.endswith("-off")}
        for case in sorted(common):
            grid[case] = grid_convergence([float(payloads[level]["friction_factors"][case]) for level in LEVELS])
    coverage = len(payloads) == 3 and all(len(payloads[level]["physical_classification"]) == 6 for level in LEVELS)
    solver_status = "pass" if len(payloads) == 3 and all(
        next(item for item in payloads[level]["acceptance_criteria"] if item["id"] == "solver-and-finite")["status"] == "pass"
        for level in LEVELS
    ) else "fail" if len(payloads) == 3 else "not-verified"
    grid_status = "pass" if grid and all(item["pass"] for item in grid.values()) else ("fail" if len(payloads) == 3 else "not-verified")
    medium_criteria = {item["id"]: item["status"] for item in payloads.get("medium", {}).get("acceptance_criteria", [])}
    criteria = [
        {"id": "solver-and-finite", "critical": True, "status": solver_status},
        {"id": "three-grid-reference-coverage", "critical": True, "status": "pass" if coverage else "not-verified"},
        {"id": "grid-independence", "critical": True, "status": grid_status},
        {"id": "correction-off", "critical": True, "status": medium_criteria.get("correction-off", "not-verified")},
        {"id": "scheme-sensitivity", "critical": True, "status": medium_criteria.get("scheme-sensitivity", "not-verified")},
        {"id": "physical-performance-classified", "critical": False, "status": "pass" if coverage else "not-verified"},
    ]
    result = {"schema_version": 2, "benchmark": "xiao-2024-rotating-pipe-validation", "levels": list(payloads), "grid_convergence": grid, "physical_classification": payloads.get("medium", {}).get("physical_classification", []), "acceptance_criteria": criteria}
    (root / "metrics.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path)
    parser.add_argument("--validation-root", type=Path)
    parser.add_argument("--reference-models", type=Path, required=True)
    parser.add_argument("--level", choices=LEVELS)
    parser.add_argument("--mode", choices=("smoke", "validation"))
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    if args.validation_root:
        payload = analyse_validation(args.validation_root.resolve())
    elif args.results and args.level and args.mode:
        payload = analyse_level(args.results.resolve(), args.reference_models.resolve(), args.level, args.mode)
    else:
        parser.error("use --validation-root or --results with --level and --mode")
    for criterion in payload["acceptance_criteria"]:
        print(f"{criterion['id']}: {criterion['status']}")
    failed = any(item["critical"] and item["status"] == "fail" for item in payload["acceptance_criteria"])
    incomplete = args.require_complete and any(item["critical"] and item["status"] != "pass" for item in payload["acceptance_criteria"])
    return int(failed or incomplete)


if __name__ == "__main__":
    raise SystemExit(main())
