#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Aggregate inertial/MRF errors across rotating-pipe mesh levels."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


MODES = ("smoke", "quick", "medium")
CELL_COUNTS = {"smoke": 98163, "quick": 238160, "medium": 854568}
TABLES = (
    ("native", "native_frame_control"),
    ("corrected", "frame_equivalence"),
)
FRAME_RMS_LIMIT = 0.01
FRAME_LINF_LIMIT = 0.02
SAFETY_FACTOR = 1.25


def collect_rows(runs: Path) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    present = []
    for mode in MODES:
        path = runs / mode / "metrics.json"
        if not path.is_file():
            continue
        present.append(mode)
        payload = json.loads(path.read_text())
        for kind, key in TABLES:
            for source in payload.get(key, []):
                rows.append(
                    {
                        "mode": mode,
                        "cells": CELL_COUNTS[mode],
                        "kind": kind,
                        "pair": source["pair"],
                        "quantity": source["quantity"],
                        "relative_rms": float(source["relative_rms"]),
                        "relative_linf": float(source["relative_linf"]),
                        "threshold_pass": bool(source["pass"]),
                    }
                )
    return rows, present


def contraction_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    indexed = {
        (str(row["kind"]), str(row["pair"]), str(row["quantity"]), str(row["mode"])): row
        for row in rows
    }
    output = []
    identities = sorted(
        {(str(row["kind"]), str(row["pair"]), str(row["quantity"])) for row in rows}
    )
    for kind, pair, quantity in identities:
        for coarse, fine in zip(MODES, MODES[1:]):
            first = indexed.get((kind, pair, quantity, coarse))
            second = indexed.get((kind, pair, quantity, fine))
            if first is None or second is None:
                continue
            coarse_rms = float(first["relative_rms"])
            fine_rms = float(second["relative_rms"])
            output.append(
                {
                    "kind": kind,
                    "pair": pair,
                    "quantity": quantity,
                    "coarse": coarse,
                    "fine": fine,
                    "rms_contraction": fine_rms / max(coarse_rms, 1.0e-30),
                    "contracts": fine_rms < coarse_rms,
                }
            )
    return output


def grid_convergence_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Compute an observed order and fine-grid GCI for frame discrepancies."""
    indexed = {
        (str(row["kind"]), str(row["pair"]), str(row["quantity"]), str(row["mode"])): row
        for row in rows
    }
    identities = sorted(
        {(str(row["kind"]), str(row["pair"]), str(row["quantity"])) for row in rows}
    )
    output = []
    for kind, pair, quantity in identities:
        samples = [indexed.get((kind, pair, quantity, mode)) for mode in MODES]
        if any(sample is None for sample in samples):
            continue
        values_by_metric = {
            "rms": [float(sample["relative_rms"]) for sample in samples],
            "linf": [float(sample["relative_linf"]) for sample in samples],
        }
        record: dict[str, object] = {"kind": kind, "pair": pair, "quantity": quantity}
        accepted = True
        for metric, values in values_by_metric.items():
            monotonic = values[2] < values[1] < values[0] or max(values) <= 1.0e-14
            h = [CELL_COUNTS[mode] ** (-1.0 / 3.0) for mode in MODES]
            orders = []
            for coarse_error, fine_error, coarse_h, fine_h in zip(
                values, values[1:], h, h[1:]
            ):
                if coarse_error > fine_error > 0:
                    orders.append(math.log(coarse_error / fine_error) / math.log(coarse_h / fine_h))
            order = sum(orders) / len(orders) if orders else None
            if max(values) <= 1.0e-14:
                gci = 0.0
            elif order is None or order <= 0:
                gci = None
            else:
                refinement = h[1] / h[2]
                gci = SAFETY_FACTOR * values[2] / max(refinement ** order - 1.0, 1.0e-30)
            floor = FRAME_RMS_LIMIT if metric == "rms" else FRAME_LINF_LIMIT
            uncertainty = max(floor, gci or 0.0)
            metric_pass = monotonic and gci is not None and values[2] <= uncertainty
            record.update({
                f"{metric}_coarse": values[0],
                f"{metric}_medium": values[1],
                f"{metric}_fine": values[2],
                f"{metric}_observed_order": order,
                f"{metric}_fine_gci": gci,
                f"{metric}_acceptance_limit": uncertainty,
                f"{metric}_contracts": monotonic,
                f"{metric}_pass": metric_pass,
            })
            accepted &= metric_pass
        record["pass"] = accepted
        output.append(record)
    return output


def analyse(runs: Path) -> dict[str, object]:
    rows, present = collect_rows(runs)
    contractions = contraction_rows(rows)
    gci = grid_convergence_rows(rows)
    missing = [mode for mode in MODES if mode not in present]
    expected_controls = all(
        any(row["kind"] == kind and f"-{family}-" in str(row["pair"]) for row in gci)
        for kind in ("native", "corrected")
        for family in ("sst", "sa")
    )
    mode_payloads = {
        mode: json.loads((runs / mode / "metrics.json").read_text())
        for mode in present
    }
    criterion_maps = {
        mode: {
            item["id"]: item["status"]
            for item in payload.get("acceptance_criteria", [])
        }
        for mode, payload in mode_payloads.items()
    }
    solver_status = (
        "pass" if not missing and all(
            criteria.get("solver-and-finite") == "pass"
            for criteria in criterion_maps.values()
        ) else "fail" if any(
            criteria.get("solver-and-finite") == "fail"
            for criteria in criterion_maps.values()
        ) else "not-verified"
    )

    def any_level_status(identifier: str) -> str:
        statuses = [criteria.get(identifier) for criteria in criterion_maps.values()]
        if "pass" in statuses:
            return "pass"
        if "fail" in statuses:
            return "fail"
        return "not-verified"

    off_status = any_level_status("correction-off")
    scheme_status = any_level_status("scheme-sensitivity")
    if missing or not expected_controls or not gci:
        frame_status = "not-verified"
    else:
        frame_status = "pass" if all(bool(row["pass"]) for row in gci) else "fail"
    criteria = [
        {"id": "solver-and-finite", "critical": True, "status": solver_status},
        {"id": "three-grid-frame-equivalence", "critical": True, "status": frame_status},
        {"id": "correction-off", "critical": True, "status": off_status},
        {"id": "scheme-sensitivity", "critical": True, "status": scheme_status},
    ]
    if any(item["status"] == "fail" for item in criteria):
        status = "fail"
    elif any(item["status"] != "pass" for item in criteria):
        status = "not-verified"
    else:
        status = "pass"
    payload: dict[str, object] = {
        "schema_version": 2,
        "benchmark": "tmr-zaets-frame-resolution",
        "status": status,
        "present_modes": present,
        "missing_modes": missing,
        "expected_controls_present": expected_controls,
        "cell_counts": CELL_COUNTS,
        "frame_metrics": rows,
        "contractions": contractions,
        "grid_convergence": gci,
        "acceptance_criteria": criteria,
    }
    with (runs / "frame-resolution.csv").open("w", newline="") as stream:
        fields = (
            "mode", "cells", "kind", "pair", "quantity",
            "relative_rms", "relative_linf", "threshold_pass",
        )
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    (runs / "frame-resolution.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    with (runs / "frame-gci.csv").open("w", newline="") as stream:
        if gci:
            writer = csv.DictWriter(stream, fieldnames=tuple(gci[0]))
            writer.writeheader()
            writer.writerows(gci)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    payload = analyse(args.runs.resolve())
    print(
        "frame-resolution: " + str(payload["status"])
        + "; present=" + ",".join(payload["present_modes"])
        + "; missing=" + ",".join(payload["missing_modes"])
    )
    return int(args.require_complete and payload["status"] != "pass")


if __name__ == "__main__":
    raise SystemExit(main())
