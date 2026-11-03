#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Unit tests for deterministic rotating-pipe post-processing."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "tutorials/incompressible/simpleFoam/tmrRotatingPipe/tools/analyse.py"
)
SPEC = importlib.util.spec_from_file_location("rotating_pipe_analysis", SCRIPT)
assert SPEC and SPEC.loader
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)
MESH_SCRIPT = (
    ROOT
    / "tutorials/incompressible/simpleFoam/tmrRotatingPipe/tools/generate_mesh.py"
)
MESH_SPEC = importlib.util.spec_from_file_location(
    "rotating_pipe_mesh", MESH_SCRIPT
)
assert MESH_SPEC and MESH_SPEC.loader
MESH = importlib.util.module_from_spec(MESH_SPEC)
MESH_SPEC.loader.exec_module(MESH)
RESOLUTION_SCRIPT = (
    ROOT
    / "tutorials/incompressible/simpleFoam/tmrRotatingPipe/tools/analyse_frame_resolution.py"
)
RESOLUTION_SPEC = importlib.util.spec_from_file_location(
    "rotating_pipe_frame_resolution", RESOLUTION_SCRIPT
)
assert RESOLUTION_SPEC and RESOLUTION_SPEC.loader
RESOLUTION = importlib.util.module_from_spec(RESOLUTION_SPEC)
RESOLUTION_SPEC.loader.exec_module(RESOLUTION)


class RotatingPipeAnalysisTests(unittest.TestCase):
    def test_frame_resolution_aggregation_and_missing_modes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            runs = Path(directory)
            for mode, rms in (("smoke", 0.2), ("quick", 0.1)):
                target = runs / mode
                target.mkdir()
                (target / "metrics.json").write_text(
                    json.dumps(
                        {
                            "native_frame_control": [
                                {
                                    "pair": "native",
                                    "quantity": "ux",
                                    "relative_rms": rms,
                                    "relative_linf": 2 * rms,
                                    "pass": False,
                                }
                            ],
                            "frame_equivalence": [],
                        }
                    )
                )
            payload = RESOLUTION.analyse(runs)
            self.assertEqual(payload["status"], "not-verified")
            self.assertEqual(payload["missing_modes"], ["medium"])
            self.assertEqual(len(payload["contractions"]), 1)
            self.assertAlmostEqual(
                payload["contractions"][0]["rms_contraction"], 0.5
            )
            self.assertTrue(payload["contractions"][0]["contracts"])
            self.assertTrue((runs / "frame-resolution.csv").is_file())

    def test_frame_qualification_requires_a_valid_native_control(self) -> None:
        passing = [{"pass": True}]
        failing = [{"pass": False}]
        status, _ = ANALYSIS.qualify_frame_equivalence(passing, [])
        self.assertEqual(status, "not-verified")
        status, _ = ANALYSIS.qualify_frame_equivalence(passing, failing)
        self.assertEqual(status, "not-verified")
        status, _ = ANALYSIS.qualify_frame_equivalence([], passing)
        self.assertEqual(status, "not-verified")
        status, _ = ANALYSIS.qualify_frame_equivalence(passing, passing)
        self.assertEqual(status, "pass")
        status, _ = ANALYSIS.qualify_frame_equivalence(failing, passing)
        self.assertEqual(status, "fail")
        status, _ = ANALYSIS.qualify_frame_equivalence(failing, failing)
        self.assertEqual(status, "not-verified")
        status, _ = ANALYSIS.qualify_frame_equivalence(failing, [])
        self.assertEqual(status, "not-verified")

    def test_performance_classification_uses_uncertainty_band(self) -> None:
        self.assertEqual(ANALYSIS.classify_performance(0.20, 0.15, 0.01), "improved")
        self.assertEqual(ANALYSIS.classify_performance(0.20, 0.205, 0.01), "neutral")
        self.assertEqual(ANALYSIS.classify_performance(0.20, 0.25, 0.01), "degraded")

    def test_three_grid_gci_requires_contraction(self) -> None:
        rows = []
        for mode, rms in (("smoke", 0.04), ("quick", 0.02), ("medium", 0.01)):
            rows.append({
                "mode": mode,
                "kind": "native",
                "pair": "n06-sst-native-inertial:n06-sst-native-mrf",
                "quantity": "ux",
                "relative_rms": rms,
                "relative_linf": 2 * rms,
            })
        gci = RESOLUTION.grid_convergence_rows(rows)
        self.assertEqual(len(gci), 1)
        self.assertTrue(gci[0]["rms_contracts"])
        self.assertGreater(gci[0]["rms_observed_order"], 0)
        self.assertIsNotNone(gci[0]["rms_fine_gci"])

    def test_aggregated_openfoam_sample_columns(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            sample_dir = case / "postProcessing/sampleProfiles/300"
            sample_dir.mkdir(parents=True)
            sample = sample_dir / (
                "station_k_nut_OmegaRot_U_turbulenceProperties:R.xy"
            )
            sample.write_text(
                "0.0 0.25 0.001 1 2 3 8 9 10 11 12 13 14 15 16\n"
            )
            radii, velocity = ANALYSIS.sample_components(
                case, "station", "U", 3
            )
            self.assertAlmostEqual(radii[0], 1.0e-6 / ANALYSIS.RADIUS)
            self.assertEqual(velocity, [[8.0, 9.0, 10.0]])
            _, kinetic_energy = ANALYSIS.sample_components(
                case, "station", "k", 1
            )
            self.assertEqual(kinetic_energy, [[0.25]])

    def test_ogrid_generator_counts_and_zones(self) -> None:
        self.assertEqual(MESH.expected_cells("smoke"), 98163)
        self.assertEqual(MESH.expected_cells("quick"), 238160)
        self.assertEqual(MESH.expected_cells("medium"), 854568)
        quick = MESH.mesh_text("quick")
        self.assertEqual(quick.count("    hex ("), 20)
        self.assertEqual(quick.count(") rotating ("), 15)
        for patch in ("inlet", "outlet", "stationaryWall", "rotatingWall"):
            self.assertIn(f"    {patch}\n", quick)
        self.assertNotIn("type wedge", quick)

    def test_tecplot_zones_and_bulk_ratio(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "profile.dat"
            source.write_text(
                '# header\n'
                'ZONE T="Exp, N=0"\n'
                '0.0 1.0 1.0\n'
                '0.9 0.7 0.7\n'
                'ZONE T="Exp, N=0.6"\n'
                '0.0 1.0 1.0\n'
            )
            zones = ANALYSIS.parse_tecplot_zones(source)
            self.assertEqual(set(zones), {"Exp, N=0", "Exp, N=0.6"})
            ratio = ANALYSIS.bulk_ratio_from_reference(zones["Exp, N=0"])
            self.assertGreater(ratio, 0.6)
            self.assertLess(ratio, 1.0)

    def test_relative_errors_identical_and_scaled(self) -> None:
        self.assertEqual(
            ANALYSIS.relative_errors([1.0, 2.0], [1.0, 2.0]),
            (0.0, 0.0),
        )
        rms, maximum = ANALYSIS.relative_errors([1.0, 1.0], [1.01, 0.99])
        self.assertAlmostEqual(rms, 0.01)
        self.assertAlmostEqual(maximum, 0.01)

    def test_interpolation_and_annular_weights(self) -> None:
        curve = [(0.0, 0.0), (1.0, 2.0)]
        self.assertEqual(ANALYSIS.interpolate(curve, 0.25), 0.5)
        weights = ANALYSIS.annular_weights([0.0, 0.5, 1.0])
        self.assertTrue(all(weight > 0 for weight in weights))
        self.assertAlmostEqual(sum(weights), 1.0)

    def test_nonfinite_sample_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "sample.raw"
            source.write_text("0 1 nan\n")
            with self.assertRaises(RuntimeError):
                ANALYSIS.numeric_rows(source)


if __name__ == "__main__":
    unittest.main()
