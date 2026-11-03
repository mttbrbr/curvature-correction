#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Deterministic tests for the independent rotating-pipe DNS branch."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = ROOT / "tutorials/incompressible/simpleFoam/xiaoRotatingPipe"


def load(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    assert specification and specification.loader
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


MESH = load("xiao_mesh", TUTORIAL / "tools/generate_mesh.py")
ANALYSIS = load("xiao_analysis", TUTORIAL / "tools/analyse.py")


class XiaoRotatingPipeTests(unittest.TestCase):
    def test_nested_periodic_ogrid(self) -> None:
        self.assertEqual(MESH.expected_cells("smoke"), 7488)
        self.assertEqual(MESH.expected_cells("quick"), 59904)
        self.assertEqual(MESH.expected_cells("medium"), 479232)
        mesh = MESH.mesh_text("quick")
        self.assertEqual(mesh.count("    hex ("), 5)
        self.assertEqual(mesh.count(") rotating ("), 5)
        self.assertIn("neighbourPatch periodicOut", mesh)
        self.assertIn("neighbourPatch periodicIn", mesh)
        self.assertNotIn("type wedge", mesh)

    def test_published_profile_fits(self) -> None:
        axial, azimuthal = ANALYSIS.fit_values(1.0, [0.0, 0.5, 0.9])
        self.assertEqual(axial[0], 0.0)
        self.assertEqual(azimuthal[0], 0.0)
        self.assertGreater(axial[-1], axial[1])
        self.assertGreater(azimuthal[-1], azimuthal[1])

    def test_classification_has_neutral_uncertainty_band(self) -> None:
        self.assertEqual(ANALYSIS.classify(0.2, 0.1, 0.01), "improved")
        self.assertEqual(ANALYSIS.classify(0.2, 0.205, 0.01), "neutral")
        self.assertEqual(ANALYSIS.classify(0.2, 0.3, 0.01), "degraded")

    def test_grid_convergence_uses_solution_differences(self) -> None:
        result = ANALYSIS.grid_convergence([0.030, 0.021, 0.020])
        self.assertTrue(result["contracts"])
        self.assertTrue(result["pass"])
        self.assertGreater(result["observed_order"], 3.0)
        self.assertGreater(result["fine_gci"], 0)
        divergent = ANALYSIS.grid_convergence([0.021, 0.022, 0.025])
        self.assertFalse(divergent["pass"])


if __name__ == "__main__":
    unittest.main()
