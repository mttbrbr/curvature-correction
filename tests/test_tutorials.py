#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Static consistency checks for the distributed OpenFOAM tutorials.

These checks need no OpenFOAM installation. They guard the conventions every
tutorial relies on: the runtime library is loaded, the material-derivative
convection scheme is declared, the explicit RC runtime names are used, and each
case ships run/clean scripts and a README.
"""

from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
TUTORIALS = ROOT / "tutorials/incompressible/simpleFoam"
LEGACY_NAMES = re.compile(r"\b(kOmegaSSTCC|SpalartAllmarasCC)(?:Coeffs)?\b")


def tracked(pattern: str) -> list[Path]:
    output = subprocess.run(
        ["git", "ls-files", "--", pattern],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return [ROOT / line for line in output.splitlines() if line]


def tutorial_dirs() -> list[Path]:
    return sorted(path for path in TUTORIALS.iterdir() if path.is_dir())


class TutorialConsistencyTests(unittest.TestCase):
    def test_tutorials_exist(self) -> None:
        self.assertGreater(len(tutorial_dirs()), 0)

    def test_each_tutorial_has_run_clean_scripts_and_readme(self) -> None:
        for case in tutorial_dirs():
            with self.subTest(case=case.name):
                for script in ("Allrun", "Allclean"):
                    path = case / script
                    self.assertTrue(path.is_file(), f"missing {script}")
                    self.assertTrue(os.access(path, os.X_OK), f"{script} not executable")
                self.assertTrue((case / "README.md").is_file(), "missing README.md")

    def test_control_dicts_load_the_library(self) -> None:
        dicts = tracked("tutorials/**/controlDict")
        self.assertGreater(len(dicts), 0)
        for path in dicts:
            with self.subTest(path=path.relative_to(ROOT)):
                self.assertIn("libcurvatureCorrectionModels", path.read_text())

    def test_fv_schemes_declare_material_derivative_convection(self) -> None:
        schemes = tracked("tutorials/**/fvSchemes")
        self.assertGreater(len(schemes), 0)
        for path in schemes:
            with self.subTest(path=path.relative_to(ROOT)):
                self.assertIn("div(phi,curvatureStrain)", path.read_text())

    def test_tutorials_use_explicit_rc_runtime_names(self) -> None:
        for path in tracked("tutorials"):
            with self.subTest(path=path.relative_to(ROOT)):
                try:
                    text = path.read_text()
                except UnicodeDecodeError:
                    continue
                self.assertIsNone(
                    LEGACY_NAMES.search(text),
                    "use kOmegaSST2003RC/SpalartAllmarasNoFt2RC, not legacy aliases",
                )


if __name__ == "__main__":
    unittest.main()
