#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Regression checks for repository licence and data boundaries."""

from __future__ import annotations

import json
from pathlib import Path
import re
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
ADAPTER = ROOT / "src/curvatureCorrectionModels"
MANIFEST = ROOT / "data/benchmark-sources.json"
TMR_REVISION = "ddf1840b6e928d1de60d920cafe431ed476f36a1"


class LicensingTests(unittest.TestCase):
    def test_root_notice_describes_mixed_licensing(self) -> None:
        notice = (ROOT / "LICENSE").read_text()
        self.assertIn("LICENSES/MIT.txt", notice)
        self.assertIn("LICENSES/GPL-3.0-or-later.txt", notice)
        self.assertTrue((ROOT / "LICENSES/MIT.txt").read_text().startswith("MIT License\n"))
        gpl = (ROOT / "LICENSES/GPL-3.0-or-later.txt").read_text()
        self.assertIn("GNU GENERAL PUBLIC LICENSE", gpl)
        self.assertIn("Version 3, 29 June 2007", gpl)

    def test_openfoam_adapter_points_to_gpl_and_has_spdx(self) -> None:
        self.assertIn("GPL-3.0-or-later", (ADAPTER / "LICENSE").read_text())
        sources = [
            *ADAPTER.rglob("*.C"),
            *ADAPTER.rglob("*.H"),
            ADAPTER / "Make/files",
            ADAPTER / "Make/options",
        ]
        self.assertGreaterEqual(len(sources), 10)
        for source in sources:
            with self.subTest(source=source.relative_to(ROOT)):
                self.assertIn("SPDX-License-Identifier: GPL-3.0-or-later", source.read_text())

    def test_download_manifest_is_pinned_and_complete(self) -> None:
        data = json.loads(MANIFEST.read_text())
        self.assertEqual(data["schema_version"], 1)
        self.assertEqual(data["tmr_revision"], TMR_REVISION)
        downloadable = [asset for asset in data["assets"] if asset["status"] == "downloadable"]
        held = [asset for asset in data["assets"] if asset["status"] == "hold"]
        self.assertEqual(len(downloadable), 24)
        self.assertGreaterEqual(len(held), 2)
        for asset in downloadable:
            with self.subTest(path=asset["local_path"]):
                self.assertEqual(asset["license"], "CC0-1.0")
                self.assertIn(TMR_REVISION, asset["url"] )
                self.assertRegex(asset["sha256"], r"^[0-9a-f]{64}$")
                result = subprocess.run(
                    ["git", "check-ignore", "--no-index", "-q", asset["local_path"]],
                    cwd=ROOT, check=False,
                )
                self.assertEqual(result.returncode, 0, "benchmark materialization must be ignored")
        for asset in held:
            self.assertNotIn("local_path", asset)
            self.assertNotIn("url", asset)
            self.assertIn("reason", asset)

    def test_provenance_register_records_release_boundaries(self) -> None:
        register = (ROOT / "THIRD_PARTY.md").read_text()
        for marker in (TMR_REVISION, "10.18419/DARUS-4345", "CC BY 4.0", "GPL-3.0-or-later", "user-supplied"):
            self.assertIn(marker, register)

    def test_fetcher_is_mit_and_resource_dirs_only_document_policy(self) -> None:
        self.assertIn("SPDX-License-Identifier: MIT", (ROOT / "tools/fetch_benchmark_data.py").read_text())
        for case in ("tmrFlatPlate", "smitsCurve", "tmrNaca0012"):
            readme = ROOT / f"tutorials/incompressible/simpleFoam/{case}/resources/README.md"
            self.assertIn(TMR_REVISION, readme.read_text())


if __name__ == "__main__":
    unittest.main()
