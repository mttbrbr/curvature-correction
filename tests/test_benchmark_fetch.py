#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for the checksum-verifying external benchmark fetcher."""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("benchmark_fetch", ROOT / "tools/fetch_benchmark_data.py")
assert SPEC and SPEC.loader
FETCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FETCH)


class BenchmarkFetchTests(unittest.TestCase):
    def asset(self, source: Path, digest: str) -> dict:
        return {
            "case": "demo",
            "license": "CC0-1.0",
            "local_path": "tutorials/demo/resources/input.dat",
            "sha256": digest,
            "status": "downloadable",
            "upstream_path": "demo/input.dat",
            "url": source.as_uri(),
        }

    def test_fetch_verify_and_materialize(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "upstream.dat"
            source.write_bytes(b"pinned benchmark\n")
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            manifest = {"schema_version": 1, "assets": [self.asset(source, digest)]}
            cache = root / "cache"
            repository = root / "repository"
            ready = FETCH.fetch_assets(
                manifest, cases={"demo"}, cache_dir=cache, repository_root=repository,
                verify_only=False, create_symlinks=True,
            )
            self.assertEqual(len(ready), 1)
            materialized = repository / "tutorials/demo/resources/input.dat"
            self.assertTrue(materialized.is_symlink())
            self.assertEqual(materialized.read_bytes(), source.read_bytes())
            FETCH.fetch_assets(
                manifest, cases={"demo"}, cache_dir=cache, repository_root=repository,
                verify_only=True, create_symlinks=True,
            )

    def test_checksum_mismatch_is_rejected_and_removed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "upstream.dat"
            source.write_bytes(b"unexpected")
            manifest = {"schema_version": 1, "assets": [self.asset(source, "0" * 64)]}
            with self.assertRaises(FETCH.FetchError):
                FETCH.fetch_assets(
                    manifest, cases={"demo"}, cache_dir=root / "cache",
                    repository_root=root / "repository", verify_only=False,
                    create_symlinks=False,
                )
            self.assertFalse((root / "cache/tutorials/demo/resources/input.dat").exists())

    def test_manifest_paths_cannot_escape_repository(self) -> None:
        for path in ("/absolute", "../escape", "safe/../../escape"):
            with self.subTest(path=path), self.assertRaises(FETCH.FetchError):
                FETCH.safe_relative_path(path)


if __name__ == "__main__":
    unittest.main()
