#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 mttbrbr
"""Download pinned benchmark inputs into a Git-ignored cache.

The manifest records upstream URLs, licences and SHA-256 digests. Downloadable
assets can be materialised as ignored symlinks at the paths expected by the
OpenFOAM tutorials. Assets marked ``hold`` are never downloaded automatically.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data/benchmark-sources.json"
DEFAULT_CACHE = ROOT / ".cache/curvature-benchmarks"


class FetchError(RuntimeError):
    """Raised when an asset cannot be fetched or verified safely."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema_version") != 1 or not isinstance(data.get("assets"), list):
        raise FetchError("benchmark manifest must use schema_version 1 and contain assets")
    return data


def safe_relative_path(value: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise FetchError(f"unsafe manifest path: {value}")
    return path


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "curvature-correction-benchmark-fetcher"})
    temporary_name: str | None = None
    try:
        with urllib.request.urlopen(request) as response, tempfile.NamedTemporaryFile(
            dir=destination.parent, prefix=destination.name + ".", delete=False
        ) as temporary:
            temporary_name = temporary.name
            while chunk := response.read(1024 * 1024):
                temporary.write(chunk)
        os.replace(temporary_name, destination)
        temporary_name = None
    except OSError as error:
        raise FetchError(f"cannot download {url}: {error}") from error
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


def materialize(cache_path: Path, repository_path: Path) -> None:
    repository_path.parent.mkdir(parents=True, exist_ok=True)
    if repository_path.is_symlink():
        if repository_path.resolve() == cache_path.resolve():
            return
        repository_path.unlink()
    elif repository_path.exists():
        raise FetchError(
            f"refusing to replace existing non-symlink benchmark file: {repository_path}"
        )
    repository_path.symlink_to(os.path.relpath(cache_path, repository_path.parent))


def fetch_assets(
    manifest: dict,
    *,
    cases: set[str],
    cache_dir: Path,
    repository_root: Path,
    verify_only: bool,
    create_symlinks: bool,
) -> list[Path]:
    selected = []
    for asset in manifest["assets"]:
        if asset.get("status") != "downloadable":
            continue
        if cases != {"all"} and asset.get("case") not in cases:
            continue
        selected.append(asset)
    if not selected:
        raise FetchError("no downloadable assets match the selected cases")

    ready = []
    for asset in selected:
        relative = safe_relative_path(str(asset["local_path"]))
        cache_path = cache_dir / relative
        expected = str(asset["sha256"])
        if not cache_path.is_file() or sha256(cache_path) != expected:
            if verify_only:
                raise FetchError(f"missing or invalid cached asset: {cache_path}")
            download(str(asset["url"]), cache_path)
        actual = sha256(cache_path)
        if actual != expected:
            cache_path.unlink(missing_ok=True)
            raise FetchError(
                f"checksum mismatch for {asset['url']}: expected {expected}, got {actual}"
            )
        if create_symlinks:
            materialize(cache_path, repository_root / relative)
        ready.append(cache_path)
    return ready


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", dest="cases", help="benchmark case to fetch; repeatable (default: all)")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--no-materialize", action="store_true", help="do not create ignored tutorial symlinks")
    parser.add_argument("--list", action="store_true", help="list manifest entries without downloading")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        manifest = load_manifest(args.manifest)
        if args.list:
            for asset in manifest["assets"]:
                print(f"{asset['case']}: {asset['status']}: {asset['upstream_path']} ({asset['license']})")
            return 0
        ready = fetch_assets(
            manifest,
            cases=set(args.cases or ["all"]),
            cache_dir=args.cache_dir.resolve(),
            repository_root=ROOT,
            verify_only=args.verify_only,
            create_symlinks=not args.no_materialize,
        )
    except (OSError, json.JSONDecodeError, KeyError, FetchError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(f"Verified {len(ready)} benchmark assets in {args.cache_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
