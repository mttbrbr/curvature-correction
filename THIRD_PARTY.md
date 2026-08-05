# Third-party software and benchmark provenance

Status date: 2026-08-05. This is an engineering compliance record, not legal advice.

## Release decision

The repository is intentionally multi-licensed. Original project files use MIT;
the OpenFOAM adapter under `src/curvatureCorrectionModels/` uses
GPL-3.0-or-later. External benchmark datasets are not committed: reproducible,
checksum-verified downloads are used only where redistribution rights are clear.

| Component | Treatment |
| --- | --- |
| Project-authored tools, tests, docs, and tutorial automation | MIT |
| OpenFOAM adapter | GPL-3.0-or-later |
| TMBWG computational inputs at the pinned revision | external download, CC0-1.0 |
| Smits experimental series and medium/fine grid archive | hold; user-supplied only |
| Appelbaum--Kloker--Wenzel DNS data | candidate external import, CC BY 4.0 |

## OpenFOAM adapter

The adapter targets OpenFOAM.com v2412 and links against GPL OpenFOAM libraries.
`SpalartAllmarasCC.C` also modifies the transport-equation structure from the
OpenFOAM Spalart--Allmaras implementation. Those files therefore carry
`SPDX-License-Identifier: GPL-3.0-or-later`, applicable upstream copyright and
modification notices. The full licence is `LICENSES/GPL-3.0-or-later.txt`.
OpenFOAM itself is a system/development dependency and is not bundled.

Official licensing information: <https://www.openfoam.com/licensing>.

## TMBWG Turbulence Modeling Resource

- repository: <https://github.com/TMBWG/turbmodels>
- pinned revision: `ddf1840b6e928d1de60d920cafe431ed476f36a1`
- upstream licence: CC0-1.0
- licence at the pinned revision: <https://github.com/TMBWG/turbmodels/blob/ddf1840b6e928d1de60d920cafe431ed476f36a1/LICENSE>

The downloadable scope is limited to TMBWG computational grids, neutral-map
files, and CFL3D reference curves for the flat-plate, NACA 0012, and quick Smits
cases. Exact URLs and SHA-256 digests are in `data/benchmark-sources.json`. Run
`python3 tools/fetch_benchmark_data.py --case CASE` to place them in the ignored
`.cache/` tree and create ignored tutorial symlinks. No downloaded bytes are
part of a Git release.

CC0 only covers rights held by the affirmer. For that reason the following are
not automatically downloaded or redistributed:

- `Smits_curvature_expdata_reworked2012.dat`, whose underlying experiment has
  third-party authorship;
- the medium/fine files from the NASA `smitscurve-grids.zip` archive, for which
  no explicit file-level redistribution licence was verified.

The three-grid study accepts medium and fine grids only through explicit
user-supplied paths. Public report generation uses the CC0 CFL3D curves, not the
held experimental series.

## Candidate DNS dataset

Appelbaum, Kloker, and Wenzel, “Data for: A systematic DNS approach to isolate
wall-curvature effects in spatially developing boundary layers,” DaRUS V3.1,
<https://doi.org/10.18419/DARUS-4345>, is licensed CC BY 4.0. It is not currently
used or downloaded. A future importer must preserve attribution, licence link,
version, transformation record, and change notice; raw data must remain
external.

## Publications

The implementation uses model equations described by Shur et al.
(<https://doi.org/10.2514/2.1058>) and Smirnov and Menter
(<https://doi.org/10.1115/1.3070573>). Papers are cited, not copied.

## Release checklist

- [x] Mixed MIT/GPL licensing is explicit at repository and component level.
- [x] OpenFOAM-derived adapter files carry GPL SPDX and provenance notices.
- [x] External computational inputs are pinned and checksum verified.
- [x] External inputs, derived reference-series tables, and plots are excluded
  from Git releases.
- [x] Ambiguous Smits materials are held and never auto-downloaded.
- [x] CC BY DNS data remains external until an attribution-aware importer exists.
