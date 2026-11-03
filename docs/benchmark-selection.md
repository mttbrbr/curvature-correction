# Benchmark selection for rotation/curvature validation

This document records candidate benchmarks, source provenance, reuse terms, and
the validation question answered by each case. It is an engineering provenance
record, not legal advice.

## Reuse policy

External numerical data and grids are not vendored. Authorised inputs are
downloaded to an ignored cache from immutable URLs and materialised as ignored
symlinks. Papers are cited but not redistributed; publisher figures and
substantial text are not copied. Every downloadable asset records its source,
revision, upstream licence, SHA-256 checksum, and local path in
`data/benchmark-sources.json`. Ambiguous files remain on hold and require an
explicit user-supplied path.

The current NASA/AIAA Turbulence Modeling Resource (TMR) repository is
distributed under CC0-1.0:

- <https://github.com/TMBWG/turbmodels>
- <https://github.com/TMBWG/turbmodels/blob/main/LICENSE>

NASA also notes that a restrictive notice attached to third-party material must
still be respected. For that reason the licence attached to the concrete
repository or dataset, rather than NASA hosting alone, is the acceptance gate.

## Selected sequence

| Priority | Case | Physics and validation question | Reference comparison | Reuse status | Disposition |
| --- | --- | --- | --- | --- | --- |
| P0 | TMR 2DZP flat plate | Neutral attached boundary layer: does enabling RC leave a nominally straight, non-rotating flow unchanged? | Karman-Schoenherr `Cf`, Coles `u+`, and independent CFL3D/FUN3D SA/SST results | TMR repository CC0-1.0 | Implemented on 69x49 and 137x97 grids; neutral criterion passes |
| P1 | TMR NACA 0012 | Curved streamlines and near-stall behaviour; direct independent check of SA-RC | CFL3D, FUN3D and NTS SA-RC plus experimental lift, drag and pressure data | TMR repository CC0-1.0 | Implemented on 113x33 and 225x65; 0/10-degree gates pass, 15-degree branch remains diagnostic |
| P1 | TMR axisymmetric rotating pipe (Zaets et al.) | Physical rotation, swirl and partial relaminarisation; closes the largest MRF/rotation evidence gap | Experimental axial/circumferential velocity, shear stress and turbulent kinetic energy | TMR repository CC0-1.0; page warns that the case and BCs remain under construction | Prototype after 2DZP; do not make a pass/fail merge gate until BC sensitivity is resolved |
| P1/P2 | Appelbaum-Kloker-Wenzel curved/flat DNS pairs | Isolates curvature from pressure-gradient effects using paired curved and flat turbulent boundary layers | Four DNS datasets: ZPG-C, ZPG-F, PG-C and PG-F | CC BY 4.0, DOI `10.18419/DARUS-4345` | Strongest independent curvature dataset, but about 1.8 GB; download on demand and commit only reduced profiles |
| P2 | TMR NASA wall-mounted hump | Curved separated flow and reattachment robustness | Experiment and several baseline RANS results | TMR repository CC0-1.0 | Useful after neutral and direct RC cases; not a correction-specific reference |
| P2 | TMR NACA 4412 trailing-edge separation | Curved separated airfoil, with published SA-RC comparisons | Two independent SA-RC codes and experiment | TMR repository CC0-1.0 | Qualitative only: TMR explicitly identifies wind-tunnel and scaling inconsistencies |

## Completed TMR 2DZP result

The `tutorials/incompressible/simpleFoam/tmrFlatPlate` case uses a
pinned manifest for the two smallest nested TMR grids and computational reference curves. Inputs
are fetched on demand under CC0-1.0 and verified with SHA-256. Six native/on/off cases were run on both grids. Correction-off
fields are identical to their parents. The largest correction-on `Cf` deviation
across the two grids is 0.773% for SST-RC and 0.459% for SA-RC, below the 1%
case limit and the report-wide 2% neutral tolerance.

On the 137x97 grid, SST CFL3D `Cf` nRMSE changes from 0.01238 to 0.01274 and SA
from 0.009053 to 0.008631. The Coles `u+` comparison changes only slightly.
This passes the neutral-baseline question; it does not demonstrate a general
accuracy improvement or formal three-grid independence.


## Completed TMR NACA 0012 result

The `tutorials/incompressible/simpleFoam/tmrNaca0012` case uses a pinned download manifest for two nested
TMR grids and CC0 CFL3D SA/SA-RC force, pressure and skin-friction curves. On 225x65 at 10 degrees, OpenFOAM SA-RC
moves CL from 1.06195 to 1.06707 toward 1.09474 and CD from 0.012981 to
0.011972 toward 0.011794. The SA-RC errors contract from 8.39%/34.45% on the
113x33 grid to 2.53%/1.51% on 225x65 for CL/CD.

The 15-degree calculation is not accepted as validation: even with continuation
from 10 degrees, the medium grid converges to a strongly separated,
under-resolved branch. A 449x129-or-finer grid and a numerical-scheme study are
required before using that point as a gate. Experimental curves were not
vendored because their concrete redistribution provenance is ambiguous.

## Cases not vendored

The University of Tokyo/Nagoya rotating-channel DNS pages contain relevant
spanwise-rotation profiles, but the public page carries an “all rights
reserved” notice and no dataset-specific open licence was found. Links may be
recorded for literature review, but their files must not be committed without
written permission or a clearer licence.

Likewise, plots digitised from publisher PDFs are excluded. A numeric fact may
be used in analysis with a citation, but digitisation is not treated as a
substitute for an openly licensed dataset.

## Acceptance metrics

### TMR 2DZP

- `Cf(Re_theta)` against the Karman-Schoenherr curve over the available
  overlap, reporting RMSE, bias and maximum relative error.
- `u+(y+)` against the Coles profile near `Re_theta=10000`.
- Native versus custom-correction-off field equality.
- Native versus correction-on `Cf` and velocity-profile differences. These
  quantify neutrality; they are not expected to be bitwise zero.
- `fr1` range and production-weighted departure from one.
- At least two nested grids before using the case as a merge gate.

### TMR NACA 0012

- SA-RC `CL` and `CD` at 0, 10 and 15 degrees against CFL3D, FUN3D and NTS.
- Surface `Cp` and `Cf` against the TMR curves.
- Baseline SA, SA-RC and experiment must remain separate comparisons.
- The initial acceptance band is the published inter-code spread, augmented by
  the local two-grid discretisation difference.

### TMR rotating pipe

- Axial and circumferential mean velocity profiles at the non-rotating and
  rotating measurement stations.
- Turbulent shear stress and turbulent kinetic energy profiles.
- Trend with rotation ratio `N`, beginning with `N=0` and `N=0.6`.
- Inertial rotating-wall and MRF representations should agree within their
  numerical uncertainty before claiming frame-independent validation.

### Curved/flat DNS pairs

- Compare the *difference* between paired curved and flat cases, not only each
  absolute profile. This is the principal defence against confusing pressure
  gradient and curvature effects.
- Use wall shear, mean velocity, Reynolds shear stress and boundary-layer
  integral quantities.
- Keep raw and derived profiles in ignored local/CI artifacts with attribution,
  transformation records, and checksums; do not redistribute them by default.
