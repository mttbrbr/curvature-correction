# TMR axisymmetric rotating pipe

This tutorial is a staged physical check of the OpenFOAM.com v2412 MRF adapter
using the Zaets et al. rotating-pipe profiles hosted by the Turbulence Modeling
Resource. It does not close the physical-validation gate: the TMR case and
boundary conditions remain under construction and reference profiles exist
only for `N=0` and `N=0.6`.

The incompressible case uses `D=0.06 m`, `Uref=10 m/s`, `Re=37000`, a
non-rotating development length of `150D`, and a measurement station `25D`
into the rotating section. The `N=0.3` cases are sensitivity calculations,
not experimental validation. The fixed bulk flow is `0.8159892677 Uref`,
obtained by area-integrating the published `N=0` profile and applying an
anchored one-seventh-power extrapolation over its unmeasured wall tail.

Run after loading OpenFOAM.com v2412:

```bash
python3 tools/fetch_benchmark_data.py --case rotatingpipe
cd tutorials/incompressible/simpleFoam/tmrRotatingPipe
./Allrun smoke
./Allrun quick
./Allrun medium
```

`smoke` runs fixed 1000-iteration native and SST-RC inertial/MRF pairs on a
dedicated 98k-cell full-cylinder O-grid. It checks completion, finite values,
continuity and emits qualified frame metrics; it is not a RANS-convergence
claim.

`quick` uses about 238k cells and runs the full two-family matrix including
`N=0.3`, correction-off checks and both frame representations. `medium` uses
about 855k cells and repeats the experimental `N=0/0.6` comparisons. The
five-block cross-section avoids a collapsed axis and is checked with
`checkMesh -allTopology -allGeometry`. Generated cases, profiles, metrics and
plots are stored below `runs/` and remain ignored.

Both quick and medium include native MRF controls for SST and SA at `N=0.6`.
Every run also refreshes `runs/frame-resolution.{csv,json}`. The summary
records frame-error contraction across the available mesh levels and remains
`not-verified` until smoke, quick and medium metrics all exist; it is evidence
for a later grid-independence decision, not a grid gate by itself.

Deterministic gates are solver/finite checks, correction-off equivalence and
qualified frame equivalence (relative RMS at most 1% and relative maximum at
most 2%). A corrected-frame mismatch is a failure only when the matching
native SST/SA frame control passes; otherwise it is reported as
`not-verified`, with both metric tables retained.

Both representations prescribe `rotatingWallVelocity`. OpenFOAM's MRF
formulation stores absolute `U`, and `MRF.correctBoundaryVelocity` imposes the
same `Omega x r` value on a wall included in the zone; retaining the explicit
boundary condition also preserves that absolute value during turbulence
correction and standalone sampling.
Physical accuracy, formal grid independence, scheme sensitivity and a second
non-zero reference remain explicitly `not-verified`.
