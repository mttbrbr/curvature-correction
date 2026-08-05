# TMR NACA 0012 SA-RC benchmark

This tutorial compares native OpenFOAM SA and the local SA-RC implementation
with CC0 CFL3D results for the TMR NACA 0012 validation case.

The case uses incompressible `simpleFoam`, chord Reynolds number 6 million,
freestream `nuTilda/nu=3`, a farfield approximately 500 chords away, and angles
of attack 0, 10, and 15 degrees. The 113x33 and 225x65 grids are exactly nested.

```bash
# From the repository root, after loading OpenFOAM.com v2412:
python3 tools/fetch_benchmark_data.py --case naca0012
cd tutorials/incompressible/simpleFoam/tmrNaca0012
./Allrun quick
./Allrun medium
```

Set `NACA_ANGLES="0 10"` to run a subset during development. The workflow
imports the PLOT3D mesh, identifies the airfoil/farfield/outlet patches, runs
SA and SA-RC, and compares averaged final `CL` and `CD` with CFL3D.

The CFL3D reference used a 897x257 grid, compressible Mach 0.15 conditions, and
a farfield point-vortex correction. The local cases are incompressible and do
not apply that correction; the 500-chord farfield makes its effect small but
not identically zero. The supplied reference is representative, not truth, and
the two local grids do not establish formal grid independence.

## Acceptance and observed results

Angles 0 and 10 degrees are quantitative gates; 15 degrees is a diagnostic
near-stall point. The quick limits are deliberately broad because 113x33 is a
smoke grid. On 225x65, the limits are 25% for CD at 0 degrees and 5%/10% for
CL/CD at 10 degrees. Both grids also require finite sane coefficients, `abs(CL) <= 0.02` at 0
degrees, final-five-sample CD span below 8%, CL relative span below 1% at
10 degrees (absolute span below 0.01 at 0 degrees), and the published
10-degree SA-RC trend:
increased CL and reduced CD relative to SA.

The retained medium results pass those gates. At 10 degrees, SA-RC gives
`CL=1.06707` and `CD=0.0119715`, versus CFL3D SA-RC `1.09474` and `0.0117939`
(2.53% and 1.51% error). Refinement contracts the corresponding quick-grid
errors from 8.39% and 34.45%. At 15 degrees, angle continuation avoids a
numerical divergence but the 225x65 solution remains on an under-resolved
separated branch (`CL=0.6753`, `CD=0.1215`); it is recorded, not accepted as
validation.

A 15-degree rerun can reuse the converged 10-degree fields and mesh:

```bash
NACA_ANGLES="15" NACA_REUSE_LEVEL=1 CURVATURE_SKIP_BUILD=1 ./Allrun medium
```

Only pinned TMBWG computational files are downloaded into the ignored cache.
Experimental data and publisher figures are excluded. Provenance, audited
revision, licence, upstream paths, and checksums are in `resources/README.md`,
`data/benchmark-sources.json`, and `THIRD_PARTY.md`.
