# TMR 2DZP neutral flat-plate benchmark

This tutorial imports the NASA/AIAA TMBWG zero-pressure-gradient flat-plate
grids and runs native, correction-off, and correction-on variants of SST and
Spalart-Allmaras. Its primary purpose in this repository is to test whether the
rotation/curvature correction remains neutral in a straight, non-rotating,
attached boundary layer.

The incompressible setup uses `U_inf=1`, `nu=2e-7`, and therefore a Reynolds
number per unit length of five million. SA uses `nuTilda/nu=3`. SST uses
`k=2.25e-7` and `omega=125`, corresponding to the TMR freestream eddy-viscosity
ratio of 0.009. The upper boundary is an OpenFOAM freestream boundary so that
the displacement thickness does not artificially accelerate the outer flow.

From the repository root, fetch the pinned CC0 inputs, then run the 69x49 grid:

```bash
python3 tools/fetch_benchmark_data.py --case flatplate
cd tutorials/incompressible/simpleFoam/tmrFlatPlate
./Allrun
```

Run the 137x97 grid:

```bash
TMR_FLATPLATE_GRID=medium ./Allrun
```

The analyser reports:

- `Cf(Re_theta)` errors against Karman-Schoenherr theory;
- errors against TMR CFL3D SA and SST reference curves;
- `u+(y+)` errors against the TMR Coles profile;
- correction-on neutrality and correction-off equivalence;
- `fr1` ranges and mean absolute departure from one.

The correction-on acceptance limit is a maximum pointwise change of 1% in
plate `Cf`. This is a regression/neutrality threshold, not a claim that the
correction factor itself must be uniformly one.

Source: <https://tmbwg.github.io/turbmodels/flatplate_val.html>

The grids and computational curves are downloaded from a pinned revision of
<https://github.com/TMBWG/turbmodels> under CC0-1.0, checksum verified, and kept
out of Git. See `resources/README.md`, `data/benchmark-sources.json`, and
`THIRD_PARTY.md`.
