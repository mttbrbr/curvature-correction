# Finite-span flap tip-vortex demonstration

This is a deliberately visual 3-D demonstration, not a validation case. A
rectangular flat flap of chord `c = 1 m`, span `b = 4 m` (`AR = 4`) and 4%
thickness is set at 12 degrees to a `20 m/s` free stream. Both tips are free,
so the pressure difference between the two faces rolls up into two
counter-rotating streamwise vortices.

The same locally refined mesh is used for native `kOmegaSST` and corrected
`kOmegaSSTCC`, making diffusion, peak vorticity and core size comparable.

## Run

Load OpenFOAM.com v2412, then execute:

```bash
cd tutorials/incompressible/simpleFoam/flapTipVortex
./Allrun
```

The script builds the library, creates the snappyHexMesh mesh once, runs both
variants on four MPI ranks and exports their last times to VTK. Expect roughly
0.5--1.5 million cells after snapping, depending on the exact mesher version.

## See the vortex in ParaView

1. Open the last dataset under `runs/SSTCC/VTK/`.
2. Apply **Contour** to `Q`; start with `Q = 1000 1/s2` and adjust between 300
   and 3000 until the two cores are continuous.
3. Colour the contour by the magnitude of `vorticity`.
4. Add a **Slice** normal to `x` at `x = 2`, `4`, or `6` and colour it by
   `vorticity_X`. The two tip cores have opposite signs.
5. Load `runs/SST/VTK/` and use identical contour and colour limits.

The case also exports the `Q = 1000` isosurface and the three wake sections
directly under each run's `postProcessing/vortexSurfaces` directory.
Additionally, each run contains an empty ParaView marker (`SST.foam` or
`SSTCC.foam`) and a standalone final isosurface at
`postProcessing/Qsurface/<latest-time>/Q.vtp`.

`fr1` is written by the corrected model in the OpenFOAM time directory. Open
the case with `paraFoam -case runs/SSTCC` to inspect it. Force histories are in
`postProcessing/forces`.

A steady RANS case shows the mean rolled-up tip vortex, not instantaneous
shedding. The flat plate and wall-function mesh favour a strong, visible
vortex rather than high-fidelity loads. Check mesh and iterative convergence
before treating quantitative differences as validation evidence.

In a 100-iteration smoke run on the development machine, the mean of the two
tip-core peaks at `x/c = 6` was `Q = 1047 1/s2`, `|omega_x| = 65.2 1/s` for
native SST and `Q = 1263 1/s2`, `|omega_x| = 71.5 1/s` for SST-RC. Thus the
corrected run retained about 21% more peak Q and 9.8% more axial vorticity six
chords downstream. These early-iteration values demonstrate visibility and
sensitivity only; use the completed runs for quantitative comparisons.
