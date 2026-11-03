# Smits 2-D convex-curvature validation

This tutorial exercises both curvature-correction models on the TMBWG 2DCC
case: a turbulent boundary layer passes through a rapid 30-degree bend. Convex
curvature stabilises the boundary layer on the lower wall, so a working
correction should reduce its skin-friction coefficient relative to the native
model near the bend.

Three automated workflows cover quick validation, restart/MPI consistency, and
grid sensitivity. All use OpenFOAM.com v2412 and keep their generated data under
`runs/`.

Fetch the pinned CC0 quick-grid/reference inputs once from the repository root:

```bash
python3 tools/fetch_benchmark_data.py --case smits
```

`Allrun` and `Allrun-restartMPI` fail early with this exact instruction if the
cache is absent or a checksum differs.

## Quick six-case validation

Load the OpenFOAM environment and run:

```bash
cd tutorials/incompressible/simpleFoam/smitsCurve
./Allrun
```

This builds the library, imports the 129x49 PLOT3D grid, runs six serial cases,
checks the results, exports VTK, and generates a comparison plot. It takes about
90 seconds on the development machine and enforces a 30-minute limit.

| Directory | Turbulence model | Correction |
| --- | --- | --- |
| `sst` | native `kOmegaSST` | n/a |
| `sstcc-off` | `kOmegaSST2003RC` | off |
| `sstcc-on` | `kOmegaSST2003RC` | on |
| `sa` | native `SpalartAllmaras` | n/a |
| `sacc-off` | `SpalartAllmarasNoFt2RC` | off |
| `sacc-on` | `SpalartAllmarasNoFt2RC` | on |

Acceptance requires the custom models with the correction off to match their
native parents, corrected SST and SA to reduce convex-wall `Cf` near the bend,
and `fr1` to be valid. SST `fr1` must remain inside [0, 1.25]; the published SA
correction is deliberately unbounded. In addition to `Cf`, the runner compares
native and correction-off internal fields at a relative L-infinity tolerance of
1e-12: `U`, `p`, `k`, `omega`, and `nut` for SST; `U`, `p`, `nuTilda`, and
`nut` for SA.

Results are `runs/summary.csv`, `runs/sst-off-fields.csv`,
`runs/sa-off-fields.csv`, and `runs/cf-comparison.svg`.

## Restart and MPI consistency

Run:

```bash
./Allrun-restartMPI
```

The test creates identical SST-RC initial conditions and compares:

1. an uninterrupted serial solution;
2. a serial solution restarted from iteration 50;
3. the same iteration-50 checkpoint restarted on four MPI ranks and
   reconstructed with `reconstructPar`.

All three paths continue to the configured residual convergence. Internal
`U`, `p`, `k`, `omega`, `nut`, and `fr1` fields must agree within a relative
L-infinity tolerance of 5e-4. The test takes about 40 seconds and writes
`runs/restart-mpi/summary.csv`. Change the rank count with
`CURVATURE_MPI_PROCS`; it must be at least two.

## Three-grid convergence study

Run:

The medium/fine NASA archive files are on redistribution hold. Supply lawfully
obtained PLOT3D files explicitly (plain or gzip-compressed):

```bash
SMITS_MEDIUM_GRID=/path/to/257x97.p2dfmt.gz \
SMITS_FINE_GRID=/path/to/513x193.p2dfmt.gz \
./Allrun-gridStudy
```

This executes native and corrected SST and SA on three official nested TMBWG
grids, using four MPI ranks by default:

| Level | Points | Cells |
| --- | ---: | ---: |
| level 3 | 129x49 | 6,144 |
| level 2 | 257x97 | 24,576 |
| level 1 | 513x193 | 98,304 |

Each finer grid is initialised with `mapFields` from the converged preceding
level. Derived fields `fr1` and `wallShearStress` are discarded after mapping
and recomputed, preventing diagnostics from contaminating a restart. The
pressure residual threshold is 1e-5 to avoid parallel GAMG jitter; velocity and
turbulence fields remain at 1e-6.

The script requires every solver to print `SIMPLE solution converged`. It then
checks that the curvature response remains physical, that the 257-to-513 `Cf`
relative RMS difference is below 5%, and that successive grid differences
contract. The summary also records an observed order and approximate fine-grid
GCI based on those profile norms, solver iterations, execution time,
corrected/base cost ratios, and `fr1` extrema. Solver logs are rejected if they
contain NaN or Inf. A clean run takes about 18 minutes and has a 25-minute
limit.

Outputs include:

- `runs/grid-study/summary.csv`;
- `runs/grid-study/sst-grid-convergence.svg`;
- `runs/grid-study/sa-grid-convergence.svg`;
- per-case `cf.csv` and `.foam` markers.

A stopped study can safely reuse only completed cases with:

```bash
CURVATURE_GRID_RESUME=1 ./Allrun-gridStudy
```

## Inspect in ParaView

Open any marker such as `runs/sstcc-on/sstcc-on.foam` or its
`runs/<case>/VTK/<case>_*.vtm` export. For a useful side-by-side comparison:

1. open `sst/sst.foam` and `sstcc-on/sstcc-on.foam`;
2. select the latest time and colour by `fr1`, `nut`, `U`, or
   `wallShearStress`;
3. repeat with `sa/sa.foam` and `sacc-on/sacc-on.foam`;
4. use **Plot Over Line** along the convex lower wall for a spatial profile.

Fine-grid markers are under paths such as
`runs/grid-study/level1-513x193/sstcc-on/sstcc-on.foam`. `fr1` is present only
in custom-model cases and is written automatically.

## Mesh-quality note

All grids are wall resolved and intentionally contain very high-aspect-ratio
boundary-layer cells. `checkMesh` reports exactly one failure at each level for
this reason: maximum aspect ratios are approximately 6,122, 6,762, and 7,100.
Topology, maximum non-orthogonality (below 28.4 degrees), and skewness are
acceptable. The runners verify cell count, connectivity, skewness, and the
known single aspect-ratio failure.

## Interpretation and limitations

The quick workflow is a regression and trend-validation test. The optional
three-grid workflow adds a meaningful grid-sensitivity check but is not a
formal uncertainty quantification or an asymptotic-order study. TMBWG reference
computations are compressible at the essentially incompressible Mach number
0.093, whereas these OpenFOAM cases are incompressible, so small absolute
differences are expected.

The TMBWG pages describe the experiment, boundary conditions, nested grids, and
reference solutions:

- <https://tmbwg.github.io/turbmodels/smitscurve_val.html>
- <https://tmbwg.github.io/turbmodels/smitscurve_grids.html>
- <https://github.com/TMBWG/turbmodels>

No benchmark bytes are versioned under `resources/`. Immutable URLs, upstream
paths, licences, and SHA-256 digests are recorded in
`data/benchmark-sources.json`; `resources/README.md` documents materialisation.
The held experimental series is not plotted. Run `./Allclean` to remove results.
