# OpenFOAM rotation/curvature-correction models

Custom, runtime-selectable incompressible RANS models for OpenFOAM.com
(OpenCFD/ESI), targeting **v2412**.

> **Licensing:** project-authored files are MIT; the OpenFOAM adapter under
> `src/curvatureCorrectionModels` is GPL-3.0-or-later. Benchmark inputs are not
> redistributed and are fetched on demand under their upstream terms. See the
> root [licensing notice](LICENSE) and [third-party register](THIRD_PARTY.md).

The library provides two published rotation/curvature corrections:

- `kOmegaSSTCC` implements the Smirnov--Menter SST-RC correction. The bounded
  factor `fr1` multiplies the production used by both the `k` and `omega`
  equations, before the standard SST production limiter.
- `SpalartAllmarasCC` implements the Spalart--Shur SA-RC correction. It changes
  only the SA production term to `Cb1*(fr1-ft2)*Stilda*nuTilda`; the destruction
  term continues to use the unmodified `Stilda`.

Both implementations include the complete material derivative
`DS/Dt = dS/dt + U.grad(S)`. The time contribution is obtained from the stored
velocity time levels and is zero when a steady-state time scheme is selected.

> **Current scope:** inertial calculations and OpenFOAM.com v2412 MRF zones are
> supported. MRF angular velocity is discovered automatically from the registered
> `IOMRFZoneList`; SRF, moving-mesh, overset and 6DoF frame motion remain outside
> the validated scope.

## References

- P. R. Spalart and M. Shur, “On the Sensitization of Turbulence Models to
  Rotation and Curvature,” *Aerospace Science and Technology* 1(5), 1997,
  <https://doi.org/10.1016/S1270-9638(97)90051-1>.
- M. L. Shur, M. K. Strelets, A. K. Travin and P. R. Spalart, “Turbulence
  Modeling in Rotating and Curved Channels,” *AIAA Journal* 38(5), 2000,
  <https://doi.org/10.2514/2.1058>.
- P. E. Smirnov and F. R. Menter, “Sensitization of the SST Turbulence Model to
  Rotation and Curvature,” *Journal of Turbomachinery* 131, 2009,
  <https://doi.org/10.1115/1.3070573>.
- NASA Turbulence Modeling Resource definitions for
  [SST-RC](https://tmbwg.github.io/turbmodels/sst.html) and
  [SA-RC](https://tmbwg.github.io/turbmodels/spalart.html).

## Build

Load an OpenFOAM.com v2412 environment and run:

```bash
cd src/curvatureCorrectionModels
wmake libso
```

The library is written to `$FOAM_USER_LIBBIN/libcurvatureCorrectionModels.so`.

## Use

Add the library to `system/controlDict`:

```text
libs ("libcurvatureCorrectionModels.so");
```

Add a dedicated convection scheme for the material derivative to
`system/fvSchemes`. Keeping it separate from `div(phi,U)` avoids applying a
vector momentum scheme to the symmetric strain tensor:

```text
divSchemes
{
    // ...existing entries...
    div(phi,curvatureStrain) Gauss linear;
}
```

SST-RC configuration:

```text
RASModel kOmegaSSTCC;

kOmegaSSTCCCoeffs
{
    curvatureCorrection on;
    Cr1 1.0;
    Cr2 2.0;
    Cr3 1.0;
    frMin 0.0;
    frMax 1.25;
    writeCurvatureDiagnostics on; // optional, default off
}
```

SA-RC configuration:

```text
RASModel SpalartAllmarasCC;

SpalartAllmarasCCCoeffs
{
    curvatureCorrection on;
    Cr1 1.0;
    Cr2 12.0;
    Cr3 1.0;
    writeCurvatureDiagnostics on; // optional, default off
}
```

`curvatureCorrection` defaults to `off` in both custom models and must be
enabled explicitly. Non-finite coefficients and correction invariants fail fast.

SA-RC deliberately does not bound `fr1`: negative production in strongly
stabilising rotation/curvature is part of the published model. The SA case must
provide the standard `nuTilda` field instead of `k` and `omega`. All parent-model
coefficients remain available in the corresponding coefficient dictionary. The
custom models always write the current correction factor as `fr1` and the
detected per-cell frame angular velocity as `OmegaRot`. With
`writeCurvatureDiagnostics on`, they also write `curvatureIndicator`,
`productionBase`, `productionCorrected`, `productionRatio`, and `fr1Clipped`.
The SA clipping mask is identically zero because its published factor is
unbounded. In MRF cases `OmegaRot` is populated inside each active MRF cell zone
and is zero elsewhere.

## Deterministic verification and validation report

The mesh-independent kinematic oracle reproduces the implemented tensor
algebra and checks uniform flow, shear, pure strain, rigid rotation, opposite
curvature senses, small/large rates, near-zero denominators, clipping,
NaN/Inf rejection, disabled-correction equivalence, and orthogonal-basis
invariance:

```bash
python3 -m unittest discover -s tests -v
```

With OpenFOAM v2412 loaded, compiled manufactured fields exercise the temporal
and convective parts of `DS/Dt`, steady-state suppression and both correction
families. The automatic MRF mapping and a two-rank restart are checked
separately:

```bash
tests/run_material_derivative_regression.sh
tests/run_mrf_smoke.sh
```

The [remaining-work backlog](docs/remaining-work.md) records issue-ready acceptance gates.
The detailed [implementation audit](docs/implementation-audit.md) maps the
equations to source lines, compares them with NASA TMR, records dimensional and
asymptotic checks, and identifies numerical/model-scope risks. The MRF implementation accounts for OpenFOAM v2412 storing absolute `U` while
using relative convective fluxes. Its tensor algebra and runtime mapping are
tested, but a physical rotating-channel validation is still required. The
OpenFOAM v2412 SST parent has SST-2003 characteristics.

After completing the selected tutorials, collect validation tables into the
Git-ignored data directory and generate a report from a schema-v1 input file:

```bash
python3 tools/collect_curvature_validation_data.py
python3 tools/generate_curvature_validation_report.py \
    path/to/validation-input.json \
    --output-dir reports/curvature-validation \
    --strict-acceptance
```

The collector uses checksum-verified CC0 computational references and
trapezoidal spatial weights. The report generator validates schema and finite
values, calculates aggregate errors, and records incomplete critical gates.
Without `--strict-acceptance`, only explicit critical failures return status 1;
strict mode also fails for `not-verified` and `not-applicable` critical gates.
Generated tables, reference-series plots, and reports remain ignored local/CI
artifacts rather than release content.

Generated deliverables are:

- `reports/curvature-validation/curvature-validation-report.md`;
- `reports/curvature-validation/curvature-validation-report.html`;
- `reports/curvature-validation/aggregated-metrics.json`.

A report may include the two-grid neutral flat-plate and NACA cases, but it must
remain explicit about missing concave, physical rotating, resolved separated-flow,
Reynolds-number, and industrial validation. Neutrality plus one convex case is
not general model validation.

## Neutral flat-plate benchmark

`tutorials/incompressible/simpleFoam/tmrFlatPlate` implements the openly
licensed TMR 2DZP zero-pressure-gradient flat plate for native, correction-off
and correction-on SST/SA comparisons:

```bash
python3 tools/fetch_benchmark_data.py --case flatplate
cd tutorials/incompressible/simpleFoam/tmrFlatPlate
CURVATURE_SKIP_BUILD=1 ./Allrun quick
CURVATURE_SKIP_BUILD=1 ./Allrun medium
```

The 69x49 and 137x97 runs recover correction-off fields exactly. Across both
grids, the maximum correction-on `Cf` change is 0.773% for SST-RC and 0.459%
for SA-RC, so both pass the declared neutral tolerance. The analysis also
compares `Cf(Re_theta)` with TMR CFL3D and Karman-Schoenherr references and
`u+(y+)` with Coles. Source, licence and selection rationale are recorded in
[benchmark selection](docs/benchmark-selection.md).


## TMR NACA 0012 SA-RC benchmark

`tutorials/incompressible/simpleFoam/tmrNaca0012` compares OpenFOAM SA and the
local SA-RC model with CC0 CFL3D data on nested 113x33 and 225x65 grids at
0, 10 and 15 degrees. The quantitative 0/10-degree gates pass on both grids.
On 225x65 at 10 degrees, SA-RC reaches 2.53% CL error and 1.51% CD error and
moves both coefficients toward the independent CFL3D SA-RC result. The
15-degree branch is explicitly diagnostic and remains under-resolved.

```bash
python3 tools/fetch_benchmark_data.py --case naca0012
cd tutorials/incompressible/simpleFoam/tmrNaca0012
CURVATURE_SKIP_BUILD=1 ./Allrun quick
CURVATURE_SKIP_BUILD=1 ./Allrun medium
```

See the [NACA tutorial documentation](tutorials/incompressible/simpleFoam/tmrNaca0012/README.md)
and [third-party register](THIRD_PARTY.md) for acceptance scope, provenance,
CC0 terms and the excluded experimental data.

## Validation tutorial

`tutorials/incompressible/simpleFoam/smitsCurve` implements the TMBWG 2DCC
Smits convex-curvature case with external, checksum-verified inputs. Its `Allrun` script builds the
library, imports a 6,144-cell PLOT3D mesh, runs native/on/off variants of SST and
SA, checks `Cf` and `fr1`, exports VTK, and creates an SVG comparison against
TMBWG CFL3D computational data:

```bash
python3 tools/fetch_benchmark_data.py --case smits
cd tutorials/incompressible/simpleFoam/smitsCurve
./Allrun
```

On the development machine the full six-case workflow completes in about 80
seconds. The baseline run gave exact agreement between each native model and
its custom counterpart with the correction off. With the correction on,
convex-wall `Cf` near the bend decreased by 18.3% for SST-RC and 19.8% for
SA-RC. See the [tutorial documentation](tutorials/incompressible/simpleFoam/smitsCurve/README.md)
for acceptance criteria, limitations, and ParaView instructions.

For stronger coverage, `Allrun-restartMPI` validates serial restart and a
four-rank MPI restart, while `Allrun-gridStudy` runs native and corrected SST/SA
on nested 129x49, 257x97, and 513x193 grids. The quick grid is fetched from
the pinned CC0 source; medium/fine archive files are licence-held and must be
supplied explicitly via `SMITS_MEDIUM_GRID` and `SMITS_FINE_GRID`. The study
requires converged solutions and contracting `Cf` differences.

## Smoke-test tutorial

`tutorials/incompressible/simpleFoam/curvedChannel` is a small 2-D bend case:

```bash
cd tutorials/incompressible/simpleFoam/curvedChannel
./Allrun
```

The case verifies meshing, runtime selection and numerical execution. It is not
a physical validation case.

## Visible tip-vortex demonstration

`tutorials/incompressible/simpleFoam/flapTipVortex` is a finite-span 3-D flap
case designed to make the rolled-up tip vortices easy to inspect. It runs
native SST and SST-RC on the same locally refined mesh and writes `Q`,
`vorticity`, forces, and VTK output for direct ParaView comparison:

```bash
cd tutorials/incompressible/simpleFoam/flapTipVortex
./Allrun
```

This is a qualitative demonstration, not a validation case. See the
[case documentation](tutorials/incompressible/simpleFoam/flapTipVortex/README.md)
for visualization settings and limitations.
