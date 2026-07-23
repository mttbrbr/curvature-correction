# OpenFOAM rotation/curvature-correction models

Custom, runtime-selectable incompressible RANS models for OpenFOAM.com
(OpenCFD/ESI), targeting **v2412**.

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
}
```

SA-RC deliberately does not bound `fr1`: negative production in strongly
stabilising rotation/curvature is part of the published model. The SA case must
provide the standard `nuTilda` field instead of `k` and `omega`. All parent-model
coefficients remain available in the corresponding coefficient dictionary. The
custom models write the current correction factor as `fr1` and the detected
per-cell frame angular velocity as `OmegaRot` at every normal OpenFOAM write
time. In MRF cases `OmegaRot` is populated inside each active MRF cell zone and
is zero elsewhere; no additional model coefficient is required.

## Deterministic verification and validation report

The mesh-independent kinematic oracle reproduces the implemented tensor
algebra and checks uniform flow, shear, pure strain, rigid rotation, opposite
curvature senses, small/large rates, near-zero denominators, clipping,
NaN/Inf rejection, disabled-correction equivalence, and orthogonal-basis
invariance:

```bash
python3 -m unittest discover -s tests -v
```

With OpenFOAM v2412 loaded, the automatic MRF mapping and a two-rank restart
are exercised with:

```bash
tests/run_mrf_smoke.sh
```

The detailed [implementation audit](docs/implementation-audit.md) maps the
equations to source lines, compares them with NASA TMR, records dimensional and
asymptotic checks, and identifies numerical/model-scope risks. The MRF implementation accounts for OpenFOAM v2412 storing absolute `U` while
using relative convective fluxes. Its tensor algebra and runtime mapping are
tested, but a physical rotating-channel validation is still required. The
OpenFOAM v2412 SST parent has SST-2003 characteristics.

After completing the three-grid tutorial, refresh the lightweight committed
evidence and generate Markdown, HTML, JSON metrics, and SVG plots with:

```bash
python3 tools/collect_curvature_validation_data.py
python3 tools/generate_curvature_validation_report.py \
    data/curvature-validation/validation-input.json \
    --output-dir reports/curvature-validation
```

The collector limits profile comparison to the experimental overlap and uses
trapezoidal spatial weights. The report generator validates schema and finite
values, calculates RMSE, normalized RMSE, MAE, bias, maximum and relative
errors, classifies each dataset, lists non-converged/missing-reference cases,
and exits non-zero for invalid data or a declared critical failure.

Generated deliverables are:

- `reports/curvature-validation/curvature-validation-report.md`;
- `reports/curvature-validation/curvature-validation-report.html`;
- `reports/curvature-validation/aggregated-metrics.json`.

The report is deliberately explicit about missing neutral, concave, rotating,
separated-flow, Reynolds-number, and industrial validation. A successful
regression or one convex-curvature benchmark is not treated as general model
validation.

## Validation tutorial

`tutorials/incompressible/simpleFoam/smitsCurve` is a self-contained version of
the TMBWG 2DCC Smits convex-curvature case. Its `Allrun` script builds the
library, imports a 6,144-cell PLOT3D mesh, runs native/on/off variants of SST and
SA, checks `Cf` and `fr1`, exports VTK, and creates an SVG comparison against
TMBWG CFL3D and experimental data:

```bash
cd tutorials/incompressible/simpleFoam/smitsCurve
./Allrun
```

On the development machine the full six-case workflow completes in about 80
seconds. The baseline run gave exact agreement between each native model and
its custom counterpart with the correction off. With the correction on,
convex-wall `Cf` near the bend decreased by 18.5% for SST-RC and 18.8% for
SA-RC. See the [tutorial documentation](tutorials/incompressible/simpleFoam/smitsCurve/README.md)
for acceptance criteria, limitations, and ParaView instructions.

For stronger coverage, `Allrun-restartMPI` validates serial restart and a
four-rank MPI restart, while `Allrun-gridStudy` runs native and corrected SST/SA
on the nested 129x49, 257x97, and 513x193 TMBWG grids. The latter requires
converged solutions and contracting `Cf` differences between successive levels.
`Allrun-extended` executes the complete suite with a 30-minute total budget.

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
