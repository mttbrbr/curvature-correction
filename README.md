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

> **Current scope:** the equations are complete for an inertial, non-rotating
> reference frame (`OmegaRot=0`). Explicitly rotating reference frames/MRF need
> the additional frame-rotation terms from the papers and are not yet supported.

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
coefficients remain available in the corresponding coefficient dictionary.

## Smoke-test tutorial

`tutorials/incompressible/simpleFoam/curvedChannel` is a small 2-D bend case:

```bash
cd tutorials/incompressible/simpleFoam/curvedChannel
./Allrun
```

The case verifies meshing, runtime selection and numerical execution. It is not
a physical validation case.
