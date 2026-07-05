# OpenFOAM curvature-correction models

Custom, runtime-selectable RANS models for OpenFOAM.com (OpenCFD/ESI),
initially targeting **v2412**.

The first model, `kOmegaSSTCC`, derives from the standard `kOmegaSST` and
multiplies its already-limited turbulent kinetic-energy production term by a
bounded curvature/rotation factor `fr`.

> **Physics status:** the current algebraic `fr` expression is an intentionally
> isolated placeholder. It uses strain and rotation invariants, but it is not a
> validated implementation of a specific published correction. Replace
> `common/curvatureFunctions.C` after selecting and validating a formulation.

## Build

Load an OpenFOAM.com v2412 (or compatible later) environment and run:

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

Select the model in `constant/turbulenceProperties`:

```text
simulationType RAS;

RAS
{
    RASModel        kOmegaSSTCC;
    turbulence      on;
    printCoeffs     on;

    kOmegaSSTCCCoeffs
    {
        curvatureCorrection on;
        Cr1                 1.0;
        Cr2                 2.0;
        Cr3                 1.0;
        frMin               0.0;
        frMax               1.25;
    }
}
```

All standard `kOmegaSST` coefficients remain available in the same coefficient
dictionary.

## Smoke-test tutorial

`tutorials/incompressible/simpleFoam/curvedChannel` is a small 2-D bend case.

```bash
cd tutorials/incompressible/simpleFoam/curvedChannel
./Allrun
```

The case is intended to verify meshing, runtime selection and basic numerical
execution. It is not a physical validation case.
