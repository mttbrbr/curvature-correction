# Rotation/curvature-correction models for OPENFOAM®

Custom, runtime-selectable incompressible RANS models for OpenFOAM.com,
targeting OpenFOAM v2412:

- `kOmegaSST2003RC`: Smirnov--Menter rotation/curvature correction on the
  OpenFOAM v2412 SST-2003-like parent model;
- `SpalartAllmarasNoFt2RC`: Spalart--Shur rotation/curvature correction for
  the `SA-noft2` model.

Both implementations include the complete material derivative
`DS/Dt = dS/dt + U.grad(S)`. Inertial calculations and OpenFOAM.com v2412 MRF
zones are implemented. MRF is experimental; SRF, moving-mesh, overset and 6DoF
frame motion are outside the verified scope.

> This offering is not approved or endorsed by OpenCFD Limited, producer and
> distributor of the OpenFOAM software via www.openfoam.com, and owner of the
> OPENFOAM® and OpenCFD® trade marks.

OPENFOAM® is a registered trade mark of OpenCFD Limited, producer and
distributor of the OpenFOAM software via [www.openfoam.com](https://www.openfoam.com/).

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
`system/fvSchemes`:

```text
divSchemes
{
    div(phi,curvatureStrain) Gauss linear;
}
```

SST configuration:

```text
RASModel kOmegaSST2003RC;

kOmegaSST2003RCCoeffs
{
    curvatureCorrection on;
    Cr1 1.0;
    Cr2 2.0;
    Cr3 1.0;
    frMin 0.0;
    frMax 1.25;
    writeCurvatureDiagnostics off;
}
```

SA configuration:

```text
RASModel SpalartAllmarasNoFt2RC;

SpalartAllmarasNoFt2RCCoeffs
{
    curvatureCorrection on;
    Cr1 1.0;
    Cr2 12.0;
    Cr3 1.0;
    writeCurvatureDiagnostics off;
}
```

The explicit RC models enable `curvatureCorrection` by default. The deprecated
aliases `kOmegaSSTCC` and `SpalartAllmarasCC` retain opt-in behavior throughout
the 1.x series. `SpalartAllmarasNoFt2RC` requires `ft2 false`. Both models write
`fr1` and `OmegaRot`; optional diagnostics include `curvatureIndicator`,
`productionBase`, `productionCorrected`, `productionRatio`, and `fr1Clipped`.

## Verification

The mesh-independent oracle checks the tensor algebra, limiting, finite-value
handling, frame invariance and inertial/MRF representation equivalence:

```bash
python3 -m unittest discover -s tests -v
```

With OpenFOAM v2412 loaded, the compiled regressions exercise the material
derivative and MRF integration:

```bash
tests/run_material_derivative_regression.sh
tests/run_mrf_mapping_regression.sh
```

These are implementation-verification tests, not a claim of general physical
validation. See [the implementation audit](docs/implementation-audit.md).

Runnable cases are available in [tutorials/incompressible/simpleFoam](tutorials/incompressible/simpleFoam):
curved channel, flap-tip vortex, Smits curved wall, TMR flat plate, NACA 0012,
TMR rotating pipe and Xiao rotating pipe. Each benchmark directory documents
its setup and run commands. Load OpenFOAM.com v2412 before running `./Allrun`.
The cases use the explicit `kOmegaSST2003RC` and `SpalartAllmarasNoFt2RC`
runtime names; benchmark inputs are fetched separately with
`python3 tools/fetch_benchmark_data.py`. Generated results and
external reference files remain local.

The repository also provides benchmark-data and report-generation utilities
for analysing externally supplied results. The preregistered v1.0 gates and private-runner contract are in
[the validation protocol](docs/validation-protocol.md).

## Repository layout

| Path | Contents |
| --- | --- |
| `src/curvatureCorrectionModels/` | OpenFOAM library: shared RC algebra in `common/`, the SST and SA model classes |
| `tests/` | Python oracle and contract tests; `tests/openfoam/` holds the compiled regressions run by `tests/run_*.sh` |
| `tools/` | Benchmark-data fetcher, validation-data collector, report generator and the `curvature_validation` oracle package |
| `tutorials/incompressible/simpleFoam/` | Runnable cases, each with `Allrun`/`Allclean` and a README |
| `data/` | Pinned benchmark-source manifest and transcribed reference models |
| `docs/` | Audit, benchmark provenance, validation protocol and backlog (see below) |
| `reports/` | Published aggregate validation evidence |

Documentation:

- [Implementation audit](docs/implementation-audit.md): equations, code mapping and literature comparison.
- [Benchmark selection](docs/benchmark-selection.md): candidate cases, provenance and results so far.
- [Validation protocol](docs/validation-protocol.md): preregistered v1.0 gates.
- [Remaining work](docs/remaining-work.md): open issues and acceptance gates.

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
- TMBWG definitions for [SST-RC](https://tmbwg.github.io/turbmodels/sst.html)
  and [SA-RC](https://tmbwg.github.io/turbmodels/spalart.html).

## Licence

The entire repository is licensed under `GPL-3.0-or-later`. See [LICENSE](LICENSE)
and [THIRD_PARTY.md](THIRD_PARTY.md).
