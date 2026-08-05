/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2011-2016 OpenFOAM Foundation
 * Copyright (C) 2019-2021 OpenCFD Ltd.
 * Copyright (C) 2026 mttbrbr
 *
 * Derived from the OpenFOAM.com v2412 SpalartAllmaras interface and modified
 * to expose the Spalart-Shur rotation/curvature correction.
 */

#ifndef SpalartAllmarasCC_H
#define SpalartAllmarasCC_H

#include "SpalartAllmaras.H"
#include "curvatureFunctions.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
class SpalartAllmarasCC
:
    public SpalartAllmaras<BasicTurbulenceModel>
{
    typedef SpalartAllmaras<BasicTurbulenceModel> Parent;

    Switch curvatureCorrection_;
    dimensionedScalar Cr1_;
    dimensionedScalar Cr2_;
    dimensionedScalar Cr3_;
    Switch writeCurvatureDiagnostics_;
    volScalarField fr1_;
    volVectorField OmegaRot_;
    volScalarField curvatureIndicator_;
    volScalarField productionBase_;
    volScalarField productionCorrected_;
    volScalarField productionRatio_;
    volScalarField fr1Clipped_;

    SpalartAllmarasCC(const SpalartAllmarasCC&) = delete;
    void operator=(const SpalartAllmarasCC&) = delete;

    void validateCorrectionCoefficients() const;
    void updateDiagnosticWriteOptions();

public:
    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    TypeName("SpalartAllmarasCC");

    SpalartAllmarasCC
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    );

    virtual ~SpalartAllmarasCC() = default;
    virtual bool read() override;
    virtual void correct() override;
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "SpalartAllmarasCC.C"
#endif

#endif
