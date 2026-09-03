/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2011-2016 OpenFOAM Foundation
 * Copyright (C) 2019 OpenCFD Ltd.
 * Copyright (C) 2026 mttbrbr
 *
 * Derived from the OpenFOAM.com v2412 kOmegaSST interface and modified to
 * expose the Smirnov-Menter rotation/curvature correction.
 */

#ifndef kOmegaSSTCC_H
#define kOmegaSSTCC_H

#include "kOmegaSST.H"
#include "curvatureFunctions.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
class kOmegaSSTCC
:
    public kOmegaSST<BasicTurbulenceModel>
{
    typedef kOmegaSST<BasicTurbulenceModel> Parent;

    // Model controls. Dimensioned scalars enforce dimensionless dictionary
    // input and make accidental dimensional coefficients fail early.
    Switch curvatureCorrection_;
    dimensionedScalar Cr1_;
    dimensionedScalar Cr2_;
    dimensionedScalar Cr3_;
    dimensionedScalar frMin_;
    dimensionedScalar frMax_;
    Switch writeCurvatureDiagnostics_;

    // GbyNu0() is const, therefore cached diagnostics must be mutable.
    mutable volScalarField fr1_;
    mutable volVectorField OmegaRot_;
    mutable volScalarField curvatureIndicator_;
    mutable volScalarField productionBase_;
    mutable volScalarField productionCorrected_;
    mutable volScalarField productionRatio_;
    mutable volScalarField fr1Clipped_;

    kOmegaSSTCC(const kOmegaSSTCC&) = delete;
    void operator=(const kOmegaSSTCC&) = delete;

    void validateCorrectionCoefficients() const;
    void updateDiagnosticWriteOptions();

protected:

    //- Apply SST-RC to the raw production before either SST limiter is used.
    virtual tmp<volScalarField::Internal> GbyNu0
    (
        const volTensorField& gradU,
        const volScalarField& S2
    ) const override;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    TypeName("kOmegaSSTCC");

    kOmegaSSTCC
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

    virtual ~kOmegaSSTCC() = default;

    virtual bool read() override;
};


//- Scientifically explicit runtime name for the OpenFOAM v2412 SST parent.
template<class BasicTurbulenceModel>
class kOmegaSST2003RC
:
    public kOmegaSSTCC<BasicTurbulenceModel>
{
public:
    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    TypeName("kOmegaSST2003RC");

    kOmegaSST2003RC
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    )
    :
        kOmegaSSTCC<BasicTurbulenceModel>
        (
            alpha, rho, U, alphaRhoPhi, phi, transport,
            propertiesName, type
        )
    {}
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "kOmegaSSTCC.C"
#endif

#endif
