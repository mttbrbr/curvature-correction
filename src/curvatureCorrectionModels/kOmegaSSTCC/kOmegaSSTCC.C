/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2011-2016 OpenFOAM Foundation
 * Copyright (C) 2019 OpenCFD Ltd.
 * Copyright (C) 2026 mttbrbr
 *
 * Derived from the OpenFOAM.com v2412 kOmegaSST implementation and modified
 * to add the Smirnov-Menter rotation/curvature correction.
 */

#include "kOmegaSSTCC.H"
#include "fvcGrad.H"
#include <cmath>

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
void kOmegaSSTCC<BasicTurbulenceModel>::validateCorrectionCoefficients() const
{
    const dimensionedScalar* coefficients[] =
    {
        &Cr1_, &Cr2_, &Cr3_, &frMin_, &frMax_
    };

    for (const dimensionedScalar* coefficient : coefficients)
    {
        if (!std::isfinite(coefficient->value()))
        {
            FatalIOErrorInFunction(this->coeffDict())
                << coefficient->name() << " must be finite"
                << exit(FatalIOError);
        }
    }

    if (frMin_.value() < 0)
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "frMin must be non-negative" << exit(FatalIOError);
    }

    if (frMin_.value() > frMax_.value())
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "frMin must not exceed frMax" << exit(FatalIOError);
    }
}


template<class BasicTurbulenceModel>
void kOmegaSSTCC<BasicTurbulenceModel>::updateDiagnosticWriteOptions()
{
    const IOobject::writeOption option
    (
        writeCurvatureDiagnostics_
      ? IOobject::AUTO_WRITE
      : IOobject::NO_WRITE
    );

    curvatureIndicator_.writeOpt(option);
    productionBase_.writeOpt(option);
    productionCorrected_.writeOpt(option);
    productionRatio_.writeOpt(option);
    fr1Clipped_.writeOpt(option);
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTCC<BasicTurbulenceModel>::GbyNu0
(
    const volTensorField& gradU,
    const volScalarField& S2
) const
{
    tmp<volScalarField::Internal> tProduction
    (
        Parent::GbyNu0(gradU, S2)
    );
    productionBase_.primitiveFieldRef() = tProduction();

    if (curvatureCorrection_)
    {
        const bool rotatingFrame =
            curvatureCorrection::updateMRFAngularVelocity(OmegaRot_);

        fr1_ = curvatureCorrection::sstCorrectionFactor
        (
            gradU,
            this->U_,
            this->phi_,
            OmegaRot_,
            rotatingFrame,
            this->omega_,
            Cr1_,
            Cr2_,
            Cr3_,
            frMin_,
            frMax_,
            &curvatureIndicator_,
            &fr1Clipped_
        );
        fr1_.correctBoundaryConditions();

        forAll(fr1Clipped_, celli)
        {
            const scalar raw = fr1Clipped_[celli];
            fr1Clipped_[celli] =
                raw < frMin_.value() || raw > frMax_.value();
        }
        fr1Clipped_.correctBoundaryConditions();
        tProduction.ref() *= fr1_.internalField();
    }
    else
    {
        fr1_ = scalar(1);
        OmegaRot_ = dimensionedVector
        (
            "zeroOmegaRot", dimless/dimTime, vector::zero
        );
        curvatureIndicator_ = scalar(0);
        fr1Clipped_ = scalar(0);
    }

    productionCorrected_.primitiveFieldRef() = tProduction();
    productionRatio_ = scalar(1);

    forAll(productionRatio_, celli)
    {
        const scalar base = productionBase_[celli];

        if (mag(base) > SMALL)
        {
            productionRatio_[celli] = productionCorrected_[celli]/base;
        }
    }

    productionBase_.correctBoundaryConditions();
    productionCorrected_.correctBoundaryConditions();
    productionRatio_.correctBoundaryConditions();

    return tProduction;
}


template<class BasicTurbulenceModel>
kOmegaSSTCC<BasicTurbulenceModel>::kOmegaSSTCC
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    Parent
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),
    curvatureCorrection_
    (
        Switch::getOrAddToDict
        (
            "curvatureCorrection",
            this->coeffDict_,
            type == word("kOmegaSST2003RC")
        )
    ),
    Cr1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cr1", this->coeffDict_, dimless, 1.0
        )
    ),
    Cr2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cr2", this->coeffDict_, dimless, 2.0
        )
    ),
    Cr3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cr3", this->coeffDict_, dimless, 1.0
        )
    ),
    frMin_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "frMin", this->coeffDict_, dimless, 0.0
        )
    ),
    frMax_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "frMax", this->coeffDict_, dimless, 1.25
        )
    ),
    writeCurvatureDiagnostics_
    (
        Switch::getOrAddToDict
        (
            "writeCurvatureDiagnostics",
            this->coeffDict_,
            false
        )
    ),
    fr1_
    (
        IOobject
        (
            "fr1",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, scalar(1))
    ),
    OmegaRot_
    (
        IOobject
        (
            "OmegaRot",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedVector(dimless/dimTime, vector::zero)
    ),
    curvatureIndicator_
    (
        IOobject
        (
            "curvatureIndicator",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, scalar(0))
    ),
    productionBase_
    (
        IOobject
        (
            "productionBase",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(sqr(dimless/dimTime), scalar(0))
    ),
    productionCorrected_
    (
        IOobject
        (
            "productionCorrected",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(sqr(dimless/dimTime), scalar(0))
    ),
    productionRatio_
    (
        IOobject
        (
            "productionRatio",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, scalar(1))
    ),
    fr1Clipped_
    (
        IOobject
        (
            "fr1Clipped",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, scalar(0))
    )
{
    if (type == typeName)
    {
        static bool warned = false;
        if (!warned)
        {
            WarningInFunction
                << "Runtime model " << typeName << " is deprecated; use "
                << "kOmegaSST2003RC. The legacy name will remain available "
                << "through the 1.x release series." << nl;
            warned = true;
        }
    }

    validateCorrectionCoefficients();
    updateDiagnosticWriteOptions();
    this->printCoeffs(type);
}


template<class BasicTurbulenceModel>
bool kOmegaSSTCC<BasicTurbulenceModel>::read()
{
    if (!Parent::read())
    {
        return false;
    }

    curvatureCorrection_.readIfPresent
    (
        "curvatureCorrection",
        this->coeffDict()
    );
    Cr1_.readIfPresent(this->coeffDict());
    Cr2_.readIfPresent(this->coeffDict());
    Cr3_.readIfPresent(this->coeffDict());
    frMin_.readIfPresent(this->coeffDict());
    frMax_.readIfPresent(this->coeffDict());
    writeCurvatureDiagnostics_.readIfPresent
    (
        "writeCurvatureDiagnostics",
        this->coeffDict()
    );
    validateCorrectionCoefficients();
    updateDiagnosticWriteOptions();

    return true;
}

} // End namespace RASModels
} // End namespace Foam
