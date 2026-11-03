/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2011-2016 OpenFOAM Foundation
 * Copyright (C) 2019-2021 OpenCFD Ltd.
 * Copyright (C) 2026 mttbrbr
 *
 * Derived from the OpenFOAM.com v2412 SpalartAllmaras implementation and
 * modified to add the Spalart-Shur rotation/curvature correction.
 */

#include "SpalartAllmarasCC.H"
#include "bound.H"
#include "fvOptions.H"
#include "fvcGrad.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include <cmath>

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
void SpalartAllmarasCC<BasicTurbulenceModel>::
validateCorrectionCoefficients() const
{
    const dimensionedScalar* coefficients[] = {&Cr1_, &Cr2_, &Cr3_};

    for (const dimensionedScalar* coefficient : coefficients)
    {
        if (!std::isfinite(coefficient->value()))
        {
            FatalIOErrorInFunction(this->coeffDict())
                << coefficient->name() << " must be finite"
                << exit(FatalIOError);
        }
    }
}


template<class BasicTurbulenceModel>
void SpalartAllmarasCC<BasicTurbulenceModel>::
updateDiagnosticWriteOptions()
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
SpalartAllmarasCC<BasicTurbulenceModel>::SpalartAllmarasCC
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
        alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName, type
    ),
    curvatureCorrection_
    (
        Switch::getOrAddToDict
        (
            "curvatureCorrection",
            this->coeffDict_,
            type == word("SpalartAllmarasNoFt2RC")
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
            "Cr2", this->coeffDict_, dimless, 12.0
        )
    ),
    Cr3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cr3", this->coeffDict_, dimless, 1.0
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
    noFt2Variant_(type == word("SpalartAllmarasNoFt2RC")),
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
        dimensionedScalar(sqr(dimLength/dimTime), scalar(0))
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
        dimensionedScalar(sqr(dimLength/dimTime), scalar(0))
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
                << "SpalartAllmarasNoFt2RC. The legacy name will remain "
                << "available through the 1.x release series." << nl;
            warned = true;
        }
    }

    if
    (
        noFt2Variant_
     && this->coeffDict_.template getOrDefault<Switch>("ft2", false)
    )
    {
        FatalIOErrorInFunction(this->coeffDict_)
            << "SpalartAllmarasNoFt2RC requires ft2 false"
            << exit(FatalIOError);
    }

    validateCorrectionCoefficients();
    updateDiagnosticWriteOptions();
    this->printCoeffs(type);
}

template<class BasicTurbulenceModel>
bool SpalartAllmarasCC<BasicTurbulenceModel>::read()
{
    if (!Parent::read())
    {
        return false;
    }

    curvatureCorrection_.readIfPresent
    (
        "curvatureCorrection", this->coeffDict()
    );
    Cr1_.readIfPresent(this->coeffDict());
    Cr2_.readIfPresent(this->coeffDict());
    Cr3_.readIfPresent(this->coeffDict());
    writeCurvatureDiagnostics_.readIfPresent
    (
        "writeCurvatureDiagnostics",
        this->coeffDict()
    );
    if
    (
        noFt2Variant_
     && this->coeffDict().template getOrDefault<Switch>("ft2", false)
    )
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "SpalartAllmarasNoFt2RC requires ft2 false"
            << exit(FatalIOError);
    }
    validateCorrectionCoefficients();
    updateDiagnosticWriteOptions();
    return true;
}

template<class BasicTurbulenceModel>
void SpalartAllmarasCC<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    typedef eddyViscosity<RASModel<BasicTurbulenceModel>> EVM;
    EVM::correct();

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    volScalarField ft2(this->ft2(chi));
    if (noFt2Variant_)
    {
        ft2 = scalar(0);
    }

    const tmp<volTensorField> tgradU(fvc::grad(U));
    const volScalarField dTilda(this->dTilda(chi, fv1, tgradU()));
    const volScalarField Stilda
    (
        this->Stilda(chi, fv1, tgradU(), dTilda)
    );

    if (curvatureCorrection_)
    {
        const bool rotatingFrame =
            curvatureCorrection::updateMRFAngularVelocity(OmegaRot_);

        fr1_ = curvatureCorrection::saCorrectionFactor
        (
            tgradU(),
            U,
            this->phi_,
            OmegaRot_,
            rotatingFrame,
            Cr1_,
            Cr2_,
            Cr3_,
            &curvatureIndicator_,
            &fr1Clipped_
        );
        fr1_.correctBoundaryConditions();
        fr1Clipped_ = scalar(0);
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

    productionBase_ =
        this->Cb1_*Stilda*this->nuTilda_*(scalar(1) - ft2);
    productionCorrected_ =
        this->Cb1_*Stilda*this->nuTilda_*(fr1_ - ft2);
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
    fr1Clipped_.correctBoundaryConditions();

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, this->nuTilda_)
      + fvm::div(alphaRhoPhi, this->nuTilda_)
      - fvm::laplacian
        (
            alpha*rho*this->DnuTildaEff(), this->nuTilda_
        )
      - this->Cb2_/this->sigmaNut_*alpha()*rho()
       *magSqr(fvc::grad(this->nuTilda_)()())
     ==
        this->Cb1_*alpha()*rho()*Stilda()*this->nuTilda_()
       *(fr1_() - ft2())
      - fvm::Sp
        (
            (
                this->Cw1_*this->fw(Stilda, dTilda)
              - this->Cb1_/sqr(this->kappa_)*ft2()
            )
           *alpha()*rho()*this->nuTilda_()/sqr(dTilda()),
            this->nuTilda_
        )
      + fvOptions(alpha, rho, this->nuTilda_)
    );

    nuTildaEqn.ref().relax();
    fvOptions.constrain(nuTildaEqn.ref());
    solve(nuTildaEqn);
    fvOptions.correct(this->nuTilda_);
    bound
    (
        this->nuTilda_,
        dimensionedScalar(this->nuTilda_.dimensions(), Zero)
    );
    this->nuTilda_.correctBoundaryConditions();
    this->correctNut();
}

} // End namespace RASModels
} // End namespace Foam
