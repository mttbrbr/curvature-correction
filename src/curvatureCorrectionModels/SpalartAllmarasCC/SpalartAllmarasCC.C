#include "SpalartAllmarasCC.H"
#include "bound.H"
#include "fvOptions.H"
#include "fvcGrad.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"

namespace Foam
{
namespace RASModels
{

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
            "curvatureCorrection", this->coeffDict_, true
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
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
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
    const volScalarField ft2(this->ft2(chi));

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
            Cr3_
        );
        fr1_.correctBoundaryConditions();
    }
    else
    {
        fr1_ = scalar(1);
        OmegaRot_ = dimensionedVector
        (
            "zeroOmegaRot", dimless/dimTime, vector::zero
        );
    }

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
