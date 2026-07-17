#include "kOmegaSSTCC.H"
#include "fvcGrad.H"

namespace Foam
{
namespace RASModels
{

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

    if (curvatureCorrection_)
    {
        fr1_ = curvatureCorrection::sstCorrectionFactor
        (
            gradU,
            this->U_,
            this->phi_,
            this->omega_,
            Cr1_,
            Cr2_,
            Cr3_,
            frMin_,
            frMax_
        );
        fr1_.correctBoundaryConditions();
        tProduction.ref() *= fr1_.internalField();
    }
    else
    {
        fr1_ = scalar(1);
    }

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
            true
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
    )
{
    if (frMin_.value() < 0)
    {
        FatalIOErrorInFunction(this->coeffDict_)
            << "frMin must be non-negative" << exit(FatalIOError);
    }

    if (frMin_.value() > frMax_.value())
    {
        FatalIOErrorInFunction(this->coeffDict_)
            << "frMin must not exceed frMax" << exit(FatalIOError);
    }

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
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

    return true;
}

} // End namespace RASModels
} // End namespace Foam
