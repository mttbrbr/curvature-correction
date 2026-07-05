#include "kOmegaSSTCC.H"
#include "fvcGrad.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTCC<BasicTurbulenceModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    // Preserve the production limiter implemented by the standard SST model.
    tmp<volScalarField::Internal> tProduction(Parent::Pk(G));

    if (curvatureCorrection_)
    {
        const tmp<volTensorField> tgradU(fvc::grad(this->U_));
        const tmp<volScalarField> tfr
        (
            curvatureCorrection::correctionFactor
            (
                tgradU(),
                Cr1_,
                Cr2_,
                Cr3_,
                frMin_,
                frMax_
            )
        );

        tProduction.ref() *= tfr().primitiveField();
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
    )
{
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

    if (frMin_.value() > frMax_.value())
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "frMin must not exceed frMax" << exit(FatalIOError);
    }

    return true;
}

} // End namespace RASModels
} // End namespace Foam
