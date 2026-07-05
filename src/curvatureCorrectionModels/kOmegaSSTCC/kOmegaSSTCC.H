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

    kOmegaSSTCC(const kOmegaSSTCC&) = delete;
    void operator=(const kOmegaSSTCC&) = delete;

protected:

    //- Standard SST production limiter multiplied by bounded fr.
    virtual tmp<volScalarField::Internal> Pk
    (
        const volScalarField::Internal& G
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

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "kOmegaSSTCC.C"
#endif

#endif
