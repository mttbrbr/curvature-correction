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
    volScalarField fr1_;

    SpalartAllmarasCC(const SpalartAllmarasCC&) = delete;
    void operator=(const SpalartAllmarasCC&) = delete;

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
