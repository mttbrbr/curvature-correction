#ifndef curvatureFunctions_H
#define curvatureFunctions_H

#include "volFields.H"
#include "dimensionedScalar.H"

namespace Foam
{
namespace curvatureCorrection
{

//- Return a bounded, dimensionless curvature/rotation production multiplier.
//
//  NOTE: This is a deliberately simple placeholder based on local invariants.
//  It is isolated here so a validated Spalart-Shur, Smirnov-Menter, or other
//  formulation can replace it without changing the turbulence model class.
tmp<volScalarField> correctionFactor
(
    const volTensorField& gradU,
    const dimensionedScalar& Cr1,
    const dimensionedScalar& Cr2,
    const dimensionedScalar& Cr3,
    const dimensionedScalar& frMin,
    const dimensionedScalar& frMax
);

} // End namespace curvatureCorrection
} // End namespace Foam

#endif
