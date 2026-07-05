#include "curvatureFunctions.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace curvatureCorrection
{

tmp<volScalarField> correctionFactor
(
    const volTensorField& gradU,
    const dimensionedScalar& Cr1,
    const dimensionedScalar& Cr2,
    const dimensionedScalar& Cr3,
    const dimensionedScalar& frMin,
    const dimensionedScalar& frMax
)
{
    // Strain and rotation tensors both have dimensions [1/s].
    const volSymmTensorField strain(symm(gradU));
    const volTensorField rotation(skew(gradU));

    const volScalarField strainRate(sqrt(2.0*magSqr(strain)));
    const volScalarField rotationRate(sqrt(2.0*magSqr(rotation)));

    // Avoid division by zero in uniform or stagnant regions.
    const dimensionedScalar smallRate
    (
        "smallRate",
        dimless/dimTime,
        SMALL
    );

    const volScalarField rStar
    (
        rotationRate/max(strainRate, smallRate)
    );

    // TODO: Replace this bounded local-invariant placeholder with the selected
    // published correction. A full rotation/curvature model commonly requires
    // a material derivative of the strain tensor and frame-rotation terms.
    const volScalarField frRaw
    (
        (scalar(1) + Cr1)
       *(
            (2.0*rStar/(scalar(1) + sqr(rStar)))
           *(scalar(1) - Cr3*atan(Cr2*(rStar - scalar(1))))
        )
      - Cr1
    );

    return min(max(frRaw, frMin), frMax);
}

} // End namespace curvatureCorrection
} // End namespace Foam
