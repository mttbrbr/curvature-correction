#include "curvatureFunctions.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

namespace Foam
{
namespace curvatureCorrection
{

namespace
{

tmp<volSymmTensorField> materialDerivativeOfStrain
(
    const volSymmTensorField& strain,
    const volVectorField& U,
    const surfaceScalarField& phi
)
{
    // dS/dt is evaluated from the stored velocity time levels.  For a steady
    // ddt scheme this contribution is identically zero.
    return
        symm(fvc::grad(fvc::ddt(U)))
      + fvc::div(phi, strain, "div(phi,U)")
      - strain*fvc::div(phi);
}

} // End anonymous namespace


tmp<volScalarField> sstCorrectionFactor
(
    const volTensorField& gradU,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& omega,
    const dimensionedScalar& Cr1,
    const dimensionedScalar& Cr2,
    const dimensionedScalar& Cr3,
    const dimensionedScalar& frMin,
    const dimensionedScalar& frMax
)
{
    const volSymmTensorField strain(symm(gradU));
    // OpenFOAM stores grad(U)_ij = dU_j/dx_i; the published W_ij uses dU_i/dx_j.
    const volTensorField rotation(-skew(gradU));
    const volScalarField strainRate(sqrt(2.0*magSqr(strain)));
    const volScalarField rotationRate(sqrt(2.0*magSqr(rotation)));
    const dimensionedScalar smallRate("smallRate", dimless/dimTime, SMALL);

    const volScalarField rStar
    (
        strainRate/max(rotationRate, smallRate)
    );

    // Smirnov-Menter: D^2=max(S^2, 0.09*omega^2).
    const volScalarField D2(max(sqr(strainRate), 0.09*sqr(omega)));
    const volTensorField WS(rotation & strain);
    const volSymmTensorField DSDt(materialDerivativeOfStrain(strain, U, phi));
    const dimensionedScalar smallFourthRate
    (
        "smallFourthRate", pow4(dimless/dimTime), SMALL
    );
    const volScalarField rHat
    (
        (2.0*WS && DSDt)
       /max(rotationRate*D2*sqrt(D2), smallFourthRate)
    );

    const volScalarField fRotation
    (
        (scalar(1) + Cr1)
       *(2.0*rStar/(scalar(1) + rStar))
       *(scalar(1) - Cr3*atan(Cr2*rHat))
      - Cr1
    );

    return min(max(fRotation, frMin), frMax);
}


tmp<volScalarField> saCorrectionFactor
(
    const volTensorField& gradU,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dimensionedScalar& Cr1,
    const dimensionedScalar& Cr2,
    const dimensionedScalar& Cr3
)
{
    const volSymmTensorField strain(symm(gradU));
    // OpenFOAM stores grad(U)_ij = dU_j/dx_i; the published W_ij uses dU_i/dx_j.
    const volTensorField rotation(-skew(gradU));
    const volScalarField strainRate(sqrt(2.0*magSqr(strain)));
    const volScalarField rotationRate(sqrt(2.0*magSqr(rotation)));
    const dimensionedScalar smallRate("smallRate", dimless/dimTime, SMALL);

    const volScalarField rStar
    (
        strainRate/max(rotationRate, smallRate)
    );

    // Spalart-Shur: D^2=(S^2+W^2)/2 and denominator D^4.
    const volScalarField D2
    (
        0.5*(sqr(strainRate) + sqr(rotationRate))
    );
    const volTensorField WS(rotation & strain);
    const volSymmTensorField DSDt(materialDerivativeOfStrain(strain, U, phi));
    const dimensionedScalar smallFourthRate
    (
        "smallFourthRate", pow4(dimless/dimTime), SMALL
    );
    const volScalarField rTilde
    (
        (2.0*WS && DSDt)/max(sqr(D2), smallFourthRate)
    );

    // SA-RC is deliberately not bounded: negative production is part of the
    // published model in strongly stabilising rotation/curvature.
    return
        (scalar(1) + Cr1)
       *(2.0*rStar/(scalar(1) + rStar))
       *(scalar(1) - Cr3*atan(Cr2*rTilde))
      - Cr1;
}

} // End namespace curvatureCorrection
} // End namespace Foam
