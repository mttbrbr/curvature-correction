#include "fvCFD.H"
#include "curvatureFunctions.H"

using namespace Foam;

namespace
{

void setLinearVelocity(volVectorField& field, const tensor& gradient)
{
    const volVectorField& centres = field.mesh().C();
    vectorField& values = field.primitiveFieldRef();

    forAll(values, celli)
    {
        values[celli] = centres[celli] & gradient;
    }

    forAll(field.boundaryField(), patchi)
    {
        const vectorField patchCentres(centres.boundaryField()[patchi]);
        vectorField patchValues(patchCentres.size(), vector::zero);

        forAll(patchValues, facei)
        {
            patchValues[facei] = patchCentres[facei] & gradient;
        }

        field.boundaryFieldRef()[patchi] == patchValues;
    }
}


void setConstantVelocity(volVectorField& field, const vector& value)
{
    field.primitiveFieldRef() = value;

    forAll(field.boundaryField(), patchi)
    {
        field.boundaryFieldRef()[patchi] == value;
    }
}


void setLinearStrain
(
    volSymmTensorField& field,
    const symmTensor& base,
    const symmTensor& xGradient,
    const symmTensor& yGradient
)
{
    const volVectorField& centres = field.mesh().C();
    symmTensorField& values = field.primitiveFieldRef();

    forAll(values, celli)
    {
        const vector& centre = centres[celli];
        values[celli] =
            base + centre.x()*xGradient + centre.y()*yGradient;
    }

    forAll(field.boundaryField(), patchi)
    {
        const vectorField patchCentres(centres.boundaryField()[patchi]);
        symmTensorField patchValues
        (
            patchCentres.size(),
            symmTensor::zero
        );

        forAll(patchValues, facei)
        {
            const vector& centre = patchCentres[facei];
            patchValues[facei] =
                base + centre.x()*xGradient + centre.y()*yGradient;
        }

        field.boundaryFieldRef()[patchi] == patchValues;
    }
}


scalar maxError
(
    const volSymmTensorField& field,
    const symmTensor& expected
)
{
    scalar error = 0;

    forAll(field, celli)
    {
        error = max(error, mag(field[celli] - expected));
    }

    return returnReduce(error, maxOp<scalar>());
}


scalar maxError(const volScalarField& field, const scalar expected)
{
    scalar error = 0;

    forAll(field, celli)
    {
        error = max(error, mag(field[celli] - expected));
    }

    return returnReduce(error, maxOp<scalar>());
}


void requirePass(const string& label, const scalar error, const scalar tolerance)
{
    if (!std::isfinite(error) || error > tolerance)
    {
        FatalErrorInFunction
            << label.c_str() << " failed: error=" << error
            << ", tolerance=" << tolerance << exit(FatalError);
    }

    Info<< label.c_str() << ": pass (error=" << error << ')' << nl;
}


scalar expectedSstFactor
(
    const tensor& gradU,
    const symmTensor& materialDerivative,
    const scalar omega
)
{
    const symmTensor strain(symm(gradU));
    const tensor rotation(-skew(gradU));
    const scalar strainRate(std::sqrt(2.0*magSqr(strain)));
    const scalar rotationRate(std::sqrt(2.0*magSqr(rotation)));
    const scalar rStar(strainRate/max(rotationRate, scalar(SMALL)));
    const scalar D2(max(sqr(strainRate), 0.09*sqr(omega)));
    const scalar denominator
    (
        max(rotationRate*D2*std::sqrt(D2), scalar(SMALL))
    );
    const scalar indicator
    (
        (2.0*(rotation & strain) && materialDerivative)/denominator
    );
    const scalar raw
    (
        2.0*(2.0*rStar/(1.0 + rStar))
       *(1.0 - std::atan(2.0*indicator))
      - 1.0
    );

    return min(max(raw, scalar(0)), scalar(1.25));
}


scalar expectedSaFactor
(
    const tensor& gradU,
    const symmTensor& materialDerivative
)
{
    const symmTensor strain(symm(gradU));
    const tensor rotation(-skew(gradU));
    const scalar strainRate(std::sqrt(2.0*magSqr(strain)));
    const scalar rotationRate(std::sqrt(2.0*magSqr(rotation)));
    const scalar rStar(strainRate/max(rotationRate, scalar(SMALL)));
    const scalar D2(0.5*(sqr(strainRate) + sqr(rotationRate)));
    const scalar indicator
    (
        (2.0*(rotation & strain) && materialDerivative)
       /max(sqr(D2), scalar(SMALL))
    );

    return
        2.0*(2.0*rStar/(1.0 + rStar))
       *(1.0 - std::atan(12.0*indicator))
      - 1.0;
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "expectSteady",
        "expect the configured ddt scheme to suppress the temporal term"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const bool expectSteady(args.found("expectSteady"));
    const scalar tolerance = 2e-12;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    const tensor currentGradient
    (
        0, 2, 0,
        0, 0, 0,
        0, 0, 0
    );
    const tensor temporalGradient
    (
        1, 0, 0,
        0, -1, 0,
        0, 0, 0
    );
    const tensor oldGradient(currentGradient - temporalGradient);

    setLinearVelocity(U, oldGradient);
    U.oldTime();
    ++runTime;
    setLinearVelocity(U, currentGradient);

    const volTensorField gradU(fvc::grad(U));
    const volSymmTensorField strain(symm(gradU));
    const surfaceScalarField zeroPhi
    (
        IOobject("zeroPhi", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("zero", dimVolume/dimTime, 0)
    );
    const volVectorField OmegaRot
    (
        IOobject("OmegaRot", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("zero", dimless/dimTime, vector::zero)
    );

    const volSymmTensorField temporal
    (
        curvatureCorrection::materialDerivativeOfStrain
        (
            strain, U, zeroPhi, OmegaRot, false
        )
    );
    const symmTensor expectedTemporal
    (
        expectSteady ? symmTensor::zero : symm(temporalGradient)
    );
    requirePass
    (
        expectSteady
      ? "Steady temporal contribution"
      : "Temporal contribution",
        maxError(temporal, expectedTemporal),
        tolerance
    );

    volVectorField convectiveU
    (
        IOobject("convectiveU", runTime.timeName(), mesh),
        U
    );
    convectiveU.clearOldTimes();
    const vector convectionVelocity(0.75, -0.4, 0);
    setConstantVelocity(convectiveU, convectionVelocity);
    convectiveU.oldTime();

    volSymmTensorField convectiveStrain
    (
        IOobject("convectiveStrain", runTime.timeName(), mesh),
        strain
    );
    const symmTensor base(0.3, 0.1, 0, -0.2, 0, 0.4);
    const symmTensor xGradient(1.0, 0.2, 0, -0.3, 0, 0.1);
    const symmTensor yGradient(-0.4, 0.3, 0, 0.8, 0, -0.2);
    setLinearStrain
    (
        convectiveStrain, base, xGradient, yGradient
    );
    const surfaceScalarField convectivePhi(fvc::flux(convectiveU));
    const volSymmTensorField convective
    (
        curvatureCorrection::materialDerivativeOfStrain
        (
            convectiveStrain,
            convectiveU,
            convectivePhi,
            OmegaRot,
            false
        )
    );
    const symmTensor expectedConvective
    (
        convectionVelocity.x()*xGradient
      + convectionVelocity.y()*yGradient
    );
    requirePass
    (
        "Convective contribution",
        maxError(convective, expectedConvective),
        tolerance
    );

    const volScalarField omega
    (
        IOobject("omega", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("omega", dimless/dimTime, 1)
    );
    const dimensionedScalar Cr1("Cr1", dimless, 1);
    const dimensionedScalar sstCr2("Cr2", dimless, 2);
    const dimensionedScalar saCr2("Cr2", dimless, 12);
    const dimensionedScalar Cr3("Cr3", dimless, 1);
    const dimensionedScalar frMin("frMin", dimless, 0);
    const dimensionedScalar frMax("frMax", dimless, 1.25);
    const volScalarField sst
    (
        curvatureCorrection::sstCorrectionFactor
        (
            gradU,
            U,
            zeroPhi,
            OmegaRot,
            false,
            omega,
            Cr1,
            sstCr2,
            Cr3,
            frMin,
            frMax
        )
    );
    const volScalarField sa
    (
        curvatureCorrection::saCorrectionFactor
        (
            gradU,
            U,
            zeroPhi,
            OmegaRot,
            false,
            Cr1,
            saCr2,
            Cr3
        )
    );

    requirePass
    (
        "SST correction",
        maxError
        (
            sst,
            expectedSstFactor(currentGradient, expectedTemporal, 1)
        ),
        tolerance
    );
    requirePass
    (
        "SA correction",
        maxError(sa, expectedSaFactor(currentGradient, expectedTemporal)),
        tolerance
    );

    Info<< "End" << nl;
    return 0;
}
