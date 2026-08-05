/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2026 mttbrbr
 *
 * Regression application for the OpenFOAM MRF registry adapter.
 */

#include "fvCFD.H"
#include "IOMRFZoneList.H"
#include "curvatureFunctions.H"

#include <cmath>

using namespace Foam;

namespace
{

bool inCellZone
(
    const fvMesh& mesh,
    const word& zoneName,
    const label celli
)
{
    const label zonei = mesh.cellZones().findZoneID(zoneName);

    if (zonei < 0)
    {
        FatalErrorInFunction
            << "Required test cellZone " << zoneName << " is missing"
            << exit(FatalError);
    }

    return mesh.cellZones()[zonei].found(celli);
}


void requireCondition(const bool condition, const string& description)
{
    if (!condition)
    {
        FatalErrorInFunction
            << description.c_str() << exit(FatalError);
    }
}


void requireNear
(
    const string& description,
    const scalar error,
    const scalar tolerance = 1e-12
)
{
    if (!std::isfinite(error) || error > tolerance)
    {
        FatalErrorInFunction
            << description.c_str() << " failed: error=" << error
            << ", tolerance=" << tolerance << exit(FatalError);
    }

    Info<< description.c_str() << ": pass (error=" << error << ')' << nl;
}


scalar boundaryConsistencyError(const volVectorField& field)
{
    scalar error = 0;

    forAll(field.boundaryField(), patchi)
    {
        const vectorField expected
        (
            field.boundaryField()[patchi].patchInternalField()
        );
        const fvPatchVectorField& actual = field.boundaryField()[patchi];

        forAll(actual, facei)
        {
            error = max(error, mag(actual[facei] - expected[facei]));
        }
    }

    return returnReduce(error, maxOp<scalar>());
}


void verifyZero(const volVectorField& field)
{
    scalar error = 0;

    forAll(field, celli)
    {
        error = max(error, mag(field[celli]));
    }

    requireNear("No-MRF zero field", returnReduce(error, maxOp<scalar>()));
    requireNear("No-MRF boundary consistency", boundaryConsistencyError(field));
}


void verifyMapping
(
    const volVectorField& field,
    const vector& positiveOmega,
    const vector& negativeOmega
)
{
    const fvMesh& mesh = field.mesh();
    scalar error = 0;
    label positiveCount = 0;
    label negativeCount = 0;
    label inactiveCount = 0;
    label outsideCount = 0;

    forAll(field, celli)
    {
        vector expected(vector::zero);

        if (inCellZone(mesh, "positive", celli))
        {
            expected = positiveOmega;
            ++positiveCount;
        }
        else if (inCellZone(mesh, "negative", celli))
        {
            expected = negativeOmega;
            ++negativeCount;
        }
        else if (inCellZone(mesh, "inactive", celli))
        {
            ++inactiveCount;
        }
        else
        {
            ++outsideCount;
        }

        error = max(error, mag(field[celli] - expected));
    }

    positiveCount = returnReduce(positiveCount, sumOp<label>());
    negativeCount = returnReduce(negativeCount, sumOp<label>());
    inactiveCount = returnReduce(inactiveCount, sumOp<label>());
    outsideCount = returnReduce(outsideCount, sumOp<label>());

    requireCondition(positiveCount == 1, "positive zone must contain one cell");
    requireCondition(negativeCount == 1, "negative zone must contain one cell");
    requireCondition(inactiveCount == 1, "inactive zone must contain one cell");
    requireCondition(outsideCount == 1, "outside region must contain one cell");
    requireNear("Disjoint-zone mapping", returnReduce(error, maxOp<scalar>()));
    requireNear("Mapping boundary consistency", boundaryConsistencyError(field));
}


void verifyTimedMapping
(
    const volVectorField& field,
    const scalar expectedOmega
)
{
    const fvMesh& mesh = field.mesh();
    scalar error = 0;
    label positiveCount = 0;
    label outsideCount = 0;

    forAll(field, celli)
    {
        vector expected(vector::zero);

        if (inCellZone(mesh, "positive", celli))
        {
            expected.z() = expectedOmega;
            ++positiveCount;
        }
        else
        {
            ++outsideCount;
        }

        error = max(error, mag(field[celli] - expected));
    }

    positiveCount = returnReduce(positiveCount, sumOp<label>());
    outsideCount = returnReduce(outsideCount, sumOp<label>());

    requireCondition(positiveCount == 1, "timed zone must contain one cell");
    requireCondition(outsideCount == 3, "timed outside region must contain three cells");
    requireNear
    (
        "Function1 mapping at omega=" + Foam::name(expectedOmega),
        returnReduce(error, maxOp<scalar>())
    );
    requireNear("Timed boundary consistency", boundaryConsistencyError(field));
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addOption
    (
        "scenario",
        "name",
        "test scenario: no-mrf, mapping, time, overlap, duplicate-registry"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word scenario
    (
        args.getOrDefault<word>("scenario", "mapping")
    );

    IOMRFZoneList MRF(mesh);
    autoPtr<IOMRFZoneList> secondMRF;

    if (scenario == "duplicate-registry")
    {
        secondMRF.reset(new IOMRFZoneList(mesh, "Second"));
    }

    volVectorField OmegaRot
    (
        IOobject
        (
            "OmegaRot",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless/dimTime, vector::zero)
    );

    if (scenario == "no-mrf")
    {
        requireCondition
        (
            !curvatureCorrection::updateMRFAngularVelocity(OmegaRot),
            "empty MRF list must report inactive"
        );
        verifyZero(OmegaRot);
    }
    else if (scenario == "mapping")
    {
        requireCondition
        (
            curvatureCorrection::updateMRFAngularVelocity(OmegaRot),
            "disjoint active MRF zones must report active"
        );
        verifyMapping
        (
            OmegaRot,
            vector(0, 0, 100),
            vector(0, 0, -40)
        );
    }
    else if (scenario == "time")
    {
        const scalar times[] = {0, 1, 2};
        const scalar expected[] = {0, 50, 100};

        for (label i = 0; i < 3; ++i)
        {
            runTime.setTime(instant(times[i]), i);
            requireCondition
            (
                curvatureCorrection::updateMRFAngularVelocity(OmegaRot),
                "timed active MRF zone must report active"
            );
            verifyTimedMapping(OmegaRot, expected[i]);
        }
    }
    else if (scenario == "overlap" || scenario == "duplicate-registry")
    {
        curvatureCorrection::updateMRFAngularVelocity(OmegaRot);
        FatalErrorInFunction
            << "Expected scenario " << scenario << " to fail"
            << exit(FatalError);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown scenario " << scenario << exit(FatalError);
    }

    Info<< "MRF scenario " << scenario << ": pass" << nl;
    Info<< "End" << nl;
    return 0;
}
