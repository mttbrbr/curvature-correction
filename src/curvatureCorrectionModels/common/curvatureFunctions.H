/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2026 mttbrbr
 *
 * OpenFOAM adapter interface for rotation/curvature invariants.
 */

#ifndef curvatureFunctions_H
#define curvatureFunctions_H

#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"

namespace Foam
{
namespace curvatureCorrection
{

// Material derivative used by both correction families. Exposed so the
// discretised temporal and convective contributions can be regression-tested.
tmp<volSymmTensorField> materialDerivativeOfStrain
(
    const volSymmTensorField& strain,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volVectorField& OmegaRot,
    const bool rotatingFrame
);

// Populate the per-cell angular velocity of any registered OpenFOAM MRF
// zones. Returns true when at least one active MRF zone is present.
bool updateMRFAngularVelocity(volVectorField& OmegaRot);

// Smirnov-Menter SST-RC correction (J. Turbomach. 131, 041010, 2009).
tmp<volScalarField> sstCorrectionFactor
(
    const volTensorField& gradU,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volVectorField& OmegaRot,
    const bool rotatingFrame,
    const volScalarField& omega,
    const dimensionedScalar& Cr1,
    const dimensionedScalar& Cr2,
    const dimensionedScalar& Cr3,
    const dimensionedScalar& frMin,
    const dimensionedScalar& frMax,
    volScalarField* indicator = nullptr,
    volScalarField* rawFactor = nullptr
);

// Spalart-Shur SA-RC correction (AIAA J. 38(5), 784-792, 2000).
tmp<volScalarField> saCorrectionFactor
(
    const volTensorField& gradU,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volVectorField& OmegaRot,
    const bool rotatingFrame,
    const dimensionedScalar& Cr1,
    const dimensionedScalar& Cr2,
    const dimensionedScalar& Cr3,
    volScalarField* indicator = nullptr,
    volScalarField* rawFactor = nullptr
);

} // End namespace curvatureCorrection
} // End namespace Foam

#endif
