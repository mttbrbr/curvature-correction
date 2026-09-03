// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 mttbrbr

#include "turbulentTransportModels.H"
#include "kOmegaSSTCC.H"

// Instantiate and register kOmegaSSTCC in the incompressible RAS selection
// table used by simpleFoam and other incompressible solvers.
makeRASModel(kOmegaSSTCC);
makeRASModel(kOmegaSST2003RC);
