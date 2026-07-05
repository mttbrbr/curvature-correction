#include "turbulentTransportModels.H"
#include "kOmegaSSTCC.H"

// Instantiate and register kOmegaSSTCC in the incompressible RAS selection
// table used by simpleFoam and other incompressible solvers.
makeRASModel(kOmegaSSTCC);
