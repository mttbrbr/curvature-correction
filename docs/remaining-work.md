# Remaining implementation and validation work

Status date: 2026-08-05. This backlog is issue-ready because the repository
GitHub Project cannot be updated with the currently available credentials.

## Completed in the current implementation

- Published SST-RC and SA-RC algebra, transient material derivative, and
  OpenFOAM.com v2412 MRF basis term.
- Explicit opt-in correction (`curvatureCorrection false` by default), finite
  coefficient/invariant guards, and optional production/clipping diagnostics.
- Python kinematic/source contracts, compiled transient/steady regression, and
  SST serial/MPI plus SA serial MRF smoke coverage.
- Pinned checksum-verifying external-data fetch, hold policy, mixed MIT/GPL
  licensing, strict report mode, and automatic quick CI for trusted changes.

## Proposed GitHub issues

| Priority | Issue title | Missing evidence or implementation | Acceptance gate |
| --- | --- | --- | --- |
| P0 | Validate RC response in a physical rotating-flow benchmark | MRF algebra and runtime mapping are verified, but no measured/DNS rotating-channel or rotating-pipe response is accepted | At least two rotation numbers plus the non-rotating baseline; profiles and uncertainty; inertial/MRF agreement within declared numerical tolerance |
| P0 | Add MRF topology and time-dependence integration matrix | Runtime test currently covers one active constant-speed zone; inactive, multiple disjoint zones, overlap rejection, and time-varying `Function1` omega lack runtime tests | Automated cases verify zero outside zones, correct value in every disjoint zone, inactive-zone exclusion, fatal overlap, and updated omega over time |
| P1 | Add an independent concave-curvature benchmark | Current quantitative evidence is neutral flat plate and convex/airfoil response | Licensed reference with mesh study; both correction families; declared trend and error gates |
| P1 | Resolve separated-flow and NACA 0012 alpha=15 evidence | The 225x65 branch is under-resolved and diagnostic only | 449x129-or-finer sequence, scheme study, converged forces/profiles, and comparison with independent computational references |
| P1 | Complete numerical and model-constant sensitivity matrix | Scheme/mesh sensitivity beyond current cases, Reynolds/rotation number, `delta/Rc`, RC constants, clips, and wall resolution are incomplete | Versioned matrix of ranges and aggregate metrics; no non-finite fields; conclusions bounded to sampled domain |
| P2 | Support additional rotating/moving frame formulations | SRF, dynamic mesh, overset, and 6DoF are outside the implemented/validated scope | One representation-equivalence derivation and runtime regression per supported formulation |
| P2 | Add attribution-aware DNS importer if DNS data is adopted | Appelbaum--Kloker--Wenzel CC BY data remains external by design | Version-pinned importer, checksum, attribution/change record, ignored raw/derived data, and no redistribution by default |
| P2 | Resolve held Smits archive provenance | Medium/fine archive and experimental series have no verified file-level redistribution clearance | Written licence evidence or continued user-supplied-only workflow; never weaken the existing hold by inference |

The first two issues are the remaining engineering gates for a strong claim of
MRF support. The next three are validation gates for broader predictive claims.
The final issues are optional scope/licensing work and do not block the current
experimental, default-off release.
