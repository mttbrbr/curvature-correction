# Physical-validation protocol

This protocol is preregistered for the v1.0 release. Runnable cases and raw
reference data live on a read-only private archive; only aggregate evidence is
published here.

The mount exposes a private `manifest.json`, a disclosure-safe
`public-manifest.json`, and an executable `bin/run-campaign`. Each campaign runs
native, correction-off and RC variants of both turbulence
families on three systematically refined meshes. The archive manifest records
the immutable case revision, input SHA-256, dataset citation, usage terms,
OpenFOAM release and mesh identifiers. Only the public manifest and the digest
of the private manifest enter the pull request.

The required inertial-frame benchmarks are the Appelbaum--Kloker--Wenzel
curved/flat DNS pairs and the NASA Ames NACA 0012 wingtip-vortex experiment.
Reference profiles with unclear redistribution rights are not committed.

## v1.0 gates

For each model family and benchmark:

- every critical quantity is finite and stable within 0.2% over the final
  averaging window;
- the medium-to-fine change is at most 2%;
- RC improves normalized RMSE by at least 10% for at least two primary
  quantities;
- no critical quantity worsens by more than 5%;
- correction-off agrees with the native parent to normalized L-infinity
  tolerance `1e-10`.

All statuses are computed by the schema-v2 report generator. Missing
observations are `not-verified` and block strict acceptance. MRF remains
experimental and physical MRF validation is not a v1.0 gate.
