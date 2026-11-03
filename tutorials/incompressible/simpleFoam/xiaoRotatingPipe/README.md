# Xiao et al. fully developed rotating pipe

This OpenFOAM.com v2412 tutorial evaluates native SST/SA and the local
curvature-corrected variants against the fully developed rotating-pipe DNS of
Xiao, Ceci, Orlandi and Pirozzoli (JFM 996, 2024, A24).  It uses a periodic
full-cylinder O-grid of length `15R`, `Re_b=44000`, constant bulk velocity and
rotation ratios `N=0`, `0.5` and `1`.

The smoke run checks one representative corrected/native pair.  The validation
run executes both turbulence families on three nested grids, adds correction-off
controls on the finest grid, and produces schema-v2 metrics:

```bash
./Allrun smoke
CURVATURE_DNS_RESUME=1 CURVATURE_MPI_PROCS=4 ./Allrun validation
```

DNS accuracy is reported as `improved`, `neutral` or `degraded`; it is not a
model-constant tuning gate.  Evidence completeness, finite/converged solutions,
correction-off recovery and the three-grid study are critical gates.
