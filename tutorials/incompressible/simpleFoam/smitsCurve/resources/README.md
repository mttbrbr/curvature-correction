# TMBWG source data

These files are unmodified inputs and reference data from the public
[Turbulence Modeling Resource 2DCC case](https://tmbwg.github.io/turbmodels/smitscurve_val.html),
maintained in the [TMBWG turbmodels repository](https://github.com/TMBWG/turbmodels).
They were retrieved on 2026-07-17.

The [official grid page](https://tmbwg.github.io/turbmodels/smitscurve_grids.html)
defines five exactly nested grids. This tutorial retains the three levels used
by its automated workflows:

- `smits_3levdn_129x49.p2dfmt` and `.nmf`: quick validation grid;
- `smits_2levdn_257x97.p2dfmt.gz` and `.nmf`: medium grid;
- `smits_1levdn_513x193.p2dfmt.gz` and `.nmf`: fine grid.

The compressed medium and fine PLOT3D files came from the official archive:
<https://www.nasa.gov/wp-content/uploads/2026/02/smitscurve-grids.zip>.
`Allrun-gridStudy` decompresses them only inside its generated `runs/` tree.

Reference data retained under `reference/` are:

- `smits_cf_cfl3d_sst.dat`: CFL3D SST convex-wall skin friction;
- `smits_cf_cfl3d_sstrc.dat`: CFL3D SST-RC convex-wall skin friction;
- `smits_cf_cfl3d_sa.dat`: CFL3D SA convex-wall skin friction;
- `smits_cf_cfl3d_sarc.dat`: CFL3D SA-RC convex-wall skin friction;
- `Smits_curvature_expdata_reworked2012.dat`: corrected experimental data from
  Smits, Young, and Bradshaw.

The PLOT3D boundary numbering is converted to named OpenFOAM patches by
`system/createPatchDict`. Neutral-map files are retained for traceability.
Verify all files from this directory with `sha256sum -c SHA256SUMS`.
