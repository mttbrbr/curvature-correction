# Curved channel (U-bend smoke case)

A minimal two-dimensional 180-degree U-bend used to check that
`kOmegaSST2003RC` builds, runs inside `simpleFoam` and writes a bounded `fr1`
field. It is a smoke test, not a validation case: there is no reference data
and the mesh (912 cells) is far too coarse for quantitative conclusions.

| Quantity | Value |
| --- | --- |
| Channel width `h` | 50 mm |
| Inner / outer bend radius | 25 mm / 75 mm (`h/Rc = 1`) |
| Inlet velocity | 10 m/s (uniform) |
| Kinematic viscosity | 1.5e-5 m²/s (`Re_h ≈ 3.3e4`) |
| Iterations | 100 steady `simpleFoam` iterations |

Run after loading OpenFOAM.com v2412 and building the library:

```bash
./Allrun
./Allclean
```

Expected outcome: the run completes without warnings, and `100/fr1` is finite
and lies in `[frMin, frMax] = [0, 1.25]`. In the bend the correction reduces
production near the convex (inner) wall and increases it near the concave
(outer) wall; on OpenFOAM.com v2412 the mean `fr1` in the inner and outer
quarter of the bend is about 0.14 and 1.25 (clipped), respectively.
