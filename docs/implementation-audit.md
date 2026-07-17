# Rotation/curvature-correction implementation audit

Audit target: OpenFOAM.com v2412, incompressible RANS library in
`src/curvatureCorrectionModels`. The parent implementations used for the
comparison are the locally installed v2412 `kOmegaSSTBase` and
`SpalartAllmarasBase`. This document distinguishes code verification from model
validation; agreement with a published equation does not establish predictive
accuracy.

## Common kinematics

The code forms

\[
S_{ij}=\tfrac12(U_{i,j}+U_{j,i}),\qquad
W_{ij}=\tfrac12(U_{i,j}-U_{j,i}),
\]

and the magnitudes

\[
S=\sqrt{2S_{ij}S_{ij}},\qquad W=\sqrt{2W_{ij}W_{ij}}.
\]

OpenFOAM stores `grad(U)_ij = dU_j/dx_i`; consequently the implementation uses
`-skew(gradU)` for the published `W`. The sign and index conversion are
consistent. The Lagrangian derivative is discretised as

\[
\frac{DS}{Dt}=\operatorname{sym}\nabla\left(\frac{\partial U}{\partial t}\right)
 + \nabla\cdot(\phi S)-S\,\nabla\cdot\phi .
\]

For a volumetric velocity flux the final two terms recover `U.grad(S)`, also on
a non-solenoidal discrete field. For `steadyState` the first term is zero. The
implementation is in `common/curvatureFunctions.C:14-27`; tensor construction
and invariants are at lines 45-69 (SST) and 93-119 (SA).

All terms are dimensionally consistent. `S`, `W`, and `omega` have dimensions
T^-1; `DS/Dt` is T^-2; both invariant numerators and denominators are T^-4;
`rStar`, the curvature indicator, and `fr1` are dimensionless.

The contractions also match the published index order:
`(rotation & strain) && DSDt` is
`W_ik S_kj DS_ij/Dt = W_ik S_jk DS_ij/Dt`, since `S` is symmetric.

## SST correction

The implemented Smirnov--Menter function is

\[
r^*=S/W,
\quad D^2=\max(S^2,0.09\omega^2),
\]

\[
\hat r=\frac{2W_{ik}S_{jk}(DS_{ij}/Dt)}{W D^3},
\]

\[
f_{rotation}=(1+C_{r1})\frac{2r^*}{1+r^*}
 [1-C_{r3}\tan^{-1}(C_{r2}\hat r)]-C_{r1},
\]

\[
f_{r1}=\min(\max(f_{rotation},f_{min}),f_{max}).
\]

Defaults are `Cr1=1`, `Cr2=2`, `Cr3=1`, `frMin=0`, and `frMax=1.25`
(`kOmegaSSTCC.C:71-113`). Only `frMin >= 0` and `frMin <= frMax` are validated.

The override calls the native raw production and multiplies it by `fr1`
(`kOmegaSSTCC.C:10-43`). In the v2412 parent, this corrected raw production is
used to form `G=nut*GbyNu0`; the omega production limiter is then applied to
`GbyNu0`, and the k production limiter to `G`. Thus both transport equations
receive the correction, before their native limiters. `k`, `omega`, and `nut`
remain bounded/corrected by the parent implementation.

There is a naming/documentation nuance: the OpenFOAM v2412 parent uses the
strain-based eddy-viscosity denominator, `gamma1=5/9`, `gamma2=0.44`, and a
factor-10 production limiter. These are SST-2003 characteristics. NASA TMR says
a correction attached to that parent should be identified as SST-2003RC (with
the appropriate `m/s/e` suffix), rather than unqualified standard SST-RC. The
correction function itself matches the NASA TMR SST-RC equation.

## Spalart--Allmaras correction

The implemented Shur--Spalart function is

\[
r^*=S/W,\qquad D^2=\tfrac12(S^2+W^2),
\]

\[
\tilde r=\frac{2W_{ik}S_{jk}(DS_{ij}/Dt)}{D^4},
\]

with the same functional form for `fr1`, using `Cr1=1`, `Cr2=12`, and `Cr3=1`.
In accordance with SA-RC, it is not clipped (`curvatureFunctions.C:121-127`).

The native production
`Cb1*(1-ft2)*Stilda*nuTilda` becomes
`Cb1*(fr1-ft2)*Stilda*nuTilda`
(`SpalartAllmarasCC.C:139-160`). Destruction, diffusion, `Stilda`, wall distance,
and `fv1` are unchanged. The solved working variable `nuTilda` is bounded at
zero and `nut=nuTilda*fv1` is recomputed (`SpalartAllmarasCC.C:164-174`). The
v2412 parent defaults `ft2=false`, so the usual effective production multiplier
is simply `fr1` unless the user opts into `ft2`.

## Correction-specific constants

| Model | Entry | Default | Dimension | Validation |
| --- | --- | ---: | --- | --- |
| both | `curvatureCorrection` | `true` | switch | parsed by OpenFOAM |
| SST | `Cr1` | 1.0 | dimensionless | no finite/range check |
| SST | `Cr2` | 2.0 | dimensionless | no finite/range check |
| SST | `Cr3` | 1.0 | dimensionless | no finite/range check |
| SST | `frMin` | 0.0 | dimensionless | non-negative, <= `frMax` |
| SST | `frMax` | 1.25 | dimensionless | >= `frMin` |
| SA | `Cr1` | 1.0 | dimensionless | no finite/range check |
| SA | `Cr2` | 12.0 | dimensionless | no finite/range check |
| SA | `Cr3` | 1.0 | dimensionless | no finite/range check |

Parent-model constants are inherited without modification. Important defaults
for traceability are SST `betaStar=0.09`, `a1=0.31`, `b1=1`, `c1=10`, and SA
`Cb1=0.1355`, `Cb2=0.622`, `kappa=0.41`, `Cv1=7.1`, `Cs=0.3`, `ft2=false`.

## Limits, continuity, and numerical risks

| State/limit | SST result | SA result | Interpretation |
| --- | ---: | ---: | --- |
| simple shear, `DS/Dt=0`, `S=W` | 1 | 1 | neutral |
| `S=W=0` after regularisation | 0 | -1 | production is still zero, but factor is not neutral |
| `S->0`, `W>0` | 0 | -1 | pure rotation has zero raw strain production |
| `W->0`, `S>0` | 1.25 | 3 | SST upper clip active in pure strain |
| indicator -> +infinity | clipped to [0,1.25] | can be negative | stabilising branch |
| indicator -> -infinity | clipped to [0,1.25] | can be large positive | destabilising branch |

The algebra is continuous for finite inputs because `max`, `min`, and `atan`
are continuous, but it is not differentiable at clipping and denominator-switch
surfaces. The `SMALL` floors prevent division by exactly zero with consistent
dimensions. They are extremely small, so near-degenerate finite fields can
still produce a very large indicator before `atan` saturates.

The principal risks are:

1. There is no explicit finite-value validation for velocity gradients,
   `omega`, coefficients, intermediate invariants, or `fr1`. NaN/Inf can
   propagate; the deterministic oracle rejects them, but production C++ does
   not fail at the point of origin.
2. SA is deliberately unbounded. The Smits case observed a broad range of
   `fr1`; stability and clipping-equivalent diagnostics must be assessed per
   application. Negative production is a model feature, not automatically a
   code defect.
3. The factor in uniform flow is 0 (SST) or -1 (SA), not 1. Since base
   production is zero, the production response is neutral. Consumers must not
   interpret `fr1` alone as a valid neutral-flow diagnostic.
4. SST clips at 1.25 in pure strain, so clipping is not restricted to pathological
   inputs. A volume/production-weighted clipping fraction should be reported;
   min/max alone cannot establish whether clipping is dominant.
5. The material derivative is scheme- and mesh-dependent. No manufactured
   solution presently quantifies its discretisation error.
6. Both models write `fr1`, but do not write the raw indicator, base/corrected
   production, or a clipping mask. These requested diagnostics require an
   additional function object or model fields.

## Frame invariance

The scalar contractions are invariant under a constant orthogonal change of
Cartesian basis; the deterministic tests verify this. The published equations
for a rotating calculation frame additionally modify `W_ij` with
`Omega_rot` and add the corresponding `Omega_rot` term to `DS_ij/Dt`. Neither
term exists here. Therefore the implementation is complete only in an inertial
frame and must not be claimed frame-invariant for MRF/rotating-frame use.
Consequently a rotating-plane-channel test performed in a rotating reference
frame would test a known-incomplete equation and has not been treated as valid
evidence.

## Literature comparison and disposition

- SA formula, constants, placement on `Cb1` production, and absence of clipping
  agree with the NASA Turbulence Modeling Resource SA-RC definition and Shur et
  al. (2000).
- SST formula, constants, `D`, denominator, clipping, and multiplication of
  both production terms agree with the NASA TMR SST-RC definition and
  Smirnov--Menter (2009).
- The rotating-frame terms explicitly present in both published definitions are
  omitted. This is documented in the repository and is a material scope limit.
- The SST correction is attached to OpenFOAM's SST-2003-like parent. The model
  name/documentation should identify this base variant precisely.

Primary references:

- <https://tmbwg.github.io/turbmodels/spalart.html>
- <https://tmbwg.github.io/turbmodels/sst.html>
- <https://doi.org/10.2514/2.1058>
- <https://doi.org/10.1115/1.3070573>

The implementation is formula-consistent within its stated inertial-frame
scope, but the audit alone cannot support merge readiness. Rotating-frame
support, neutral-flow regression evidence, independent concave/rotating
benchmarks, and broader validation remain separate acceptance gates.
