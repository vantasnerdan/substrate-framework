---
description: Accepted framework claim C-CST-009
author: framework-registry
created: '2026-09-05T03:16:39+00:00'
updated: '2026-09-05T03:16:39+00:00'
tags:
- substrate-framework
- accepted-claim
- C-CST-009
category: claims
confidence: established
status: active
---
# C-CST-009

## Statement
The accepted statement is reproduced exactly from the claim registry.

Use C-CST-008's cells in the explicit phase Cauchy--Born embedding of the
quadratic Euler material action,

    eta=U+Ez,
    pi=rho[(u0.grad)eta+V+A Ez],
    Axi=P(xi cross omega0)-curl(P(xi cross omega0))/lambda.

All ambient fluid is retained. U is the declared coherent tube-centroid
plus continuous-ambient coordinate; Phi is the registered local core-angle
field, q=Phi-curl U/2. The ensemble pairs time-reversed and reflected full
states, with one COMMON macro conjugate variable V and independently varied
microscopic reaction momenta. Unresolved nonaffine variation/relaxation is
excluded by the declared Cauchy--Born closure, not silently made into gauge.
The law and marks are fixed before deformation. The claim concerns the
linear continuum and its full second spatial jet, not all wave numbers.

Pulling back that ONE action, taking its stated ensemble average and
eliminating the complete retained momenta yields positive computed

    j=nu E_Palm[B^2/P]/3,
    kappa=nu E_Palm[Hq-N^2/P]/3, alpha=kappa/4.

The six-moment match gives b=0 in the leading physical mixed mass and
ell=g-kappa b/j=-kappa/2!=0. Thus the optical branch has nonzero actual
coarse translation/core-angle transfer
U/Phi=-j sigma |k|/(2rho)+O(|k|^3) in curl helicity sigma. The leading
mass density is the total rho, derived also by the exact common Galilean
translation of the same stationary fluid, not a tube filling fraction.

Zero-moment compact STF and neighbor-angle attachments, with finite
structurally selected amplitudes, give positive shear mu and positive
transverse/longitudinal spin curvatures C_T,C_L after retaining all added
gradient inertia. Their coefficients are full canonical action/Schur
integrals, including the negative bare material covariance stiffness,
potential reaction squares and the surviving STF rate-source norm.
The SAME derivative field normalization transforms both kinetic and
potential forms. In its stated convention the equations are

    rho U_N,tt=(mu+alpha)Delta U_N+2alpha curl Phi_N-grad p,
    div U_N=0,
    j Phi_N,tt=C_T Delta Phi_N+(C_L-C_T)grad div Phi_N
                    +2alpha curl U_N-4alpha Phi_N.

These equations and their angular-momentum/couple-stress balance are exact
coefficient identities for the declared second-gradient conditional action.
The optical gap is4alpha/j; the acoustic transverse speed squared ismu/rho.
There is no incompressible longitudinal displacement wave. The physical
gradient masses, finite-radius Fourier filters and macro affine tube-spin
row accompany the field map. In particular full tube spin is not renamed
j Phi_dot: its additional affine/shape current is explicitly retained.

Removing the core-angle and angle-gradient populations deletes Phi before
any zero-inertia division and leaves a positive incompressible Navier--Cauchy
sector when the independent STF population remains. Removing every attached
population and the vortical covariance instead leaves linearized Euler.
The explicit limiting families and unchanged-density bookkeeping are in
0105/population-limits.md.

This is an exact conditional variational continuum construction from the
Euler material action. It does NOT claim that unrestricted Euler trajectories
remain in the chosen finite phase family. Their reconstruction residual is
R=Udot-V+Ezdot-AEz;0095 supplies its exact complement/memory and the matching
physical observation correction. Its vanishing is not an unstated premise
or a conclusion of the finite KKS calculation. The stronger unrestricted
realization question remains distinct, with0101 as finite-time progress.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: C-CST-008. Assumptions: All hypotheses of C-CST-008 are retained. The explicit canonical phase Cauchy-Born embedding and fixed coherent marked law are the declared closure; nonaffine relaxation is excluded., The common macro momentum is averaged with independently varied microscopic time-reversal reactions before its elimination. Reflection pairs complete laws and reverses curl eigenvalue; it is not claimed to occur within a single fixed-helicity ergodic sample., The constitutive statement is the exact linear second-spatial-jet action in its long-wave positive-mass neighborhood. Both kinetic and potential forms undergo the same derivative field map, and physical finite-radius/shape observations remain explicit., The full reconstruction residual and free-Euler complement are not assumed zero. No unrestricted trajectory realization, global SO(3) symmetry, all-wave-number dispersion, compressional Euler wave, or old alpha=L_v*T/6 formula is promoted., Exact proof:0097/0102 and the independently reviewed construction/limits. Existing DOP853 integration verifies the reduced equations only; no numerical trajectory is relabeled as an Euler-existence proof.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.175.0` with provenance `campaigns/P251-cosserat-from-vortex-euler/adjudication.yaml`.

- `proposals/P251-cosserat-from-vortex-euler/attempts/0105/claim-transaction.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0105/population-limits.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0105/route-closure-map.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0108/C-CST-009-review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0108/shared-scientific-review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0097/canonical-phase-pullback.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0097/material-cotangent-bridge.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0097/verify.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0097/first-canonical.stdout`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0102/canonical-coefficient-join.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0102/verify.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0102/stdout.txt`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0098/stationary-compact-assembly.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0107/README.md`
- `proposals/P251-cosserat-from-vortex-euler/verify_cst004.py`
- `proposals/P251-cosserat-from-vortex-euler/verify_cst005.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0056/cst005-first.stdout.txt`
- `src/substrate_framework/euler_orbit.py`
- `tests/test_euler_orbit.py`
- `src/substrate_framework/euler_fourier.py`
- `src/substrate_framework/micropolar.py`
- `tests/test_micropolar.py`
