---
description: Accepted framework claim C-CST-008
author: framework-registry
created: '2026-09-05T03:16:39+00:00'
updated: '2026-09-05T03:16:39+00:00'
tags:
- substrate-framework
- accepted-claim
- C-CST-008
category: claims
confidence: established
status: active
---
# C-CST-008

## Statement
The accepted statement is reproduced exactly from the claim registry.

Fix rho>0 and a nonzero curl eigenvalue lambda. Within the declared
stationary finite-variance isotropic Gaussian Beltrami law, there is a
positive-density stationary marked selection of actual invariant knotted
vortex-tube neighborhoods with finite volume and a continuous ambient
complement. On a uniformly bounded good-patch event one can construct
smooth compact isovortical core/reaction fields Q,S wholly inside each
tube, whose full-space induced velocities are also compact, with:

    core_angle(Q)=1, core_angle(S)=0,
    Omega(Q,S)=B!=0,
    H=[[Hq,N],[N,P]] positive definite,
    G(Q)=(B^2/P)n, G(S)=0, L(Q)=0, L(S)=B n.

Here G=rho integral r cross xi and L=rho integral r cross v_xi are
different measured moments, Omega is the Euler KKS form, and H is the
complete isovortical second variation of Euler kinetic energy. Every
quantity is constructed from actual field jets, compact profiles and
Euler integrals. No core mass, locking constant or desired frequency is
supplied. The core angle is a normalized local vorticity-direction
observation. It is not asserted to be a global compact symmetry action.

The same event supports additional zero-G/zero-L compact pairs with
positive reduced stiffness/inertia ratio above any prescribed finite
bound. Their supports can be made disjoint, so the full H/KKS locality
used for their reaction blocks is exact. Finite radii, profile norms,
positive margins, laws and number densities are declared geometric
inputs; no universal five-scalar filament formula is claimed.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: Constant-density incompressible Euler on R^3; real smooth Beltrami samples with a fixed nonzero curl eigenvalue, finite energy density, not finite total energy., EPS tube existence and persistence and the Enciso-Peralta-Salas-Romaniega Gaussian support propositions are declared primary mathematical imports, with archived applicability receipt0060. Exact calculus, compactness/measurable selection, Poisson marking and Gaussian covariance/mixing arguments are explicit mathematical imports., The bounded good-patch geometry, compact cutoffs, finite carriers, density and field/mark law are fixed before deformation. These are geometric ensemble inputs; no target stiffness or frequency is fitted., The normalized local core-vorticity angle is not a global compact group symmetry or an unrestricted Euler invariant mode. Complete finite-core effects are evaluated by actual Euler integrals; no universal filament five-scalar formula is asserted.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.175.0` with provenance `campaigns/P251-cosserat-from-vortex-euler/adjudication.yaml`.

- `proposals/P251-cosserat-from-vortex-euler/attempts/0105/claim-transaction.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0108/C-CST-008-review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0108/shared-scientific-review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0085/compact-spin-theorem.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0085/positive-compact-pair.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0103/six-moment-construction.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0103/joint_moments.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0103/first.stdout`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0102/canonical-coefficient-join.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0098/stationary-compact-assembly.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0071/README.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0060/source-review.md`
- `src/substrate_framework/euler_compact.py`
- `tests/test_euler_compact.py`
- `src/substrate_framework/euler_fourier.py`
- `tests/test_euler_fourier.py`
