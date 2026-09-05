---
description: Accepted framework claim C-CST-012
author: framework-registry
created: '2026-09-05T09:08:47+00:00'
updated: '2026-09-05T09:08:47+00:00'
tags:
- substrate-framework
- accepted-claim
- C-CST-012
category: claims
confidence: established
status: active
---
# C-CST-012

## Statement
The accepted statement is reproduced exactly from the claim registry.

For positive Psi, lambda and density rho, let
psi=Psi sum_j cos(b_j.x), where the three planar wavevectors have
length lambda and sum to zero. On their triangular periodic cell,
u=(J grad psi,0) is an actual smooth stationary Euler array with
full-cell velocity covariance C_v=3 Psi^2 lambda^2 I/4. Its
zero-axial-wavevector horizontal linear Euler group is uniformly
bounded in L2 vorticity on the real vector isotypic sector under
sixty-degree rotations, including its translation kernel.

The entire first-shell Arnold complement has lower bound 2/3.
Bounded centered coordinates on the actual invariant separatrix
polygons give conserved adjoint translation rows with nonzero
matrix lambda^2*(-Psi)*[[0,1],[-1,0]]. Together these control
the complete excited group, without a transport-frequency gap.

For axial Bloch k approaching zero, actual common-velocity V0
and Kelvin-prepared displacement X0 Euler histories satisfy

  m(t)=cos(|k| c_b t)V0-|k| c_b sin(|k| c_b t)X0
       +O_T(|k|/lambda)(|V0|+|k X0|),
  c_b=sqrt(3)*Psi*lambda/2,

uniformly on every fixed slow window 0<=t<=T/|k|, after the
declared background nondimensionalization, with the corresponding
matrix-sine estimate for |k|X. Here m is the complete Euler
horizontal velocity mean and X_t=m, X(0)=X0. Full pressure,
horizontal return, vertical flow and separatrix contributions
remain. The compensated current p_c=m+ik<r Z> differs from m
by the controlled remainder; its derivative, not a derivative
inferred from uniform closeness of m, licenses the slow action.

The same material action has exact initial symplectic pairing
rho X0.V0, zero displacement-displacement pairing, and initial
displacement stiffness rho k^2 lambda^2 C_v/(lambda^2+k^2).
The complete moving action in the actual compensated-current
slow chart has leading mass rho and stiffness rho k^2 C_v,
with O(|k|/lambda) normalized corrections. Common-V has its
specified initial circulation variation; prepared-D is fixed-
Kelvin. Each trajectory transports its own circulation data.

This is an analytic actual-Euler acoustic-window theorem, not a
fitted oscillator, sampled spectrum or assumed group bound.
Nonlinear amplitude can be chosen last on each finite window;
no k-uniform nonlinear amplitude or infinite-time nonlinear
stability is asserted. The integrated velocity X is not silently
identified with a marked material centroid. Constant-curl lifting,
nonaxial wavevectors, isotropic optical joining, EPS topology and
the full autonomous Euler-to-Cosserat continuum remain distinct
constructions, not conclusions of this planar theorem.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: Smooth constant-density incompressible Euler on the specified full triangular cell, with axial Bloch perturbations and complete pressure. Positive Psi, lambda, rho are fixed before the small axial wavevector limit., The exact vector isotypic sector is selected by common horizontal data and preserved by the complete planar-background Euler operator. No conclusion about other representations or lifted backgrounds is imported., Standard periodic Fourier analysis, smooth transport groups, bounded perturbation and finite-time Euler/Lin well-posedness are declared mathematical imports. The actual Fourier gap, weak transport-domain primitives, dual conservation and acoustic estimate are derived in0161 and0146., Analytic proof is the primary oracle; symbolic checks expose field, boundary, kernel, pressure and action errors. General finite-spectrum API use retains its explicit commensurability and stationarity hypotheses and does not extend the acoustic theorem.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.177.0` with provenance `campaigns/P251-triangular-acoustic-array/adjudication.yaml`.

- `proposals/P251-cosserat-from-vortex-euler/attempts/0165/independent-review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0161/planar-acoustic-array.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0161/verify.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0161/repaired-run.stdout`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0161/receipt.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0146/acoustic-normal-form.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0146/verify.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0146/repaired-run.txt`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0167/receipt.md`
- `src/substrate_framework/euler_acoustic.py`
- `tests/test_euler_acoustic.py`
- `src/substrate_framework/euler_phase.py`
