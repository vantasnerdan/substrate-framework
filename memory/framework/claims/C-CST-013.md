---
description: Accepted framework claim C-CST-013
author: framework-registry
created: '2026-09-05T10:21:58+00:00'
updated: '2026-09-05T10:21:58+00:00'
tags:
- substrate-framework
- accepted-claim
- C-CST-013
category: claims
confidence: established
status: active
---
# C-CST-013

## Statement
The accepted statement is reproduced exactly from the claim registry.

Let u be a smooth stationary mean-zero periodic incompressible Euler
field, p its pressure divided by density, A=u.grad, K=k*kappa with
|kappa|=1, and D,V perpendicular to kappa. Prepare the actual Kelvin
displacement D and independently specified common velocity V. Their
zero-order material displacement is U0=D+tV. Write
eta=U0+i*k*chi+O(k²), where chi is mean-zero and solenoidal.
With microscopic mean-preserving Leray P, its actual first cell obeys

  chi_tt+2P A chi_t+P(A²+Hess p)chi
    =-2P[(kappa.u)V]+P[(kappa.grad p)U0+kappa(U0.grad p)],
  chi(0)=0, chi_t(0)=-P[(kappa.u)D+kappa(u.D)].

For the physical Eulerian velocity mean m, set X=D+integral m dt.
Put a=kappa.u and P_kappa=I-kappa tensor kappa. Its complete second
fixed-time directional spatial jet is

  m_t=k² P_kappa{<a²>U0+<a chi_t+u(kappa.chi_t)>
                   +<(kappa.grad p)chi+grad p(kappa.chi)>}+O(k³).

The corresponding mean material rate differs by
m-<eta>_t=-k²<a chi-u(kappa.chi)>+O(k³). These exact observed rows
require no inverse for the unknown second cell: its divergence
constraint and periodic integration by parts remove it. The identity
P_K grad f=-i*k*P_K kappa f retains the pressure numerator, including
the slow harmonic, rather than differentiating a singular bare
projector estimate.

In0153's bounded constant-curl insertion family u_N=u_wave+h_N,
normalized H2 defect epsilon_N tends to zero while the global C3
bound remains fixed and prescribed C-CST-011 local optical/tube
margins are retained. On each fixed time window, the actual first
cell Euler/Lin H1 estimate controls the physical observed second
coefficients by C_T*epsilon_N, uniformly in the quadrature period
and direction. Low-frequency insertion self-beats are included.

Averaging the actual coherently prepared phase actions over one
whole-field O(3), phase and time-reversal law before elimination
gives initial physical symplectic mass rho. For the actual averaged
transverse solution rows f_N,g_N and W_N=f_N*g_N,t-g_N*f_N,t,

  M_N=rho/W_N,
  K_N=rho*(f_N,t*g_N,tt-g_N,t*f_N,tt)/W_N².

A selected finite cell and then sufficiently small nonzero k have
M_N>0 and K_N/(rho*k²)=2*v_wave²/15+O_T(epsilon_N)+O_N,T(k)>0,
on the fixed window, while retaining the same local optical object.
The O(epsilon_N*k²) time-dependent mass, full moving connection and
physical/canonical current difference remain explicit in this action.

This is an actual prepared mean/action transfer theorem, not full
Bloch-operator C2 convergence, a period-independent analytic radius,
acoustic-time homogenization, or an autonomous coupled Cosserat law.
Optical density is positive at the selected finite cell; no positive
lower bound in an infinite-volume limit is asserted.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: C-CST-011. Assumptions: Constant-density incompressible Euler with complete periodic pressure; microscopic Leray preserves constants and the slow directional harmonic is P_kappa. Kelvin D and common V have independently specified circulation data., The specific constant-curl Fourier/Herglotz insertion, its fixed total variation, normalized H2 convergence and ordered finite-cell optical transfer are proved in0153/0170; no generic small L2 perturbation theorem is substituted. C-CST-011 supplies the retained local optical object at its accepted scope., Standard periodic Fourier analysis, smooth finite-time Euler/Lin well-posedness, normalized energy estimates and differentiable finite-cell Bloch dependence are declared mathematical imports. The response-specific cancellation and period-uniform bound are derived here., Select the optical target and fixed time window, then a sufficiently accurate finite cell, then the nonzero macro wavevector. Full time connections and physical observation rows remain; stationarity is not a Markov closure premise., Analytic PDE proof is the primary oracle; exact symbolic pressure/current checks and independent Cartesian API tests are exposing corroboration, not numerical continuum-existence evidence.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.178.0` with provenance `campaigns/P251-prepared-mean-transfer/adjudication.yaml`.

- `proposals/P251-cosserat-from-vortex-euler/attempts/0170/prepared-mean-jet.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0170/verify.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0170/receipt.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0173/independent-review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0153/same-field.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0151/one-wave-response.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0154/coherent-physical-action.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0177/README.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0177/first-pytest.stdout`
- `src/substrate_framework/euler_acoustic.py`
- `tests/test_euler_acoustic_response.py`
- `src/substrate_framework/euler_phase.py`
