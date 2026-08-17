---
description: Independent review and acceptance of C-WLN-003 (worldline reparametrization, local constant-e gauge, massless Weyl identity, massive obstruction)
author: vantasner-review and codex
created: '2026-08-17T20:30:58Z'
updated: '2026-08-17T20:45:00Z'
tags:
- substrate-framework
- claim-review
- P227
- worldline
- reparametrization
- weyl
category: decisions
confidence: established
status: active
---
# Review of C-WLN-003

## Claim Under Review

The conditional exact statement has three narrowly coupled parts.
Under an orientation-preserving local reparametrization
`tau'=f(tau)` with smooth positive rate `r=df/dtau`, the
transformation laws are `dot(x)'=dot(x)/r` and `e'=e/r`; the massive
or massless einbein action is exactly invariant after the Jacobian.
On an interval where `e` is smooth and positive, `df/dtau=e/e0`
locally attains any declared constant `e0>0` (a global gauge
additionally requires the integrated map to cover the intended
parameter interval). For positive conformal factor `Omega`, the
massless density is invariant under
`g -> Omega^2 g, e -> Omega^2 e`; for `m > 0` the mass term changes
by `-(Omega^2-1) e (m c0)^2 / 2`, which vanishes for `Omega=1` and
is nonzero for nontrivial positive `Omega != 1`.

## Sourced Inputs

The review read release `v0.159.0`, campaign P227, the canonical
`worldline_reparametrization_residual`,
`constant_e_gauge_rate`, `massless_worldline_weyl_residual`, and
`massive_mass_term_weyl_change` APIs in
`src/substrate_framework/relativistic_particle.py`, the canonical
tests, and the P227 verifier.

## Independence

The review transforms the density and integration measure directly
(without using the residual helpers), integrates the proposed local
gauge rate under explicit regularity assumptions, and compares the
massive and massless Weyl terms using only `sp.simplify` algebra.

### Independent derivation

1. **Reparametrization**: under `tau -> tau'(tau)` with
   `r = df/dtau > 0`, the integration measure transforms as
   `dtau' = r dtau`, and the velocity/einbein transforms as
   `dot x' = dot x / r`, `e' = e / r`. The transformed density is
   `L' = sigma'/(2 e') - e' (m c0)^2 / 2` with
   `sigma' = sigma / r^2` (because
   `sigma' = g_ij (dot x^i / r)(dot x^j / r) = sigma/r^2`). The
   transformed integrand is
   `r L' = sigma/(2 e) - e (m c0)^2 / 2 = L`. The residual is
   `r L' - L = 0` (verified). Both massive and massless densities
   are exactly invariant.

2. **Local constant-e gauge**: solving the ODE
   `df/dtau = e(tau) / e_target` gives `f(tau) = (1/e_target) ∫ e dtau + f(0)`;
   then `e / (df/dtau) = e_target`, a constant. This is a local
   statement: the ODE requires `e > 0` on the interval (assumed);
   global domain/endpoint coverage requires separate analysis.

3. **Massless Weyl identity**: under
   `g -> Omega^2 g, e -> Omega^2 e`,
   `sigma' = Omega^2 sigma`,
   `L' = (Omega^2 sigma) / (2 Omega^2 e) = sigma/(2 e) = L`. The
   residual is `L' - L = 0` (verified).

4. **Massive Weyl obstruction**: the mass term
   `T_massive = -e (m c0)^2 / 2` becomes
   `T_massive' = -(Omega^2 e) (m c0)^2 / 2`, so the change is
   `T_massive' - T_massive = -(Omega^2 - 1) e (m c0)^2 / 2`
   (verified). For `Omega = 1` this is 0 (identity); for
   `Omega = 2` it is `-3 e (m c0)^2 / 2`, nonzero as claimed.

## Verification Status

Exact symbolic substitution proves all four statements:
reparametrization invariance (both branches), local constant-e gauge
formula, massless Weyl identity, and massive Weyl obstruction.
The identity `Omega=1` and nontrivial `Omega=2` are both confirmed.
The claim earns `symbolic_verified`.

## Sensitivity and Counterexamples

- **Wrong reparametrization weight** (`e' = e/r^2` instead of
  `e/r`): the residual `r L' - L` becomes
  `c0^2 e m^2 / 2 (1 - 1/r) + (r - 1) sigma / (2e)`, which is
  nonzero. Mutation rejected.
- **Wrong constant-gauge rate** (`df/dtau=e/(2e_target)`): the
  transformed einbein is `2e_target`, not the declared target. Mutation
  rejected.
- **Wrong Weyl weight on e** (`e -> Omega e` instead of
  `Omega^2 e`): the residual
  `L_massive' - L` becomes
  `(Omega sigma + c0^2 e^2 m^2 (1 - Omega) - sigma) / (2 e)`,
  nonzero. Mutation rejected.
- **Conflating `Omega=1` with nontrivial `Omega`**: the residual is 0
  for `Omega=1` (identity, expected) and nonzero for `Omega=2`
  (mutation expected). The verifier explicitly separates the two
  cases.

## Framework Compatibility

The claim is a one-dimensional worldline action identity only. It is
not Maxwell/Yang-Mills conformal symmetry, a photon count, a graviton
claim, or a spacetime Weyl-gauge theory. The constant-e result is
local unless global domain and endpoint conditions are separately
supplied; the massive-mass-term change is a one-dimensional
algebraic identity, not a field-theory conformal claim.

## Dependency and Consumer Replay

The direct dependencies are `sympy.simplify` and the repository's
`exact_real` / `positive_exact` helpers (used by the canonical
APIs). Consumers are the four named helpers, the canonical worldline
tests, the P227 verifier (all 14 checks pass), and the merged
tutorial suites. No accepted field-theory claim may consume this as a
conformal-symmetry derivation.

## Competing Candidate Audit

P227's candidate B (coordinate-function transformation oracle)
provides the direct verification of the candidate-A residual APIs.
The review above is a candidate-B-style rederivation, validating the
chosen ledger without using its named helpers. Document-local
duplication remains the rejected integration baseline.

## Four-Axis Decision

The four governance axes are assigned as follows. Verification is symbolic_verified, review is accepted,
`compatibility=compatible_extension`, `epistemic=active`. No
supersession. Scope ceiling is claim content: one-dimensional
worldline action identity, not a field-theory conformal statement.

## Cross-References

Source campaign: `campaigns/P227-einbein-worldline-harvest/proposal.yaml`.
Source PR: PR #60 (merged). Review request:
`campaigns/P227-einbein-worldline-harvest/review-requests/C-WLN-003.md`.
Executable corrective review:
`campaigns/P227-einbein-worldline-harvest/reviews/independent_worldline_lorentz_review.py`.
