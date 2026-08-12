---
description: Review request for proposed claim C-WLN-003
author: codex
created: '2026-08-12T11:36:38+02:00'
updated: '2026-08-12T11:36:38+02:00'
tags:
- substrate-framework
- claim-review-request
category: decisions
confidence: working
status: active
---

## Claim Under Review

The requested statement has three narrowly coupled parts. Under an
orientation-preserving local reparametrization `tau'=f(tau)` with smooth
positive rate `r=df/dtau`, `dot(x)'=dot(x)/r` and `e'=e/r`, the massive or
massless einbein action is exactly invariant after the Jacobian. On an interval
where `e` is smooth and positive, `df/dtau=e/e0` locally attains any declared
constant `e0>0` (a global gauge additionally requires the integrated map to
cover the intended parameter interval). Finally, for positive conformal factor
`Omega`, the massless density is invariant under
`g->Omega^2*g, e->Omega^2*e`; for `m>0` the mass term changes by
`-(Omega^2-1)*e*(m*c0)^2/2`, which vanishes for the identity `Omega=1` and is
nonzero for nontrivial positive `Omega!=1`.

## Sourced Inputs

The review inputs are release `v0.159.0`, P227,
`worldline_reparametrization_residual`, `constant_e_gauge_rate`,
`massless_worldline_weyl_residual`, `massive_mass_term_weyl_change`, their
tests, and the P227 verifier.

## Independence Requested

The reviewer should transform the density and integration measure directly,
integrate the proposed local gauge rate under explicit regularity assumptions,
and compare massive and massless Weyl terms without calling the residual
helpers.

## Verification and Sensitivity

Author evidence is exact SymPy algebra. Wrong einbein reparametrization weight,
wrong Weyl weight, wrong rate, and conflating `Omega=1` with a nontrivial
rescaling must fail. The new verifier explicitly separates the identity
transformation from the massive obstruction.

## Framework Compatibility

The Weyl statement is only a one-dimensional worldline action identity. It is
not Maxwell/Yang-Mills conformal symmetry, a photon count, a graviton claim, or
a spacetime Weyl-gauge theory. The constant-e result is local unless global
domain and endpoint conditions are separately supplied.

## Dependency and Consumer Replay

Replay the canonical worldline tests, P227 verifier, and tutorial consumers.
No accepted field-theory claim may consume this as a conformal-symmetry
derivation.

## Competing Candidate Audit

P227's coordinate route provides the direct transformation oracle for the
selected small residual APIs; document-local duplication remains the rejected
integration baseline.

## Four-Axis Decision Requested

The proposed axes are verification `symbolic_verified`, review `unaudited`,
compatibility `unassessed`, epistemic `proposed`, and relationship
`challenges: []`. The reviewer must assign the actual compatibility verdict.

## Promotion Boundary and Continuation

Promotion is reviewer-owned and must preserve every regularity and nontrivial-
transformation qualifier. Overbroad field-theory or global-gauge language is a
reason to narrow, not accept.

## Cross-References

The canonical goal is issue #59 and the source proposal is P227.
