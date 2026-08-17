---
description: Independent review and acceptance of C-LOR-001 (future unit timelike-vector orbit metrics, exact H^2 and H^3 pullbacks)
author: vantasner-review
created: '2026-08-17T20:50:00Z'
updated: '2026-08-17T20:50:00Z'
tags:
- substrate-framework
- claim-review
- P227
- lorentz-orbit
- hyperboloid
category: decisions
confidence: established
status: active
---
# Review of C-LOR-001

## Claim Under Review

The conditional exact statement is only the exact induced geometry of
the future unit timelike-vector orbit in mostly-plus flat spacetime.
In 2+1 dimensions the standard hyperboloid pullback is
`diag(1, sinh(eta)^2)` on `H^2`; in 3+1 it is
`diag(1, sinh(eta)^2, sinh(eta)^2 sin(theta)^2)` on `H^3`.

## Sourced Inputs

The review read release `v0.159.0`, proposal P227, the canonical
`unit_timelike_vector_orbit_metric` API in
`src/substrate_framework/lorentz_orbits.py`, the canonical tests,
and the P227 verifier.

## Independence

The review differentiates the displayed hyperboloid embeddings and
pulls back the ambient metrics by hand (no import of
`unit_timelike_vector_orbit_metric`). The unit-hyperboloid constraint
and induced signature are verified independently.

### Independent derivation

1. **2+1D embedding** (`H^2`):
   `X(t, eta, theta) = (cosh eta, sinh eta cos theta, sinh eta sin theta)`.
   - Unit-timelike norm (mostly-plus `eta = diag(-1, 1, 1)`):
     `X^T eta X = -cosh^2(eta) + sinh^2(eta) = -1` (verified).
   - Tangent vectors:
     `dX/deta = (sinh eta, cosh eta cos theta, cosh eta sin theta)`,
     `dX/dtheta = (0, -sinh eta sin theta, sinh eta cos theta)`.
   - Pullback `g_ab = dX/dq^a eta dX/dq^b`:
     `g_{eta,eta} = -sinh^2 eta + cosh^2 eta = 1`,
     `g_{eta,theta} = 0`,
     `g_{theta,theta} = sinh^2 eta (sin^2 theta + cos^2 theta) = sinh^2 eta`.
   - Result: `g_H2 = diag(1, sinh^2 eta)`, matching the claim.

2. **3+1D embedding** (`H^3`):
   `X(t, eta, theta, phi) = (cosh eta, sinh eta sin theta cos phi, sinh eta sin theta sin phi, sinh eta cos theta)`.
   - Unit-timelike norm:
     `X^T eta X = -cosh^2 eta + sinh^2 eta (sin^2 theta + cos^2 theta) = -1` (verified).
   - Pullback along `(eta, theta, phi)`:
     `g_{eta,eta} = 1`,
     `g_{theta,theta} = sinh^2 eta`,
     `g_{phi,phi} = sinh^2 eta sin^2 theta`,
     cross terms zero (verified).
   - Result: `g_H3 = diag(1, sinh^2 eta, sinh^2 eta sin^2 theta)`,
     matching the claim.

3. **Comparison with canonical helper**: the independent pullback is
   subtracted from `unit_timelike_vector_orbit_metric(3)` and
   `unit_timelike_vector_orbit_metric(4)`; both differences are
   identically zero matrices (verified by `sp.simplify`).

## Verification Status

Exact symbolic substitution proves the unit-timelike norm (`-1`), the
3+3 (2+1D) and 4+4 (3+1D) pullback entries, and the off-diagonal
zeros. The independent rederivation equals the canonical helper
without using it. The claim earns `symbolic_verified`.

## Sensitivity and Counterexamples

- **Sign-flipped ambient metric** (`eta = diag(+1, -1, -1)` instead
  of `diag(-1, +1, +1)`): the H^2 pullback becomes
  `diag(-1, -sinh^2 eta)`, which does NOT match the canonical
  `diag(1, sinh^2 eta)`. Mutation rejected.
- **Wrong dimension** (`unit_timelike_vector_orbit_metric(5)`):
  raises `ValueError("spacetime_dimension must be 3 or 4")`; no
  metric is defined for unsupported dimensions. Mutation rejected.
- **Swapped angular embedding**: an embedding that exchanges `theta`
  and `phi` for the spatial direction preserves the unit-timelike
  norm `-1` (because the angular factors are still a unit vector),
  but yields a different pullback metric; the canonical helper is
  not recovered. Mutation rejected on exact pullback.

## Framework Compatibility

The claim describes a single future unit timelike vector orbit. It
is NOT the orbit of a timelike plane, a tube worldsheet, a
Lorentz-invariant probability measure, or a physical ensemble. No
field, action, or measure is identified merely by dimension.

## Dependency and Consumer Replay

The direct dependency is `sympy.diff`, `sympy.trigsimp`, and
`sympy.Matrix` (all stdlib). Consumers are
`unit_timelike_vector_orbit_metric`, the canonical
`tests/test_lorentz_orbits.py` (3 tests pass), and the P227 verifier
(all 14 checks pass).

## Competing Candidate Audit

P227's candidate A selected a dimension-aware pullback (the
canonical helper) over duplicated tutorial formulas; the direct
embedding calculation above is the candidate-B independent
rederivation, validating the chosen helper without importing it.

## Four-Axis Decision

The four governance axes are assigned as follows. Verification is symbolic_verified, review is accepted,
`compatibility=compatible_extension`, `epistemic=active`. No
supersession. Scope ceiling is claim content: induced metric on a
single future unit timelike-vector orbit, dimensions 3 and 4 only.
Any extension beyond these dimensions requires new evidence.

## Cross-References

Source proposal: `proposals/P227-einbein-worldline-harvest/proposal.yaml`.
Source PR: PR #60 (merged). Review request:
`proposals/P227-einbein-worldline-harvest/review-requests/C-LOR-001.md`.
