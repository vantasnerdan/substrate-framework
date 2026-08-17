---
description: Review request for proposed claim C-WLN-001
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

The requested decision concerns this exact classical worldline statement. For
a supplied nondegenerate pseudo-Riemannian metric in mostly-plus convention, a
timelike velocity norm `sigma<0`, and exact positive `e`, `m`, and `c0`, the
density `L=sigma/(2e)-e*(m*c0)^2/2` has auxiliary equation
`sigma+e^2*(m*c0)^2=0`. Its unique positive solution is
`e=sqrt(-sigma)/(m*c0)`, and eliminating it gives exactly
`-m*c0*sqrt(-sigma)`. With canonical covector momentum `p_mu`, the inverse
Legendre map is `dot(x)^mu=e*g^(mu nu)*p_nu` and the Hamiltonian is the pure
constraint `e*(g^(mu nu)*p_mu*p_nu+(m*c0)^2)/2`.

## Sourced Inputs

The author packet points the reviewer to release `v0.159.0`, proposal P227,
`relativistic_particle.py`, `tests/test_relativistic_particle.py`, and
`verify.py`. The merged tutorials in PRs #38 and #39 are provenance, not
authority. No accepted scientific claim is proposed as a premise.

## Independence Requested

The reviewer should rederive the `e` Euler equation, positive branch,
substitution, covector momentum, inverse Legendre map, and Hamiltonian without
calling `massive_einbein_ledger` or `einbein_hamiltonian_ledger`. The branch and
index placement are load-bearing.

## Verification and Sensitivity

The author evidence is exact SymPy algebra. Existing mutations reverse the
mass-term sign and replace inverse-metric covector raising with multiplication
by the covariant metric; both must fail. The reviewer should inspect the actual
expressions for unresolved roots or assumptions and test at least one exact
non-diagonal metric.

## Framework Compatibility

This is a local classical auxiliary-field identity. It does not select a
background metric, quantize the particle, identify a species, or create a
field-theory claim. Positive mass is essential and the zero-mass branch belongs
to C-WLN-002.

## Dependency and Consumer Replay

The direct implementation and tutorial consumers must pass together. Public
exports and the P227 verifier are the affected framework surfaces; no accepted
claim consumer currently depends on this proposed result.

## Competing Candidate Audit

P227 preregistered a scalar-contraction ledger, a larger coordinate-function
derivation, and document-local duplication. The ledger was selected for exact
branch typing and reuse; the coordinate route remains the requested independent
derivation.

## Four-Axis Decision Requested

The author requested an independent decision, not acceptance by this packet.
The proposed axes are verification `symbolic_verified`, review `unaudited`,
compatibility `unassessed`, epistemic `proposed`, with no `supersedes`
relationship.

## Promotion Boundary and Continuation

Acceptance requires a separate reviewer-owned governance transaction with
registry, release, generated docs, accepted memory, and full closure checks.
If any branch or index statement fails, retain the identifier and narrow or
repair the proposal; do not treat artifact merge as promotion.

## Cross-References

The canonical goal is issue #59. Evidence lives in P227 and the named canonical
module and tests.
