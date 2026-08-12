---
description: Review request for proposed claim C-WLN-002
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

The requested decision concerns the massless specialization made before
auxiliary-field elimination. For exact positive `e`, setting `m=0` in the
einbein density gives the regular expression `L=sigma/(2e)`; its `e` equation
is exactly the null constraint `sigma=0`. The coordinate Euler equation in a
supplied nondegenerate metric is
`ddot(x)^a=-Gamma^a_bc*dot(x)^b*dot(x)^c+(dot(e)/e)*dot(x)^a`; constant positive
`e` therefore gives the affine null-geodesic equation.

## Sourced Inputs

The reviewer should read release `v0.159.0`, P227, `massless_einbein_ledger`,
`einbein_geodesic_acceleration`, their exact tests, and the P227 verifier.
Merged tutorial prose is source provenance only.

## Independence Requested

The reviewer should vary the massless action directly with respect to `e` and
the coordinates, reconstruct Christoffel symbols from metric derivatives, and
confirm the lower-index Euler residual. The massive positive root must not be
used and then specialized to `m=0`.

## Verification and Sensitivity

Author evidence is exact SymPy algebra in dimensions three and four. Tests
reject use of the massive eliminator at zero mass and a wrong geodesic-force
sign. The review should include a nonconstant metric derivative tensor and
inspect all residual components.

## Framework Compatibility

This claim is classical, local, and conditional on supplied exact metric data.
It derives neither a field equation nor a massless particle representation.
Nullness follows as a constraint and is not assumed through a singular
square-root action.

## Dependency and Consumer Replay

The canonical worldline and pseudo-Riemannian tests, P227 verifier, and both
tutorial suites form the direct replay. No accepted claim currently imports
this proposal.

## Competing Candidate Audit

P227 registered scalar-contraction, coordinate-function, and document-local
routes. The coordinate Euler calculation is retained as the independent oracle
for the smaller reusable ledger.

## Four-Axis Decision Requested

The proposed axes before independent review are verification
`symbolic_verified`, review `unaudited`, compatibility `unassessed`, epistemic
`proposed`, with no replacement relationship.

## Promotion Boundary and Continuation

Only a reviewer may accept and promote this claim in a complete governance
transaction. Failure of the coordinate Euler equation requires repair or
narrowing while issue #59 remains open.

## Cross-References

The canonical goal is issue #59; evidence is under P227 and
`src/substrate_framework/relativistic_particle.py`.
