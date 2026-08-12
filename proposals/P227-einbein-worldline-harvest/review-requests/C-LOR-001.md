---
description: Review request for proposed claim C-LOR-001
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

The requested claim is only the exact induced geometry of the future unit
timelike-vector orbit in mostly-plus flat spacetime. In 2+1 dimensions the
standard hyperboloid pullback is `diag(1,sinh(eta)^2)` on `H^2`; in 3+1 it is
`diag(1,sinh(eta)^2,sinh(eta)^2*sin(theta)^2)` on `H^3`.

## Sourced Inputs

Review release `v0.159.0`, P227, `lorentz_orbits.py`, its tests, and the P227
verifier. The earlier H3 helper is unpromoted provenance, not authority.

## Independence Requested

Differentiate the displayed hyperboloid embeddings and pull back the ambient
metrics without calling `unit_timelike_vector_orbit_metric`. Independently
verify the unit-hyperboloid constraint and induced signature.

## Verification and Sensitivity

Author evidence is exact SymPy matrix algebra. A sign-flipped ambient
metric, exchanged angular embedding, or unsupported dimension must not satisfy
the expected matrix.

## Framework Compatibility

This describes a single future unit timelike vector. It is not the orbit of a
timelike plane, tube worldsheet, Lorentz-invariant probability measure, or
physical ensemble.

## Dependency and Consumer Replay

Replay `tests/test_lorentz_orbits.py`, P227, and the flux-tube callers that use
the older H3 surface. The exact H2/H3 result has no accepted-claim premise.

## Competing Candidate Audit

P227 selected a dimension-aware pullback over duplicated tutorial formulas;
the direct embedding calculation remains an independent rederivation.

## Four-Axis Decision Requested

The proposed axes are verification `symbolic_verified`, review `unaudited`,
compatibility `unassessed`, epistemic `proposed`, with no supersession.

## Promotion Boundary and Continuation

A distinct reviewer must keep the claim at induced-metric scope. Any extension
beyond dimensions three and four requires new evidence.

## Cross-References

The canonical goal is issue #59; the proposal and verifier are under P227.
