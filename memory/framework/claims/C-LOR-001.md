---
description: Accepted framework claim C-LOR-001
author: framework-registry
created: '2026-08-17T20:45:00Z'
updated: '2026-08-17T20:45:00Z'
tags:
- substrate-framework
- accepted-claim
- C-LOR-001
category: claims
confidence: established
status: active
---
# C-LOR-001

## Statement
The accepted statement is reproduced exactly from the claim registry.

In mostly-plus flat spacetime, the future-directed unit timelike-vector orbit has exact induced metric diag(1,sinh(eta)^2) on H^2 in 2+1 dimensions and diag(1,sinh(eta)^2,sinh(eta)^2*sin(theta)^2) on H^3 in 3+1 dimensions, obtained by pulling back the ambient metric along the standard hyperboloid embeddings. This claim concerns one unit timelike vector only; it is not a timelike-plane or tube-worldsheet orbit, an invariant probability measure, or a physical ensemble.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: The ambient metrics are diag(-1,1,1) and diag(-1,1,1,1) in the displayed mostly-plus convention., The embeddings use the standard future sheet and the displayed hyperbolic polar coordinates., Only spacetime dimensions three and four and only the induced metric are claimed.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.160.0` with provenance `campaigns/P227-einbein-worldline-harvest/adjudication.yaml`.

- `campaigns/P227-einbein-worldline-harvest/verify.py`
- `campaigns/P227-einbein-worldline-harvest/reviews/independent_worldline_lorentz_review.py`
- `campaigns/P227-einbein-worldline-harvest/attempts/0002/result.yaml`
- `memory/vantasner/decisions/C-LOR-001-review.md`
- `src/substrate_framework/lorentz_orbits.py`
- `tests/test_lorentz_orbits.py`
