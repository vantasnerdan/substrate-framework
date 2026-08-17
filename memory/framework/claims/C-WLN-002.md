---
description: Accepted framework claim C-WLN-002
author: framework-registry
created: '2026-08-17T20:45:00Z'
updated: '2026-08-17T20:45:00Z'
tags:
- substrate-framework
- accepted-claim
- C-WLN-002
category: claims
confidence: established
status: active
---
# C-WLN-002

## Statement
The accepted statement is reproduced exactly from the claim registry.

Setting m=0 before auxiliary-field elimination gives the regular worldline density L=sigma/(2e), whose einbein equation is exactly the null constraint sigma=0. For a supplied nondegenerate metric, the coordinate Euler equation is ddot(x)^a=-Gamma^a_bc*dot(x)^b*dot(x)^c+(dot(e)/e)*dot(x)^a. Constant positive e therefore gives the affine null-geodesic equation. This is a local conditional classical identity and derives no field equation, particle representation, or singular massless square-root action.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: The massless specialization is made before solving the massive auxiliary equation., The metric is supplied, exact, nondegenerate, differentiable to the required local order, and uses the displayed convention., The einbein is smooth and positive; its derivative and the coordinate velocity are local exact data., Affine parametrization is concluded only where the einbein is constant and the null constraint is also imposed.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.160.0` with provenance `campaigns/P227-einbein-worldline-harvest/adjudication.yaml`.

- `campaigns/P227-einbein-worldline-harvest/verify.py`
- `campaigns/P227-einbein-worldline-harvest/reviews/independent_worldline_lorentz_review.py`
- `campaigns/P227-einbein-worldline-harvest/attempts/0002/result.yaml`
- `memory/vantasner/decisions/C-WLN-002-review.md`
- `src/substrate_framework/relativistic_particle.py`
- `tests/test_relativistic_particle.py`
