---
description: Accepted framework claim C-WLN-001
author: framework-registry
created: '2026-08-17T20:45:00Z'
updated: '2026-08-17T20:45:00Z'
tags:
- substrate-framework
- accepted-claim
- C-WLN-001
category: claims
confidence: established
status: active
---
# C-WLN-001

## Statement
The accepted statement is reproduced exactly from the claim registry.

For a supplied nondegenerate pseudo-Riemannian metric in mostly-plus convention, timelike velocity norm sigma<0, and positive exact einbein e, mass m, and signal speed c0, the worldline density L=sigma/(2e)-e*(m*c0)^2/2 has auxiliary equation sigma+e^2*(m*c0)^2=0 with unique positive solution e=sqrt(-sigma)/(m*c0). Eliminating e gives exactly -m*c0*sqrt(-sigma). With canonical covector momentum p_mu, the inverse Legendre map is dot(x)^mu=e*g^(mu nu)*p_nu and the Hamiltonian is the pure constraint e*(g^(mu nu)*p_mu*p_nu+(m*c0)^2)/2. This is a conditional local classical identity; it selects no background, particle species, quantization, or field theory.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: The metric is supplied, exact, nondegenerate, and uses the displayed mostly-plus convention., The velocity norm is provably negative and e, m, and c0 are provably positive., Momentum is a covector and index raising uses the inverse supplied metric., The statement is local classical worldline mechanics and includes no background-selection or quantization premise.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.160.0` with provenance `campaigns/P227-einbein-worldline-harvest/adjudication.yaml`.

- `campaigns/P227-einbein-worldline-harvest/verify.py`
- `campaigns/P227-einbein-worldline-harvest/reviews/independent_worldline_lorentz_review.py`
- `campaigns/P227-einbein-worldline-harvest/attempts/0002/result.yaml`
- `memory/vantasner/decisions/C-WLN-001-review.md`
- `src/substrate_framework/relativistic_particle.py`
- `tests/test_relativistic_particle.py`
