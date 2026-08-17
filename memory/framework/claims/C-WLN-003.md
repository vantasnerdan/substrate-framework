---
description: Accepted framework claim C-WLN-003
author: framework-registry
created: '2026-08-17T20:45:00Z'
updated: '2026-08-17T20:45:00Z'
tags:
- substrate-framework
- accepted-claim
- C-WLN-003
category: claims
confidence: established
status: active
---
# C-WLN-003

## Statement
The accepted statement is reproduced exactly from the claim registry.

Under an orientation-preserving local reparametrization tau'=f(tau) with smooth positive rate r=df/dtau, dot(x)'=dot(x)/r and e'=e/r make the massive or massless einbein action exactly invariant after the integration Jacobian. On an interval where e is smooth and positive, df/dtau=e/e0 locally attains any declared constant e0>0; global domain and endpoint coverage remain separate requirements. Under a positive pointwise conformal factor Omega, the massless density is invariant for g->Omega^2*g and e->Omega^2*e, while for m>0 the mass term changes by -(Omega^2-1)*e*(m*c0)^2/2 and is nonzero at each point where positive Omega differs from one. These are worldline action identities, not spacetime field-theory conformal symmetry.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: Reparametrizations are smooth, local, orientation preserving, and have strictly positive rate., The einbein and target constant are positive, and local gauge construction does not imply global interval coverage., The conformal factor is positive and the massive obstruction assumes positive m and c0., No Maxwell, Yang-Mills, gravitational, helicity, or spacetime Weyl-gauge conclusion is imported.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.160.0` with provenance `campaigns/P227-einbein-worldline-harvest/adjudication.yaml`.

- `campaigns/P227-einbein-worldline-harvest/verify.py`
- `campaigns/P227-einbein-worldline-harvest/reviews/independent_worldline_lorentz_review.py`
- `campaigns/P227-einbein-worldline-harvest/attempts/0002/result.yaml`
- `memory/vantasner/decisions/C-WLN-003-review.md`
- `src/substrate_framework/relativistic_particle.py`
- `tests/test_relativistic_particle.py`
