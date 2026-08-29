---
description: Accepted framework claim C-GOT-006
author: framework-registry
created: '2026-08-29T07:06:51+00:00'
updated: '2026-08-29T07:06:51+00:00'
tags:
- substrate-framework
- accepted-claim
- C-GOT-006
category: claims
confidence: established
status: active
---
# C-GOT-006

## Statement
The accepted statement is reproduced exactly from the claim registry.

C-GOT-006 (symbolic-verified synthesized core optical action theorem): On the C-GOT-002 domain, pull any differentiable covariant metric action S[g,fields] back through the complete pointwise optical map q=(N,v,gamma) -> g(q), where v=V/c. Its Euler derivatives obey E_q=J^T E_g; because det(J)=-2*s*N*det(gamma) is nonzero, E_q=0 if and only if E_g=0. Specializing S to the accepted mostly-plus Einstein-scalar action C-STG-001 therefore gives exactly its Einstein and scalar equations in the complete optical variables, without losing a metric equation. Supplying the same reconstructed metric to accepted C-WLN-001 and C-WLN-002 simultaneously gives the massive square-root/ Hamiltonian sector and the massless null/affine-geodesic sector at their accepted local scopes. This is a complete conditional infrared optical presentation, not a derivation of a substrate: it derives no microscopic constitutive action, massless spin-two emergence, universal material coupling, preferred foliation, Newton constant, impedance law, ether ontology, or empirical equivalence.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Theorem Classification
Category is `synthesized`; layer is `core`.
Structural gap: Closes the missing tenth-field and variational bridges between a complete optical metric chart, accepted covariant field dynamics, and accepted massive and massless supplied-metric worldline sectors.

Glue proof: `sympy` at `campaigns/P248-sciama-gothic-metric-audit/companion/exact_metric_checks.py` entrypoint `run`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: C-STG-001, C-WLN-001, C-WLN-002, C-GOT-002. Assumptions: The complete optical point map is used on C-GOT-002's exact domain with q=(N,v,gamma), compact metric variations or the corresponding gravitational boundary term, and all accepted dependency hypotheses retained., The same reconstructed supplied metric is used independently in the accepted massive and massless worldline identities; no field-to-particle realization is inferred.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.170.0` with provenance `campaigns/P248-sciama-gothic-metric-audit/adjudication.yaml`.

- `campaigns/P248-sciama-gothic-metric-audit/attempts/0003/result.yaml`
- `campaigns/P248-sciama-gothic-metric-audit/attempts/0003/lean-correction-output.txt`
- `campaigns/P248-sciama-gothic-metric-audit/companion/exact_metric_checks.py`
- `campaigns/P248-sciama-gothic-metric-audit/reviews/C-GOT-001-C-GOT-006-review.md`
- `formal/SubstrateFramework/OpticalGothic.lean`
- `src/substrate_framework/optical_gothic.py`
- `tests/test_optical_gothic.py`
