---
description: Accepted framework claim C-GOT-003
author: framework-registry
created: '2026-08-29T07:06:51+00:00'
updated: '2026-08-29T07:06:51+00:00'
tags:
- substrate-framework
- accepted-claim
- C-GOT-003
category: claims
confidence: established
status: active
---
# C-GOT-003

## Statement
The accepted statement is reproduced exactly from the claim registry.

C-GOT-003 (exact material-harmonic compatibility): For positive nbar, material residual M=partial_t nbar+V dot grad(nbar)+nbar div(V), and the exact gothic time-gauge residual derived from C-GOT-001, H_s=partial_t(nbar^2)-s div(nbar^2 V), one has H_s=2*nbar*M-2*nbar*(1+s)*V dot grad(nbar) -nbar^2*(2+s)*div(V). For the material-flow sign s=-1 this reduces to H_-=2*nbar*M-nbar^2 div(V), so M=H_-=0 is equivalent to material continuity plus incompressibility; it does not require spatial optical homogeneity, and V=0 with arbitrary static nbar is an exact counterexample to that inference. For the paper's printed s=+1 sign, imposing M=0 gives H_+=-nbar*(4 V dot grad(nbar)+3 nbar div(V)), exposing the paper's inconsistent sign change between its metric and continuity discussion.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: C-GOT-001. Assumptions: The determinant mean is positive and differentiable, V is differentiable, x0=c*t is fixed, and s is exactly plus or minus one.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.170.0` with provenance `campaigns/P248-sciama-gothic-metric-audit/adjudication.yaml`.

- `campaigns/P248-sciama-gothic-metric-audit/attempts/0002/result.yaml`
- `campaigns/P248-sciama-gothic-metric-audit/companion/exact_compatibility_checks.py`
- `campaigns/P248-sciama-gothic-metric-audit/reviews/C-GOT-001-C-GOT-006-review.md`
- `formal/SubstrateFramework/OpticalGothic.lean`
- `src/substrate_framework/optical_gothic.py`
- `tests/test_optical_gothic.py`
