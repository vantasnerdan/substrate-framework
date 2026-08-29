---
description: Accepted framework claim C-GOT-002
author: framework-registry
created: '2026-08-29T07:06:51+00:00'
updated: '2026-08-29T07:06:51+00:00'
tags:
- substrate-framework
- accepted-claim
- C-GOT-002
category: claims
confidence: established
status: active
---
# C-GOT-002

## Statement
The accepted statement is reproduced exactly from the claim registry.

C-GOT-002 (exact complete optical ADM bijection): Let N>0 be independent, gamma_ij be real symmetric positive definite, V^i be real, c>0, v=V/c, and s in {-1,+1}. The same block formulas define every mostly-plus Lorentzian metric whose chosen spatial block is positive definite and whose temporal Schur complement is negative. The inverse reconstruction is gamma=g_spatial, v=s gamma^-1 g_mixed, N=sqrt(-(g_00-g_mixed^T gamma^-1 g_mixed)), and V=c v. The metric determinant is -N^2 det(gamma), sqrt(-g)=N sqrt(det gamma), and the ten-component point-map Jacobian in coordinates (N,v^1,v^2,v^3,gamma_11,gamma_22,gamma_33,gamma_12,gamma_13,gamma_23) has determinant -2*s*N*det(gamma), hence is nonsingular on the stated domain. Eight fixed-seed nonsingular SPD numerical round trips agree at infinity norm below 2e-13 as regression evidence.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: The chart is restricted to a chosen spacelike foliation with positive-definite spatial block and negative temporal Schur complement., Coordinates use x0=c*t; N and c are positive and the flow orientation is exactly plus or minus one.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.170.0` with provenance `campaigns/P248-sciama-gothic-metric-audit/adjudication.yaml`.

- `campaigns/P248-sciama-gothic-metric-audit/attempts/0003/result.yaml`
- `campaigns/P248-sciama-gothic-metric-audit/attempts/0003/numerical-output.txt`
- `campaigns/P248-sciama-gothic-metric-audit/companion/exact_metric_checks.py`
- `campaigns/P248-sciama-gothic-metric-audit/companion/numerical_checks.py`
- `campaigns/P248-sciama-gothic-metric-audit/reviews/C-GOT-001-C-GOT-006-review.md`
- `src/substrate_framework/optical_gothic.py`
- `tests/test_optical_gothic.py`
