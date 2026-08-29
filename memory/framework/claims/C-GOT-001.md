---
description: Accepted framework claim C-GOT-001
author: framework-registry
created: '2026-08-29T07:06:51+00:00'
updated: '2026-08-29T07:06:51+00:00'
tags:
- substrate-framework
- accepted-claim
- C-GOT-001
category: claims
confidence: established
status: active
---
# C-GOT-001

## Statement
The accepted statement is reproduced exactly from the claim registry.

C-GOT-001 (exact constrained optical-gothic map): Let gamma_ij be a real symmetric positive-definite 3-by-3 optical tensor, V^i a real flow, c>0, v=V/c, s in {-1,+1}, nbar=(det gamma)^(1/3), and N=nbar^(-1/2). The mostly-plus block metric g_00=-N^2+gamma_ij v^i v^j, g_0i=s gamma_ij v^j, g_ij=gamma_ij is Lorentzian, has det(g)=-nbar^2 and sqrt(-g)=nbar, and has inverse g^00=-1/N^2, g^0i=s v^i/N^2, g^ij=gamma^-1^ij-v^i v^j/N^2. Therefore its gothic inverse density has components gothic^00=-nbar^2, gothic^0i=s nbar^2 v^i, and gothic^ij=nbar gamma^-1^ij-nbar^2 v^i v^j, reproducing the issue-184 paper for s=+1. The construction is an exact nine-field chart whose image obeys N^2*nbar=1; it is not a bijection to the ten independent components of a general foliation-adapted Lorentzian metric, with diag(-4,1,1,1) at gamma=I and V=0 an explicit excluded witness.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: Coordinates use x0=c*t; c and N are positive, gamma is real symmetric positive definite, V is real, and s is exactly plus or minus one., The positive real cube root defines nbar and fixes the determinant-slaved lapse.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.170.0` with provenance `campaigns/P248-sciama-gothic-metric-audit/adjudication.yaml`.

- `campaigns/P248-sciama-gothic-metric-audit/attempts/0002/result.yaml`
- `campaigns/P248-sciama-gothic-metric-audit/attempts/0002/sympy-output.txt`
- `campaigns/P248-sciama-gothic-metric-audit/companion/exact_metric_checks.py`
- `campaigns/P248-sciama-gothic-metric-audit/reviews/C-GOT-001-C-GOT-006-review.md`
- `src/substrate_framework/optical_gothic.py`
- `tests/test_optical_gothic.py`
