---
description: Accepted framework claim C-GOT-005
author: framework-registry
created: '2026-08-29T07:06:51+00:00'
updated: '2026-08-29T07:06:51+00:00'
tags:
- substrate-framework
- accepted-claim
- C-GOT-005
category: claims
confidence: established
status: active
---
# C-GOT-005

## Statement
The accepted statement is reproduced exactly from the claim registry.

C-GOT-005 (exact conformal/log-shear decomposition): For three positive principal optical indices n_i and nbar=(n_1 n_2 n_3)^(1/3), the paper's additive deviations n_i-nbar are not generally trace-free; the witness (4,2,1) has nbar=2 and trace 1. Define instead lbar=(log n_1+log n_2+log n_3)/3 and E_i=log n_i-lbar. Then sum_i E_i=0, product_i exp(E_i)=1, and E_i=0 in the isotropic sector. This logarithmic unimodular variable separates conformal and shear spectral modes exactly and removes the double counting present in equation (58).

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: The three principal indices are positive real numbers and logarithms use their real branches.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.170.0` with provenance `campaigns/P248-sciama-gothic-metric-audit/adjudication.yaml`.

- `campaigns/P248-sciama-gothic-metric-audit/attempts/0002/result.yaml`
- `campaigns/P248-sciama-gothic-metric-audit/companion/exact_compatibility_checks.py`
- `campaigns/P248-sciama-gothic-metric-audit/reviews/C-GOT-001-C-GOT-006-review.md`
- `formal/SubstrateFramework/OpticalGothic.lean`
- `src/substrate_framework/optical_gothic.py`
- `tests/test_optical_gothic.py`
