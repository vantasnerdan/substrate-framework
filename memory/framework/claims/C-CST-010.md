---
description: Accepted framework claim C-CST-010
author: framework-registry
created: '2026-09-05T03:16:39+00:00'
updated: '2026-09-05T03:16:39+00:00'
tags:
- substrate-framework
- accepted-claim
- C-CST-010
category: claims
confidence: established
status: active
---
# C-CST-010

## Statement
The accepted statement is reproduced exactly from the claim registry.

For an integrable finite-range orientation energy of complete unresolved
states with an independent local product-Haar law, independent coherent
frame shifts leave its averaged energy invariant. Consequently its
coherent conservative torque and couple stiffness vanish, while quadratic
fluctuation energy need not vanish. Uniform one-point marginals alone do
not imply the result: correlated phases supply an explicit counterexample.
A no-retained-coherent-angular-current closure is separately stated for
the dynamical no-spin limit; static isotropy alone does not remove memory
or an angular current.

Exact proof:0058's Haar change of variables and canonical phase integration.
Corroboration: the already executed200000-sample signed-response Monte Carlo
with its declared seed, error model and bias mutation, at its narrower
first-moment scope. This does not assert that the coherent law of C-CST-009
is product Haar or that stationary Gaussian Euler isotropy by itself
erases its positive locking. It states the separately declared contrast
ensemble premise and its consequence, without deleting fluctuations.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: Complete unresolved orientation states have an independent local product-Haar law; the finite-range orientation energy is integrable. The theorem is not inferred from uniform one-point marginals or Gaussian Euler isotropy., No retained coherent angular current is a separate explicit dynamical closure. Static conservative invariance does not remove quadratic fluctuations, viscosity, transient memory or finite-frequency response., Exact proof: Haar change of variables and the correlated-phase counterexample. Executed signed-response Monte Carlo is narrower finite-sample first-moment corroboration, not proof of a dynamical no-current law.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.175.0` with provenance `campaigns/P251-cosserat-from-vortex-euler/adjudication.yaml`.

- `proposals/P251-cosserat-from-vortex-euler/attempts/0105/claim-transaction.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0108/C-CST-010-review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0108/shared-scientific-review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0058/haar-contrast.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0058/independent-review.md`
- `proposals/P251-cosserat-from-vortex-euler/verify_cst006.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0025/stdout.txt`
- `src/substrate_framework/micropolar.py`
- `tests/test_micropolar.py`
