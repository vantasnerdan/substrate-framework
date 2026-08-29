---
description: Complete optical ADM pullback theorem for Einstein-scalar fields and massive/massless worldlines
author: codex
created: '2026-08-29T06:30:00+00:00'
updated: '2026-08-29T07:06:51+00:00'
tags:
- substrate-framework
- theorem-synthesis
- optical-metric
- gothic-metric
category: proposals
confidence: established
status: archived
---

## Exact Theorem Target

The core target `C-GOT-006` states that a positive lapse, a symmetric
positive-definite spatial optical tensor, and a real material flow define a
complete algebraic ADM metric chart; pulling the accepted Einstein-scalar
action through this bijection preserves the full Euler equations, and using the
same metric in the accepted worldline claims yields the massive square-root and
Hamiltonian sector plus the massless null and affine-geodesic sector. The
theorem closes the paper's missing
tenth-field and action-variation gaps without claiming a microscopic ether
origin, a selected coupling, or empirical equivalence.

## Accepted Composition Boundary

The pinned boundary is release `v0.169.0` at
`443f1bd6edf0c1997f42b05c3b6719142e7a11c1`. `C-STG-001` supplies the exact
mostly-plus Einstein-scalar action, its field equations, and on-shell stress
conservation. `C-WLN-001` supplies the massive einbein/square-root worldline
identity for any supplied nondegenerate metric. `C-WLN-002` supplies the
massless null constraint and affine geodesic equation at constant einbein.
`C-GOT-002`, reviewed first in this transaction, supplies the complete optical
ADM bijection. None of the accepted atoms is reopened.

## Structural Gap

No accepted theorem currently maps a complete set of continuum optical fields
onto all Lorentzian metrics adapted to a spacelike foliation and carries both
the covariant field action and the two worldline sectors through that same map.
Closing the gap makes a precise conditional infrared optical-gravity theory
available while exposing why the paper's determinant-slaved nine-field map
cannot claim equivalence to general GR.

## Assumptions and Exclusions

The lapse is positive, the spatial tensor is real symmetric positive definite,
the flow is real, and variations are compact or include the appropriate
gravitational boundary term. The theorem excludes microscopic substrate
dynamics, derivation of a massless spin-two excitation, a physical constitutive
stress, impedance selection, a value of Newton's constant, and empirical or
ontological identification of the fields with an ether. These exclusions leave
the exact conditional field-redefinition theorem intact.

## Proof Route

The target kind is `fixed_theorem`. SymPy proves the block determinant,
inverse, reconstruction, and nonsingular component Jacobian and checks the
chain-rule implication. Lean encodes the finite algebraic and logical glue,
imports no target theorem, contains no proof escapes, and records its axiom
footprint. A second mechanism is unnecessary because the statement is a fixed
field-redefinition theorem; numerical SPD round trips are regression evidence,
not the proof.

## Evidence Ledger

The evidence modalities remain separate and will be filled at the frozen proof
boundary.

| Method | Artifact | Exact scope | Verdict |
| --- | --- | --- | --- |
| SymPy | `campaigns/P248-sciama-gothic-metric-audit/companion/exact_metric_checks.py` and package tests | block bijection, Jacobian, and chain-rule glue | passed: 31 exact campaign checks across the modular exact suite |
| Lean | `formal/SubstrateFramework/OpticalGothic.lean` | finite determinant/sign/compatibility and composition logic | passed: pinned build, no proof escapes, `commonMetricComposition` has no axioms |
| SciPy/NumPy | `campaigns/P248-sciama-gothic-metric-audit/companion/numerical_checks.py` | SPD round-trip, SI applicability regression, forced-wave convergence, mutation sensitivity | passed: 18 checks, including eight SPD round trips below `2e-13` and observed second-order convergence |

## Promotion Transaction

The promotion transaction reviews `C-GOT-006` individually after `C-GOT-002`
and the unchanged accepted atoms are pinned. The accepted entry will use
`category: synthesized`, `layer: core`, and composition metadata naming the
SymPy or Lean glue entrypoint. Direct consumers, release membership, generated
documentation, and accepted memory will be synchronized in the same promotion
transaction.

## Proportional Validation

The new module and formal surface are initially additive, but the transaction
also changes claim and release governance, so the final validation boundary is
full unless GitNexus and `scripts/validate_changed.py` prove a narrower
append-only leaf case remains valid. The exact verifier, focused tests, Lean
check, accepted dependency tests, changed-file decision, and `git diff --check`
will be recorded once at the frozen head.

## Constructive Review

One substantive review will test the exact higher statement against the raw
proofs and accepted atom scopes. If it finds a blocker, the review will first
state the strongest useful theorem, then the unsupported extension, minimum
repair, and restoring evidence. One correction check will inspect only that
repair and directly affected edges.

## Attempts and Frontier

The source paper's determinant-slaved map is the first route and is expected to
establish a constrained nine-field submanifold rather than the target
bijection. Restoring an independent lapse is the representation repair. If the
pullback Jacobian or action chain fails, the next route is a first-order
Palatini optical action; it is not needed unless the second-order fixed theorem
is obstructed.

## Done Gate

The synthesis closes only when `C-GOT-006` is proved, individually accepted,
dependency-closed, registered, pinned in a release, represented in generated
documentation and memory, and free of hidden in-scope debt. A conditional
infrared theorem with explicit exclusions is success; an unproved material or
ontological interpretation is not silently included.

## Final Promotion

The fixed theorem was proved and accepted as `C-GOT-006` after independent
review identified and the bounded correction check closed two statement-level
scope defects: the point-map coordinates are `q=(N,v,gamma)` with `v=V/c`, and
`C-WLN-001` supplies only its accepted massive square-root/Hamiltonian sector.
The accepted dependencies remain unchanged inputs, the exact glue is recorded
in the campaign, and the theorem is pinned in release `v0.170.0` alongside its
five new supporting claims.
