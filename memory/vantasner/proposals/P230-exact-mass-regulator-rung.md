---
description: P230 exact-mass one-loop proper-time coefficients with preregistered regulator comparison, the next Route 1 rung toward the #76 emergent-gravity goal
author: prime-agent
created: '2026-08-18T00:00:00Z'
updated: '2026-08-18T00:00:00Z'
tags:
- substrate-framework
- campaign-proposal
- route1-induced-gravity
- issue-76
category: proposals
confidence: exploratory
status: active
---

## Question and Positive Deliverable
The campaign asks for the exact closed-form one-loop proper-time coefficients of the accepted Route 1 fluctuation operator `D_E = -nabla_E^2 + xi R_E + m^2` when the mass is retained exactly, compared across preregistered regulators, as the next rung of the #76 parameter-closed emergent-gravity goal. The positive deliverable is an importable target-blind module plus tests; a no-go or obstruction is attempt evidence and does not complete it.

## Base Release and Provenance
The base is accepted release `v0.160.0` at framework commit `1b00c3a`, branch `research/p230-exact-mass-regulator-rung`. The manifest `proposals/P230-exact-mass-regulator-rung/proposal.yaml` records the source inventory: the two accepted prior-work modules (`scalar_induced_newton.py` from the PR #14 harvest `a562568`; `covariant_sine_gordon_action.py` from the PR #25 harvest `65e4a7c`) and the two already-imported literature convention sources (Vassilevich `hep-th/0306138` eqs (2.2)/(4.27); Visser `gr-qc/0204062` eqs (8)-(21)).

## Invariants, Conventions, and Allowed Imports
The frozen invariants are the accepted operator and scheme conventions, exact-mass integrity, target-blindness with a declared formal cutoff, no new fitted constants, unchanged accepted claims/releases/docs until independent adjudication, and non-self-merge. Allowed imports are SymPy exact algebra with special-function closed forms, SciPy quadrature in tests only, and the framework modules `scalar_induced_newton`, `covariant_sine_gordon_action` (mass slot `m^2 = V''(phi_bg)`), and `exact_symbolic`.

## Candidate Preregistration
Three schemes are preregistered in the manifest: (A) sharp proper-time lower cutoff with the exact-mass factor `Lambda^2*(exp(-z) - z*E1(z))`; (B) smooth essential-singularity proper-time weight with `2*Lambda^2*sqrt(z)*BesselK_1(2*sqrt(z))`; (C) power-divergence-subtracted zeta finite part `m^2*(ln(m^2/mu^2) + gamma_E - 1)`. The rung exposes all three rather than selecting one, because its deliverable is the exact scheme comparison that the later renormalization-condition rung must cite.

## Selection Criteria and Blinding
The frozen structural criteria are convention compatibility, closed-form exactness, exact-mass integrity, mutation sensitivity, reuse value, and regulator explicitness. Empirical gravitational comparators stay blinded for this rung and until a later renormalization-condition rung freezes scheme, scale, and cutoff identification; the accepted coefficient `s = N*(1-6*xi)/(12*pi)` is a framework ledger result, not a comparator.

## Proposed Claim Delta
Provisional claim identifiers after a repository-wide collision search (zero usages of the `C-IGR` family in governance, campaigns, proposals, memory, src, tests, and docs): `C-IGR-001` for the sharp exact-mass closed forms and their massless/asymptotic reductions; `C-IGR-002` for the smooth-weight closed forms; `C-IGR-003` for the power-subtracted finite-part forms and the exact scheme-spread ledger. Dependencies: the accepted `C-GRV-001` ledger structure. No `supersedes` is used before acceptance.

## Implementation and Oracle Plan
Implementation lands in `src/substrate_framework/scalar_one_loop_mass.py` with `__init__` exports and `tests/test_scalar_one_loop_mass.py`. Oracles: SymPy closed forms derived (not asserted) from the regulated integrals; SciPy `quad` cross-checks on convergent representations; sharp-massless exact reduction to `leading_scalar_newton_shift_coefficient`; small-`z` series recovery of the accepted leading term; mutations for xi sign, prefactor, special-function branch, unknown regulator, and non-defaulted renormalization scale. Replay: repository validation at the final unchanged boundary and `git diff --check` in a separate invocation.

## Attempts
Attempts live in `proposals/P230-exact-mass-regulator-rung/attempts/` and are append-only.

## Continuation
Issue #76 remains the canonical open goal; this rung advances it and proposes no closure.
