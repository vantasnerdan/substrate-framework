---
description: Independent review and acceptance decisions for C-WLN-001 C-WLN-002 C-WLN-003 C-LOR-001 and C-LOR-002 (P227 worldline harvest)
author: vantasner and codex
created: '2026-08-17T20:30:58Z'
updated: '2026-08-17T20:45:00Z'
tags:
- substrate-framework
- claim-review
- P227
- worldline
- lorentz
category: decisions
confidence: established
status: active
---
# P227 Claim Reviews

## Claims Under Review

C-WLN-001, C-WLN-002, C-WLN-003, C-LOR-001, and C-LOR-002 from the
P227 einbein-worldline-harvest campaign
(`campaigns/P227-einbein-worldline-harvest/`). Each review packet is
the canonical conditional exact classical statement derived from the
merged PRs #38 and #39 tutorials, exposed as a small reusable
ledger in `src/substrate_framework/relativistic_particle.py`,
`lorentz_orbits.py`, and `lorentz_little_groups.py`. The author
(`codex`) supplied the implementation and the review packets; this
file closes the closure campaign by recording the distinct reviewer
decisions.

## Independence and Sensitivity

The draft reviewer (`vantasner-review`) supplied prose rederivations. The
merger audit found false printed residuals and generator products in the
C-WLN-002 and C-LOR-002 drafts, corrected them, and added
`reviews/independent_worldline_lorentz_review.py`. That executable reviewer
derives the claims before importing the canonical helper APIs and passes 51
mutation-sensitive checks:

- **C-WLN-001**: derives the `e` Euler equation, solves the auxiliary
  constraint for the positive branch, substitutes back into the
  density, identifies the canonical covector momentum and inverse
  Legendre map, and verifies the pure-constraint Hamiltonian on a
  non-diagonal metric. Rejects wrong-sign mass-term, wrong index
  placement, negative branch, and `m=0` mutations.
- **C-WLN-002**: specializes `m=0` first, varies the resulting
  density with respect to `e` and the coordinates, reconstructs the
  Christoffel symbols from metric derivatives, and confirms the
  lower-index Euler residual. Rejects massive eliminator at `m=0`
  and wrong geodesic-force sign mutations.
- **C-WLN-003**: transforms the density and integration measure
  directly under reparametrization, integrates the local gauge rate
  ODE, and compares the massless and massive Weyl terms. Rejects
  wrong reparametrization weight, wrong constant-gauge rate, wrong Weyl
  weight, and identity / nontrivial `Omega` conflation.
- **C-LOR-001**: differentiates the displayed hyperboloid embeddings
  and pulls back the ambient metric by hand for both 2+1 and 3+1.
  Rejects sign-flipped ambient metric, wrong dimension, and angular
  swap mutations.
- **C-LOR-002**: reconstructs the generators from the Lorentz-algebra
  condition `M^T eta + eta M = 0` and the fixed-vector equations,
  then multiplies matrices directly. Rejects the relative-sign mutations
  `K1-J2` and `K2+J1` plus the wrong commutator sign. An overall generator
  sign is correctly recognized as remaining in the stabilizer.

All five reviews include a comparison of the independent
rederivation against the canonical helper; in each case the
independent result agrees exactly in the displayed repository convention.

## Four-Axis Decisions

All five claims are `symbolic_verified`, `accepted`,
`compatible_extension`, and `active`. None uses a supersession
edge. The scope ceilings are claim content:

- C-WLN-001: classical auxiliary-field identity, not a quantum or
  field-theoretic statement. The zero-mass branch is reserved for
  C-WLN-002.
- C-WLN-002: classical, local, conditional on supplied exact metric
  data; derives neither a field equation nor a massless particle
  representation.
- C-WLN-003: one-dimensional worldline action identity only, not a
  Maxwell/Yang-Mills conformal symmetry, photon count, graviton
  claim, or spacetime Weyl-gauge theory.
- C-LOR-001: induced metric on a single future unit timelike-vector
  orbit, dimensions 3 and 4 only. Not a timelike two-plane, tube
  worldsheet, or Lorentz-invariant ensemble.
- C-LOR-002: matrix Lie algebra only. No global topology, unitary
  representation, spin/anyon/continuous-spin classification,
  helicity, parity pairing, photon, graviton, or interaction claim.

## Source Decision

All five review packets pass independent rederivation, mutation
sensitivity, and framework-compatibility checks. No competing
candidate (B coordinate-function derivation, C document-local
duplication) overrules the chosen candidate-A ledgers. Proposal P227
transitions from `in_review` to `accepted`. The implementation
remains mathematically unchanged from PR #60. The closure transaction also
materializes the immutable campaign, accepted registry entries, release
`v0.160.0`, generated documentation, and accepted memory required to make
the five four-axis decisions authoritative.

## Cross-References

See P227, the five per-claim reviews, the P227 verifier
(`campaigns/P227-einbein-worldline-harvest/verify.py`, 14 checks
pass), the canonical modules
`src/substrate_framework/{relativistic_particle, lorentz_orbits,
lorentz_little_groups}.py`, and the canonical worldline/Lorentz
test suites. Accepted base `v0.159.0`; promoted release `v0.160.0`; source PR #60 (merged).
Canonical goal: issue #59; closure tracker: issue #74.
