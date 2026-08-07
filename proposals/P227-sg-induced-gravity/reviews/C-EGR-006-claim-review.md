---
description: Independent review of nonlinear TEGR-to-Einstein claim C-EGR-006
author: axis-review
created: '2026-08-07T08:12:00Z'
updated: '2026-08-07T08:12:00Z'
tags:
- substrate-framework
- claim-review
- general-relativity
category: decisions
confidence: working
status: active
---
# Review of C-EGR-006

## Claim Under Review

C-EGR-006 states that the C-EGR-005 collective action is nonlinear GR, not a
linear spin-two analogy. In the fixed mostly-plus convention it constructs
Weitzenbock torsion, the Levi-Civita scalar, and the boundary divergence from
one nondegenerate coframe and requires `R+T-B=0`. Therefore `-mu*T/2` is
`+mu*R/2` modulo the displayed boundary term, fixing `kappa=1/mu` and
`G=c^4/(8*pi*mu)` for minimally coupled conserved matter.

## Sourced Inputs

The review reads C-EGR-005, C-MED-003, C-STG-001 for the accepted action/sign
convention, the primary teleparallel theorem, P227 attempts and evidence, the
complete coframe implementation and tests, the primary verifier, and the
fresh reviewer. The rejected horizon route and accepted induced-gravity
ceiling are checked only as alternatives and negative guards.

## Independence

The fresh review does not import the P227 geometry engine. It derives the
exponential-FLRW values `T=6H^2`, `B=18H^2`, and `R=12H^2` separately, then
checks `R+T-B=0`. It matches `mu/2` to `1/(2*kappa)` and derives the Newton
coefficient directly. The primary route instead constructs every tensor index
contraction and Christoffel component from a supplied coframe.

## Verification Status

The identity and coefficient earn `symbolic_verified_with_sourced_theorem`.
The general exact identity is a sourced differential-geometric theorem, while
the implementation independently constructs both sides for arbitrary supplied
coframes and closes two structurally distinct exact tests: nonlinear FLRW and
a local Lorentz frame on flat geometry. No perturbative truncation is used in
the action equivalence.

## Sensitivity and Counterexamples

Replacing `1/4` by `1/3` produces residual `-H^2/2` on exponential FLRW.
Changing the boundary sign or the teleparallel action sign breaks the positive
Einstein coefficient. Degenerate coframes and nonexact inputs fail validation.
Nonlinear functions of `T` are outside the claim because they would not differ
from Einstein-Hilbert by only the stated boundary term and may carry extra
modes.

## Framework Compatibility

The claim preserves the accepted mostly-plus Einstein sign and lands exactly
on the C-STG-001 gravitational normalization with a derived rather than
supplied `kappa`. The flat spin connection or covariant coframe/connection
pair, admissible boundary variation, minimal matter coupling, and conservation
are explicit premises. No cosmological term is derived or silently set from an
observation.

## Dependency and Consumer Replay

The graph is C-MED-003 -> C-EGR-005 -> C-EGR-006. C-EGR-007 and C-EGR-008 are
the direct positive consumers. The landed linearized response ledger receives
exactly the derived `G,c` and returns the same `1/mu` stress coefficient. All
focused and accepted dependency tests pass; no existing consumer changes.

## Competing Candidate Audit

Candidate C was selected structurally before static or radiative numbers. The
local Lorentz theorem repairs P226's precise teleparallel rejection by fixing
the ratios, while C-EGR-005 repairs the remaining scale ambiguity. Induced and
horizon routes do not enter this proof.

## Four-Axis Decision

The nonlinear equivalence is ready for maintainer review without authority
self-assignment.

- Verification: `symbolic_verified_with_sourced_theorem`
- Review: `proposed_pending_maintainer`
- Compatibility: `compatible_extension_candidate`
- Epistemic: `proposed`
- Relationship: new theorem depending on C-EGR-005 and C-MED-003

## Promotion Transaction

Acceptance must register the claim after C-EGR-005, retain the exact APIs and
tests, issue a release, regenerate governance products, and synchronize
accepted memory. The proposal does not perform those authority changes.

## Continuation if Not Accepted

An exact sign, boundary, connection, or coefficient failure returns to
candidate C. A linearized agreement or numerical metric cannot repair a failed
nonlinear identity.

## Done Gate

The positive object, exact identity, coefficient provenance, mutations,
consumers, and debt gates pass locally. Maintainer acceptance and governance
materialization remain pending.

## Cross-References

See P227, C-EGR-005, C-MED-003, C-STG-001, the literature audit, both
verifiers, `sine_gordon_teleparallel.py`, and C-EGR-007/C-EGR-008.
