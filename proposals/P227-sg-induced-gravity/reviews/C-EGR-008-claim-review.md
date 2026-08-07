---
description: Independent review of full-retardation breather radiation claim C-EGR-008
author: axis-review
created: '2026-08-07T06:59:00Z'
updated: '2026-08-07T06:59:00Z'
tags:
- substrate-framework
- claim-review
- radiation
category: decisions
confidence: working
status: active
---
# Review of C-EGR-008

## Claim Under Review

C-EGR-008 is a resolution-bounded full-retardation prediction from the same
exact stress and derived coupling used by C-EGR-007. It retains the retarded
phase across the entire line source, applies the exact TT angular norm and
weak-wave flux, and reports the normalized power, analytic small-amplitude
transform, small-width limit, harmonic tail, and one-cycle loss gate.

## Sourced Inputs

The review reads C-EGR-006/007, C-SG-012/017, the landed conditional retarded
ledger, the complete numerical module/tests, all P227 refinements, the primary
verifier, and the fresh Simpson implementation. It also checks the shared
trapezoid owner and the no-legacy-alias source audit.

## Independence

The fresh reviewer imports none of the P227 modules. It evaluates the exact
breather pressure directly, projects four harmonics with SciPy Simpson
quadrature, integrates each spatial retarded phase separately, and performs a
fresh direction integral. Its result `0.0002875595530006727` differs from the
refined package result by about `2.21e-12`. Adaptive quadrature separately
reproduces the analytic leading transform.

## Verification Status

The maximum verdict is
`numerical_refinement_and_independent_quadrature_verified`. The symbolic source
kernel is first cross-checked against the canonical C-SG-012 stress. Three
simultaneous refinements converge to `0.0002875595507880884`; the two finest
runs agree at displayed precision. The independent method, analytic leading
transform, `eta` asymptote, positive flux, harmonic decay, source-size ratio,
and backreaction fraction all pass. This is not promoted to exact verification
because the full nonlinear transform remains sampled quadrature.

## Sensitivity and Counterexamples

Spatial, temporal, angular, and harmonic refinements are independently
load-bearing. Even grid sizes, inadequate temporal sampling, invalid extent,
float frequency, and nonnumeric physical coefficients fail with context. A
slow quadrupole replacement is invalid here because the source-size frequency
ratio exceeds 33. Changing the upstream common SG scale changes absolute power
but leaves `P/(mu*c)` invariant, while changing `G` independently would fail
the shared-provenance test.

## Framework Compatibility

The source is exactly the C-EGR-007 stress component; no ad hoc quadrupole or
separate radiative source is introduced. The no-incoming boundary condition,
far-zone weak field, Isaacson averaging, and one-way approximation are explicit
premises. The 0.377 percent cycle-loss estimate bounds, but does not solve,
radiation reaction.

## Dependency and Consumer Replay

The claim depends on C-EGR-006/007 and reuses their identical gravity ledger
and source energy. Numerical ownership stays in `substrate_framework.numerics`.
The 27 focused tests, primary 20-check verifier after the link-map addition,
fresh six-check review, and 127 accepted dependency tests pass. Full repository
replay remains the final pre-PR gate.

## Competing Candidate Audit

Radiative numbers were opened after candidate C was selected without an
empirical comparator. Full retardation was retained because the structural
source-size test rejects the slow-source shortcut, not because it improved
agreement with an external datum.

## Four-Axis Decision

The numerical claim is ready for maintainer review at its bounded verdict.

- Verification: `numerical_refinement_and_independent_quadrature_verified`
- Review: `proposed_pending_maintainer`
- Compatibility: `compatible_extension_candidate`
- Epistemic: `proposed`
- Relationship: numerical consequence of C-EGR-006 and C-EGR-007

## Promotion Transaction

If accepted after its dependencies, promotion includes the numerical API,
tests, frozen resolution audit, claim registry, release, generated docs, and
accepted memory. The benchmark precision and assumptions must remain in the
registered statement.

## Continuation if Not Accepted

A failed convergence, independent method, source identity, boundary condition,
or backreaction gate returns to the radiative layer. It cannot be repaired by
rounding the benchmark or reverting to a slow quadrupole formula.

## Done Gate

The local source, refinement, independence, limit, provenance, consumer, and
debt gates pass at a numerical verdict. Maintainer authority and repository
materialization remain pending.

## Cross-References

See C-EGR-005/006/007, C-SG-012/017, P227 numerical and dependency audits, the
radiation module/tests, and both verifiers.
