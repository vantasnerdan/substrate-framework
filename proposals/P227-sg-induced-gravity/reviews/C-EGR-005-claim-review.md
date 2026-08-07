---
description: Independent review of compact sine-Gordon coframe claim C-EGR-005
author: axis-review
created: '2026-08-07T06:59:00Z'
updated: '2026-08-07T06:59:00Z'
tags:
- substrate-framework
- claim-review
- teleparallel
category: decisions
confidence: working
status: active
---
# Review of C-EGR-005

## Claim Under Review

C-EGR-005 is a constructed microscopic-to-collective map, not a statement that
the accepted one-dimensional scalar alone necessitates gravity. It declares a
nondegenerate four-link cell at the accepted SG length `ell`, maps the physical
link covectors to a coframe and mostly-plus metric, and gives every coframe
torsion channel an explicit compact spectral phase. Local Lorentz invariance
selects the TEGR quadratic ratio class; the one-channel-per-cell SG density
fixes its last multiplier. The exact chord map proves the microscopic density
is `-mu*T/2` without a regulator, observed `G`, or bare gravitational term.

## Sourced Inputs

The review reads accepted release `v0.159.0`, C-MED-003, C-GRV-001 as a
negative guard, the P227 freeze and all six attempts, the full candidate
comparison, the Itin-Hehl-Obukhov primary theorem, the complete proposed
teleparallel module and tests, the primary verifier, the fresh reviewer, the
dependency and source audits, and rejected P226 only as historical evidence.

## Independence

The fresh reviewer imports none of the proposed P227 modules. It reconstructs
the 24-component torsion tensor from antisymmetric local channels, derives the
three invariant quadratic form and its Hessian afresh, and obtains the exact
eigenvalue multiplicities `{-2:3,-1:8,-1/2:1,1/2:3,1:8,2:1}`. It separately
rederives the cosine/chord identity. The proposal's constitutive matrix,
spectral basis, helper constants, and coupling ledger are not reused.

## Verification Status

The maximum local verdict is `symbolic_verified`. The link-to-coframe metric
and measure are exact. The constitutive Hessian is symmetric and full rank,
its orthonormal spectral reconstruction is exact, every microscopic phase is
`2*pi` periodic, and the aggregate compact/action residual is zero. No
unevaluated integral, floating input, empirical comparator, or unresolved
symbolic condition enters the claim.

## Sensitivity and Counterexamples

Changing the SG cell length in the chord map makes the exact density identity
fail. Changing any TEGR invariant ratio makes the nonlinear identity fail on a
nontrivial FLRW coframe. A singular link matrix is rejected rather than called
a metric. Float phases and link covectors are rejected by exact APIs. Arbitrary
phase arrays are not claimed to be coframe torsion: admissible microscopic
link configurations must satisfy the plaquette integrability and Bianchi
constraints, which prevents the map from silently identifying unrelated
objects.

## Framework Compatibility

The construction preserves C-MED-003 coefficient meanings, mostly-plus
signature, exact-input discipline, the common SG scale direction, and the
accepted claim ceilings. Its new structural assumptions are explicit: four
cell-link covectors, one channel per independent torsion component and cell
area, parity-even local Lorentz symmetry, and coframe-integrable plaquettes.
They define a new microscopic model rather than reinterpret C-SG-011 or a
Gordon metric.

## Dependency and Consumer Replay

The only accepted positive dependency is C-MED-003. Direct consumers are
C-EGR-006, the new static/radiative composers, package exports, focused tests,
and P227 verifiers. The 127-test accepted dependency replay and 27-test P227
suite pass. Existing modules are unchanged, and no migration unit is created.

## Competing Candidate Audit

Five candidates were frozen before outputs. Entropy lacked a quantum state
law; induced gravity retained baseline and regulator freedom; Gordon supplied
kinematics; direct Einstein-scalar composition retained `kappa`. Candidate C
wins because an independent frame symmetry fixes the invariant ratios and SG
cell normalization fixes the remaining scale before any benchmark is opened.

## Four-Axis Decision

The claim is technically ready for maintainer adjudication but is not
self-promoted by this review.

- Verification: `symbolic_verified`
- Review: `proposed_pending_maintainer`
- Compatibility: `compatible_extension_candidate`
- Epistemic: `proposed`
- Relationship: new construction depending on C-MED-003

## Promotion Transaction

If accepted, promotion must add C-EGR-005 to the registry, retain the
importable implementation/tests and immutable P227 record, create a release,
regenerate docs and accepted memory, update any source dispositions, and run
the integrated governance validator. None of those authority edits are made
in this proposal branch.

## Continuation if Not Accepted

If the declared cell channel count, local Lorentz premise, or coframe
integrability is rejected as insufficiently upstream, the claim returns to
candidate C. A repair must construct that missing link dynamics explicitly;
neither a fitted overall coefficient nor candidate-B regulator may substitute.

## Done Gate

The claim-specific exact, mutation, framework-fit, and implementation gates
pass with no open code debt. Authority acceptance and materialization remain
pending, so the governance done gate is not yet claimed.

## Cross-References

See P227, C-MED-003, C-GRV-001, the candidate comparison, proposed claims,
`sine_gordon_teleparallel.py`, its tests, both verifiers, and C-EGR-006.
