---
description: Review request for proposed claim C-CMV-003
author: codex
created: '2026-08-12T11:36:38+02:00'
updated: '2026-08-12T11:36:38+02:00'
tags:
- substrate-framework
- claim-review-request
category: decisions
confidence: working
status: active
---

## Claim Under Review

The requested claim is a narrow sector and root decomposition conditional on
C-CMV-001. In the stated SU(2) background-Feynman-gauge proper-time kernel, the
spin-up and spin-down transverse sectors each contribute `11/6` to the `b2`
coefficient, the longitudinal/timelike pair contributes `-1/3`, and the
complex ghost contributes `+1/3`; the latter cancel exactly, leaving total
`11/3` in the transverse sectors. For the canonical accepted SU(3)
`T3=lambda3/2` generator, the three positive off-diagonal root magnitudes are
exactly `(1,1/2,1/2)`, so an SU(3) T3 background has three corresponding
charged sectors. Each inherits a lowest-mode negative eigenvalue for nonzero
field. This is a quadratic one-loop instability statement only.

## Sourced Inputs

Review release `v0.159.0`, accepted C-LIE-001 at source, proposed C-CMV-001,
P229, `chromomagnetic_sectors.py`, `su3.py`, focused tests, and the P229 narrow
verifier.

## Independence Requested

Expand each unsimplified sector kernel independently, keeping vector and ghost
weights separate. Derive adjoint root charges from commutators with canonical
`T3`, not from the expected tuple or literature prose.

## Verification and Sensitivity

Author evidence is exact SymPy series and canonical matrix differences. A
ghost-weight mutation changes the total; replacing T3 by canonical T8 changes
the root pattern. The reviewer should ensure sector coefficients use one common
orbital normalization and confirm that root magnitudes, charge multiplicities,
and complex-vector counting are not conflated.

## Framework Compatibility

Proposed dependencies are C-CMV-001 and C-LIE-001. C-LIE-001 supplies only the
generator normalization; it does not create a physical color sector. The claim
does not prove global stability, two-loop behavior, a ring-resummed vacuum,
helicity content, or a substrate identification.

## Dependency and Consumer Replay

Replay canonical SU(3) generator tests, chromomagnetic background/sector tests,
the P229 narrow verifier, and public imports. The mixed `stability_verdict`
report contains excluded two-loop prose and is not the promotion oracle.

## Competing Candidate Audit

P229 preregistered closed-kernel and explicit-mode routes. Sector-by-sector
expansion plus generator commutators is selected because it exposes
multiplicities and conventions; a hard-coded root tuple is the rejected island.

## Four-Axis Decision Requested

The proposed axes are verification `symbolic_verified`, review `unaudited`,
compatibility `unassessed`, epistemic `proposed`, and dependencies
`[C-CMV-001,C-LIE-001]`.

## Promotion Boundary and Continuation

A distinct reviewer must adjudicate this claim separately from C-CMV-001. If
the SU(3) charged-vector multiplicity requires more representation-algebra
evidence, narrow to the exact root statement and leave the rest active.

## Cross-References

The canonical goal is issue #59; P229 contains all proposal evidence.
