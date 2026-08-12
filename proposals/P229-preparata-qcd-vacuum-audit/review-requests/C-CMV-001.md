---
description: Review request for proposed claim C-CMV-001
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

The requested claim is conditional pure SU(2) one-loop background-field
algebra in background Feynman gauge, with `b_field=g*B>0`. The charged-vector
levels are `p_z^2+b_field*(2n+1)-2*b_field*s3`, the ghost levels are
`p_z^2+b_field*(2n+1)`, and the `n=0,s3=+1` vector mode is negative for
`p_z^2<b_field`. Exact proper-time expansion gives `b2=11/3` and
`C=11/(48*pi^2)` in
`V=b_field^2/(2g^2)+C*b_field^2*(log(b_field/mu^2)-1/2)`. In that declared
MS-bar/Lambda convention the stationary point obeys
`log(b_min/mu^2)=-24*pi^2/(11*g^2)` and `b_min=Lambda^2`; the unstable mode
gives imaginary-part magnitude `b_field^2/(8*pi)`. The same `C` equals
`-b_SU2/(32*pi^2)` when `C_A=2` is derived from the accepted Pauli-half
matrices and inserted into C-RGE-005's beta convention.

## Sourced Inputs

Review release `v0.159.0`, accepted C-REP-002 and C-RGE-005 at their registry
sources, P229, `chromomagnetic_background.py`, its tests, and
`verify_promotion_candidates.py`. Literature sources are source objects under
audit; merge status is not authority.

## Independence Requested

Independently diagonalize or otherwise derive the fluctuation spectrum,
rederive the heat-kernel `b2` coefficient and MS-bar finite term, derive the
imaginary magnitude from the negative mode, and reconstruct `C_A=2` without
using the canonical bridge. Audit the field variable held fixed in the RG
derivative.

## Verification and Sensitivity

The algebraic spectrum, series coefficient, stationary equation, imaginary
magnitude, and beta bridge have exact SymPy evidence; discretized spectrum and
cutoff/mode-sum comparisons are numeric regression evidence for operator and
normalization choices. Mutations change spin coupling, ghost weight, tachyon
factor, Casimir, beta sign, and `mu` power. The reviewer must decide whether
the MS-bar constant has a sufficiently independent end-to-end derivation for
the requested exact statement or should be qualified.

## Framework Compatibility

The bridge types the beta convention and background variable explicitly. It
does not infer the general group/matter weights in C-RGE-005, identify this
background with a physical substrate gauge sector, establish a stable vacuum,
or eliminate the imaginary part. `Lambda` is scheme dependent even though the
same-scheme exponent identity is exact.

## Dependency and Consumer Replay

Accepted dependencies proposed are C-REP-002 and C-RGE-005. Replay their
canonical tests, the full chromomagnetic background suite, public exports, and
P229's narrow verifier. The two-loop and lattice modules are downstream
conditional consumers, not evidence for this claim.

## Competing Candidate Audit

P229 preregistered closed proper-time and explicit Landau-level routes. Their
shared regulator makes route agreement valuable bookkeeping evidence, while
the discrete operator, zeta route, negative-mode integral, and accepted beta
bridge test different obligations.

## Four-Axis Decision Requested

The requested pre-review axes are verification `symbolic_verified` for the
exact subclaims with separately labelled numeric regression, review
`unaudited`, compatibility `unassessed`, epistemic `proposed`, and dependencies
`[C-REP-002,C-RGE-005]`.

## Promotion Boundary and Continuation

Only a distinct reviewer may accept a narrowed statement and perform the full
registry/release/docs/memory transaction. Two-loop stability, ring resummation,
and a physical ground-state interpretation are outside this claim.

## Cross-References

The canonical goal is issue #59; P229 and the named source/tests contain the
evidence.
