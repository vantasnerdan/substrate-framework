---
description: Independent review of derived-coupling static breather claim C-EGR-007
author: axis-review
created: '2026-08-07T08:14:00Z'
updated: '2026-08-07T08:14:00Z'
tags:
- substrate-framework
- claim-review
- weak-field
category: decisions
confidence: working
status: active
---
# Review of C-EGR-007

## Claim Under Review

C-EGR-007 composes the exact C-SG-017 breather energy and C-SG-012 conserved
stress with C-EGR-006, and nothing else, to produce the mostly-plus static weak
field. It fixes `M=E/c^2`, derives `G`, and states the exact profile-scale
compactness `4*eta^2/pi` and potential `-2*eta^2/pi`.

## Sourced Inputs

The review reads C-SG-012, C-SG-017, C-EGR-005/006, the landed fiber source and
linearized response modules, the new static composer and tests, P227 evidence,
and both verifiers. No Schwarzschild metric or empirical mass is supplied as
an input.

## Independence

The fresh reviewer uses only `E=16*eta*mu*ell`, `M=E/c^2`, and
`G=c^4/(8*pi*mu)`. It independently eliminates every intermediate scale to
obtain `r_s=4*eta*ell/pi` and compactness `4*eta^2/pi`. The proposed static API
is not imported.

## Verification Status

The maximum verdict is `symbolic_verified`. The source mass, Newton
coefficient, Schwarzschild radius, compactness, Green-function sign, `g_00`,
and spatial metric factor are exact SymPy expressions. The `eta=0.03` decimal
is regression display only; it does not select or normalize the claim.

## Sensitivity and Counterexamples

A common SG coefficient multiplier raises the source mass and lowers `G`
reciprocally, leaving the dimensionless field unchanged; this is the required
upstream scaling, not a fitted cancellation. Supplying a separate mass or `G`
would break provenance. Wrong mostly-plus signs fail the accepted Green
function ledger. Large `eta` can violate the weak-field premise, so the claim
reports compactness rather than generalizing the linearized result.

## Framework Compatibility

The composition preserves C-SG-017 energy units, the `(c*t,x)` C-SG-012 stress
convention, transverse-delta source localization, and C-EGR-006's derived
coupling. It adds only an observation radius and explicitly requires it to be
outside a weak source.

## Dependency and Consumer Replay

The claim depends on C-EGR-006, C-SG-012, and C-SG-017. C-EGR-008 consumes the
same stress and gravity ledger. The exact source, linearized response, and new
static tests pass with all accepted dependency tests; no existing caller is
modified.

## Competing Candidate Audit

The static number was opened only after candidate C had passed its geometry,
ratio, and scale gates. It therefore confirms a selected construction rather
than choosing among constructions. Direct Einstein-scalar and entropy-derived
couplings are absent.

## Four-Axis Decision

The exact static consequence is ready for maintainer adjudication.

- Verification: `symbolic_verified`
- Review: `proposed_pending_maintainer`
- Compatibility: `compatible_extension_candidate`
- Epistemic: `proposed`
- Relationship: consequence of C-EGR-006, C-SG-012, and C-SG-017

## Promotion Transaction

Promotion follows C-EGR-005/006 and includes the static API/tests, registry,
release, generated records, and accepted memory. None are self-promoted here.

## Continuation if Not Accepted

If source localization, mass conversion, sign, or weak-field control fails,
repair that exact layer. A numerical benchmark cannot substitute for it.

## Done Gate

The local exact, source, sign, scale, limit, consumer, and debt gates pass.
Authority and materialization remain pending.

## Cross-References

See C-SG-012, C-SG-017, C-EGR-005/006/008, P227, the static module/tests, and
both verifiers.
