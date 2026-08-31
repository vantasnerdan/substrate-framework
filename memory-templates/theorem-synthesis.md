# Theorem Synthesis Contract Template

Use for a campaign whose positive objective is one higher theorem composed from
accepted framework claims. Instantiate it in memory and create a matching
proposal manifest before substantive proof work.

```md
---
description: <one higher theorem to promote>
author: <agent-id>
created: '<ISO-8601>'
updated: '<ISO-8601>'
tags:
- substrate-framework
- theorem-synthesis
category: proposals
confidence: exploratory
status: active
---

## Exact Theorem Target

State one positive theorem precisely enough to encode in SymPy or Lean. Say
whether it is `core` or `interpretive`; an interpretive target names hypothesis
H in the statement itself. Explain the useful structural gap it closes. Do not
shrink the target to a tautology, definition, or isolated numeral merely to make
verification easy.

## Accepted Composition Boundary

List at least two distinct accepted claim IDs and quote only the exact portions
used. Record the pinned release and source modules. These atoms are inputs, not
objects of renewed acceptance review.

## Structural Gap

Name the implication, identity, elimination, convention map, or other glue that
is not already an accepted theorem. Explain what downstream question becomes
answerable when it closes.

## Assumptions and Exclusions

List every added assumption and explicit exclusion. For an interpretive
theorem, record the hypothesis label and statement and keep its consequences
out of the core dependency layer. A real declared hypothesis is not a defect;
prefer an honest conditional theorem to rejecting useful composition.

## Proof Route

Declare `target_kind: fixed_theorem` when the statement is fixed. One complete
proof route is sufficient; compare alternatives only when mechanisms are
actually competing or a second route materially reduces risk. Select SymPy for
exact algebraic glue or Lean for a finite formal theorem, and name the artifact
and entrypoint. For Lean, record imports, statement audit, axiom footprint, and
the mapping to the physics statement.

## Evidence Ledger

Separate each modality and its scope. Exact proof establishes the encoded
statement. Measurement, numerics, or simulation test applicability or empirical
adequacy and never retroactively become the proof. A support gap is not a
refutation of either accepted dependency.

| Method | Artifact | Exact scope | Verdict |
| --- | --- | --- | --- |
| SymPy or Lean |  | glue theorem |  |
| Measurement/numeric/simulation if useful |  | applicability or consequence |  |

## Promotion Transaction

Record the synthesized claim entry, individual review of the new theorem,
dependency closure, affected consumers, release pin, generated docs, and
accepted-memory synchronization. Do not re-review the accepted atoms.

## Proportional Validation

Record the exact theorem verifier, statement/axiom audit, affected tests and
consumers, GitNexus impact, and the `scripts/validate_changed.py` decision. An
additive leaf theorem may use scoped validation when no existing contract is
changed; use full replay for cross-cutting or uncertain impact. Do not repeat an
equivalent validation at the same unchanged boundary. Freeze the composition
boundary before review; adjacent defects become follow-up unless they directly
invalidate the new glue, and review prose never triggers another oracle run.

## Constructive Review

Use one substantive review and one correction check. Lead with the strongest
meaningful theorem the glue proves, then identify any unsupported extension and
the minimum repair. Recheck only corrected glue and affected consumers. Do not
reopen accepted atoms, audit every historical attachment, or commission another
review because prose, counts, or evidence roles changed.

## Attempts and Frontier

Preserve failed proof routes append-only and state the next materially different
route. Keep supporting implementation and locally complete lemmas as commits on
the synthesis campaign branch; do not open a rung or partial theorem PR.

## Done Gate

Close only when the exact higher claim is individually accepted, registered,
dependency-closed, pinned in a release, and free of hidden in-scope debt. An
explicit interpretive hypothesis or honest exclusion is not debt. If the
positive theorem remains unproved, continue materially different proof routes
until the scientific-exhaustion certificate in `AGENTS.md` passes. Only then may
the terminal campaign PR open, using `Advances` rather than `Fixes`.
```

Matching proposal fields:

```yaml
campaign_type: synthesis
target_kind: fixed_theorem
structural_gap: <gap closed by the higher theorem>
composition_dependencies: [C-SECTOR-001, C-OTHER-001]
candidates:
  - id: proof
    description: <complete SymPy or Lean glue route>
claims_proposed: [C-NEW-001]
```
