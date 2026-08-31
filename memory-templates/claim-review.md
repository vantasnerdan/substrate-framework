# Claim Review and Promotion Template

Use one instance for each claim proposed for acceptance or changed accepted
statement. Do not create one per evidence attachment or theorem entrypoint; use
`evidence-attachment-review.md` for a lightweight evidence-role audit.

Begin every section with a plain-prose sentence so the memory index can disclose
it correctly.

```md
---
description: Constructive review of <claim-id>
author: <reviewer-id>
created: '<ISO-8601>'
updated: '<ISO-8601>'
tags:
- substrate-framework
- claim-review
category: decisions
confidence: working
status: active
---

## Claim and Positive Role
Quote the exact statement, hypotheses, quantifiers, regime, conventions, and
the useful framework question it answers. A claim that only restates a
definition, literal, or vacuous special case has not met the meaningfulness
floor.

## Frozen Transaction
Record the base/head or tree hash, exact claim delta, changed implementation or
evidence records, accepted dependency propositions actually used, affected
consumers, and existing validation receipt. Unchanged accepted dependencies and
adjacent corpus records are outside this review.

## Strongest Supported Positive Statement
State the strongest meaningful result supported by the sourced artifacts. If it
differs from the proposal, preserve as much useful scope as the evidence permits
and identify the smallest honest change: quantifier correction, explicit
hypothesis, core/interpretive split, or evidence-role relabel.

## Evidence Map
Classify what each coherent evidence group actually contributes. One row may
cover related entrypoints proving the same proposition.

| Evidence | Proposition established | Role: exact proof / corroborating subclaim / regression / applicability / provenance | Bridge to claim | Limit |
| --- | --- | --- | --- | --- |
|  |  |  |  |  |

An attachment need not prove the parent claim in full. It must not be credited
beyond its role. Missing a bridge means “not established by this attachment,”
not “the claim is false.”

## Oracle Audit
Inspect the one strongest practical oracle for the load-bearing proposition.
For numerical or simulation evidence, first inspect the analytic-closure
receipt: exact equations and variations, scaling, bounds and limits, asymptotic
structure, strongest non-numerical conclusion, named numerical remainder, and
frozen design freedoms. Confirm that the computed predicate answers that
remainder rather than defining a different question. Then record only the
applicable statement/axiom audit, mutation or counterexample, limit, refinement,
independent rederivation, and solver/error facts needed to exclude a false
green. Pre-gate exploratory samples cannot support acceptance or refutation.
Reuse a passing receipt while its code and inputs are unchanged; do not validate
this review record.

## Findings
Classify each finding once. A current blocker must be a counterexample or
contradiction, an absent/circular load-bearing step, a declared dependency that
does not supply what is used, or an affected-consumer failure.

| Finding | Direct evidence | Blocking in boundary / minimum correction / follow-up | What would resolve or overturn it |
| --- | --- | --- | --- |
|  |  |  |  |

Reserve `refuted` for an explicit falsifier under the stated hypotheses. Use
`unverified`, `qualified`, or a narrower evidence role for missing support.
Every narrowing finding carries its upgrade path in the same row: what specific
evidence would restore or extend the stronger form. A finding without an upgrade path
is incomplete, because it hands the next agent a smaller question instead of a next
step. If prior reviews of this claim already reduced its scope without new contrary
evidence arriving in between, say so here and escalate to the user rather than narrow
a third time.

## Compatibility and Consumers
Record assumptions, imports, units, conventions, invariants, and the result of
the impact-bounded consumer replay. Declared hypotheses and unfinished parent
work are not debt. Record only defects introduced inside this transaction.

## Four-Axis Decision
State the decision without inflating one axis from another.

- Verification:
- Review:
- Compatibility:
- Epistemic:
- Relationship:
- Strongest accepted or proposed statement:

## Promotion Transaction
List only the registry, implementation/test, campaign, release, generated-doc,
memory, migration, and validation changes actually required by this claim. Use
the existing content-addressed validation receipt when the boundary is unchanged.

## Correction Check
After requested changes, check only the corrected statements, altered evidence
roles, and directly affected dependency or consumer edges. This is the final
review pass, not a new corpus audit. Record `not needed` when no correction was
requested.

## Result and Frontier
Lead with the positive result retained. If the claim is not accepted, state one
decisive missing construction or test and leave the parent objective open.
Adjacent findings appear once as follow-up and do not expand this transaction.

## Cross-References
Link proposal, claim, dependencies, evidence, consumers, release, validation
receipt, and parent research arc.
```
