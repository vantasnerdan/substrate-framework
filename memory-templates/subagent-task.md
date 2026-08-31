# Delegated Worker Contract Template

Use only when delegation is authorized. Give each worker one bounded route and a
disjoint write surface. The worker does not open a PR. Failure is returned as
attempt evidence with its mechanism and next route; it never masquerades as
completion, exhaustion, or permission for the parent campaign to stop.

```md
---
description: <single positive deliverable owned by this worker>
author: <agent-id>
created: '<ISO-8601>'
updated: '<ISO-8601>'
tags:
- substrate-framework
- delegated-task
category: efforts
confidence: working
status: active
---

## Positive Deliverable
State the exact object this worker must produce and the relevant success gates. No-go, residual, obstruction, or partial verifier is not an accepted substitute.

## Sourced Inputs
List exact release, claim ids, paths, equations, data, parent contract, and write surface. Provide facts rather than the parent agent's expected interpretation.

## Invariants and Boundaries
State conventions and accepted invariants. Name files the worker may edit and surfaces owned by others.

## Candidate and Method
Record the assigned candidate or independent review route, claim-appropriate
oracle, selection criteria, and comparator gate. Do not change foundations to
rescue the candidate. If numerical, supply the parent's passed analytic-closure
receipt, strongest exact result, named irreducible remainder, and frozen design
freedoms before naming the SciPy routine or other method, precision,
equations/data, discretization, tolerances, status gate, refinement plan, error
norm, and independent check. Without that receipt the worker may do only
explicitly exploratory hypothesis-generation or implementation debugging; it
cannot return scientific evidence or a candidate refutation.

## Steps
List dependency-ordered implementation, verification, sensitivity, and consumer checks.

## Attempts
Append failed routes with commands, diagnoses, and next materially different attempts. Continue until the deliverable and debt gates pass.

## Validation
Record exact commands for the positive claim, mutation/counterexample behavior, refinement or independent route, and affected consumers.

## Debt Ledger
Record and discharge every new import, assumption, parameter, residual, convention conflict, or broken consumer.

## Results
Return the positive result, changed paths, status axes, reproduction commands, and claim delta. Failed attempts stay here while the contract remains active.

## Done Gate
Confirm the deliverable and applicable success conditions from `AGENTS.md`; otherwise
name the next route as its positive contribution to resolution — the object it
constructs, the question it closes, or the distinction it establishes, phrased as what
becomes true in the record when it succeeds. A refutation counts when stated through
its mechanism; pure avoidance or risk reduction is an unfinished answer. A worker
may finish its assigned route boundary, but the parent executor must reconcile
the result and continue the active obligation through method repair,
representation change, and alternative concepts until campaign success or
certified exhaustion.

## Cross-References
Link parent contract, proposal, claims, artifacts, and sibling awareness paths.
```
