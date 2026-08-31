# Numerical Check Template

Instantiate whenever a numerical check produces evidence worth remembering — a
technique that worked, a failure mode worth naming, a number another attempt will
cite. It is a lab notebook page, not a compliance form, but any record cited as
scientific evidence must retain Question, Eligibility and Objective Bridge,
Method, What Was Seen, error budget/refinement, scoped Reading, and Artifacts.
Only a provenance-only note may omit those evidence sections. Store the filled
record in memory under `efforts` (or link it from the attempt's result file).
When the check leaves an open question, state the next step as its positive
contribution to resolution—the object it constructs, question it closes, or
distinction it establishes.

A production record also requires a passed analytic-closure receipt. A pre-gate
sample uses `numerical_role: exploratory_only`; it may generate a hypothesis or
debug an implementation, but it cannot support a claim or anti-claim, set later
thresholds, or be relabelled as production evidence.

Begin every section with a plain-prose sentence. Inline code, a table, or a list
does not satisfy the memory index's first-content disclosure contract.

```md
---
description: <what question this check answers and what it found, in one line>
author: <agent-id>
created: '<ISO-8601>'
updated: '<ISO-8601>'
tags:
- substrate-framework
- numerical-check
category: efforts
confidence: working
status: active
obligation: <obligation-node-id>
license_ids: []
analytic_closure_receipt: <path or exploratory_only>
numerical_remainder: <one proposition left after algebra/calculus/theorems>
numerical_role: <exploratory_only | numeric_evidence | simulation_evidence>
computed_predicate: <literal predicate evaluated>
proposition: <mathematical proposition implied>
maximum_verdict: <typed maximum verdict>
evidence_scope: <numeric evidence | unresolved | representation scoped>
---

## Question

What was being computed, and why it was delicate — usually because the quantity of interest is orders of magnitude below the dominant scale.

## Analytic Eligibility and Objective Bridge

Record the mathematical object, symmetry or conservation license, ensemble,
admissible representation, observable, and background/branch identity. State
the exact equations and complete variations, non-dimensional groups and scaling,
analytic bounds and limits, asymptotic operator or continuum threshold, and the
strongest result already fixed without discretization. Link the receipt, name
one residual proposition, explain why the current analytic ladder does not
decide it, and state the computed predicate, proposition it implies, maximum
verdict, parent obligation it advances, and nearby claims it cannot decide. If
an upstream license or analytic receipt is absent, classify the check as
`exploratory_only` or representation-scoped and return to construction rather
than using precision to force a verdict.

## Method

How it was computed: show how the numerical predicate is compiled from the
analytic remainder, then record the formulation chosen, frozen design freedoms
(chart or ansatz, box, boundary condition, mesh/basis, fit form, tolerance),
rejected alternatives when instructive, evaluator, background forward-error
estimate, branch-identity observables, constrained tangent space or gauge, full
error budget, and the execution context worth knowing later (thread pins,
invocation path, versions) if results sit near the 1e-13-relative level where
runner settings matter.

## What Was Seen

The numbers and behaviour actually observed — including how values moved under changes of mesh, domain, quadrature, or method, since that movement is often the most informative part.

## Error Budget and Refinement

Record background forward error and propagation, branch identity, representation
and constraint projection, operator construction, mesh/domain/basis crossed
refinement, boundary conditions, quadrature, solver/eigenpair residual,
roundoff/execution floor, observable normalization, and the required signal
margin. Vary mesh at fixed domain and domain at fixed mesh; co-varying both is
not a continuum test.

## Reading

What the observations mean in small-ratio terms: which cross-checks agree, whether the soft direction is bulk or boundary, what remains hypothesis versus established. Name failure modes plainly when they occurred; a named mechanism is the reusable part.

## Artifacts

Paths to scripts, data files, and logs.
```
