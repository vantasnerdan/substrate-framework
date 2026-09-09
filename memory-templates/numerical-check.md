# Numerical Evidence Achievement

Use for a reusable numerical result or diagnosed method. Read the physics oracle
reference and applicable small-ratio skill. Earn scientific evidence with the
Question, Analytic Construction and Objective Bridge, Method, Observations,
Error Budget/Refinement, Reading and Artifacts below; provenance-only notes may
omit evidence sections. Store in memory efforts or link from the attempt.
Production evidence consumes analytic closure and a frozen remainder/design;
exploratory_only samples retain hypothesis/debugging scope. Name the next step
as the construction or distinction it earns. Begin sections with prose for indexing.

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

## Analytic Construction and Objective Bridge

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

## Observations

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
