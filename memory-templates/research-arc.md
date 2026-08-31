# Physics Research Arc Template

Use with `$physics-erdos-loop`. The contract is the live plan, append-only research history, continuation state, and promotion record for one positive physics objective.

```md
---
description: <exact physics object, mechanism, or derivation sought>
author: <agent-id>
created: '<ISO-8601>'
updated: '<ISO-8601>'
tags:
- substrate-framework
- research-arc
- <topic>
category: efforts
confidence: exploratory
status: active
---

## Positive Objective and Success
State the requested result precisely and explain why it is scientifically useful.
Completion requires the result itself, accepted dependency closure, natural
framework fit, one claim-appropriate verifier audit, mechanism comparison only
when genuinely open, impact-bounded replay, importable implementation, accepted
promotion, synchronized records, and no unresolved defect inside the promoted
scope. Honest failure is attempt evidence; an explicit hypothesis is not debt.

## Authority and Prior Work
Record the accepted release, source commit, claim ids, importable modules, campaigns, memory searches, and graph/source queries actually checked. Distinguish accepted state, proposal state, and attempt state. Accepted canon controls release and promotion; it remains falsifiable, so record credible independent challenges rather than treating conflict as automatic rejection.

- Accepted release:
- Accepted claims reused:
- Source modules read:
- Memory searches:
- Campaign evidence:
- Genuine unresolved objective:

## Definitions, Invariants, and Analytic Specification
Freeze variables, domains, units, normalizations, symmetries, topology,
initial/boundary conditions, regularity assumptions, and cross-sector invariants
before candidate selection. After structural selection and before every planned
numerical obligation, work and record exact identities and elimination,
complete variations, conservation
and symmetry reductions, non-dimensionalization and scaling, analytic bounds,
virial identities, limits, asymptotics, perturbation theory, and applicable
existence/uniqueness/compactness/no-go theorems. State the strongest conclusion
already earned, one residual proposition that still needs computation, and why
the current analytic ladder does not decide it. Only then freeze precision,
chart or ansatz, boundary treatment, discretization family, domain truncation,
fit form, stability condition, error norm, thresholds, and maximum numerical
verdict before inspecting production results.

## Permitted Imports and Assumptions
List every permitted external or framework input with provenance. Distinguish a
declared hypothesis or interpretive assumption from a hidden premise. Only a
hidden or promised-but-undischarged input enters the debt ledger.

## Candidate Set
Register at least two plausible concepts only when selecting an open scientific
mechanism. For a fixed theorem or specified construction, register one complete
route and add another only when it materially reduces uncertainty.

| Candidate | Construction | New objects/parameters | Natural-fit case | Expected falsifier | Status |
| --- | --- | --- | --- | --- | --- |
| A |  |  |  |  | untested |
| B |  |  |  |  | untested |

## Selection Criteria and Comparator Gate
Freeze the ordered structural criteria before seeing comparison values: invariant compatibility, explanatory reach, assumption cost, parameter economy, symmetry/topology, dimensions, limits, cross-sector composition, robustness, then empirical prediction. State when comparator values may be opened and any justified exception.

## Claim Delta
List each claim proposed, challenged, qualified, or potentially superseded. A proposal uses `challenges`; only accepted replacements use `supersedes`.

| Claim id | Exact statement | Dependencies | Relationship | Oracle | Consumers |
| --- | --- | --- | --- | --- | --- |
|  |  |  |  |  |  |

## Frozen Review Transaction
Name the exact claim delta, changed implementation, new or changed evidence
groups, declared dependency propositions, affected consumers, base/head or tree
hash, and planned validation receipt. Unchanged accepted dependencies, adjacent
corpus items, and raw theorem counts are outside the transaction.

## Claim Ladder
Build the smallest dependency-first ladder. Each row names the strongest
analytic result available, any irreducible numerical remainder, the strongest
practical oracle, and its sensitivity test. Use SymPy or manual algebra for
exact obligations, Lean for an auditable formal theorem, and SciPy only for a
named root/spectrum/integral/optimization/ODE/BVP/PDE remainder after the
analytic-closure receipt passes. Absence of a convenient closed form is not by
itself sufficient, and no obligation is forced into one tool.

| Step | Positive intent | Requires | Pass licenses | Does not license | Analytic result / numerical remainder | Oracle | Failure scope | Unlocks | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 |  |  |  |  |  |  |  |  | pending |

## Importable Implementation
Name canonical package APIs to add or reuse. Campaign scripts must call these APIs and must not duplicate constants, solvers, profiles, convention conversions, or check helpers. For numerical work, state whether `substrate_framework.numerics` applies and identify the claim-owned equation, operator, initial/boundary data, mesh, tolerances, and error metric. Route canonical sampled trapezoidal integration through its compatibility helper; mutable standalone scripts for the current environment use `np.trapezoid`, never removed `np.trapz`, and tractable exact integrals stay symbolic. Preflight direct attributes, imported names, and dynamic `getattr` access; nested fallback defaults are eager, so use a two-step `None` fallback.
Structure the verifier corpus as one small standalone module per claim with a thin aggregator; a monolithic multi-hundred-line script localizes failures poorly.

## Campaign-Branch Checkpoints
Record every locally complete unit as a coherent commit and append-only attempt
on the campaign branch. Do not open a PR for a rung, utility, subclaim, partial
goal, or clean milestone. These checkpoints preserve value without changing the
active obligation or terminal state.

| Unit | Local proposition | Evidence | Commit | Parent obligation advanced | Downstream claims still unlicensed |
| --- | --- | --- | --- | --- | --- |

## Attempts
Append one row per scientific attempt. Preserve source, stdout/stderr, elapsed time, and exact command. A native immutable-source abort caused only by missing direct, imported, or dynamic `np.trapz`—including an eagerly evaluated nested fallback—is compatibility provenance: run an alias-only compatibility replay and use that replay for scientific adjudication rather than consuming or rejecting a candidate. Failure of the repaired scientific route triggers the next route; it never closes the arc.

| Attempt | Obligation | Candidate/method | Artifact | Route verdict and scope | Licenses earned/missing | Method repair | Representation change | Alternative concept | Routes remaining |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0001 |  |  |  |  |  |  |  |  |  |

## Framework-Fit Audit
Assess each candidate before empirical fit: invariant preservation, imports, free parameters, conventions, known limits, cross-sector consumers, and whether it demands unrelated narrative changes. Determine whether a mismatch is a candidate defect or independent evidence against canon. Reject and replace a defective candidate; route a surviving canonical inconsistency through the separate revision gate without blocking correct conditional artifacts.

## Verifier Audit
Record the one strongest practical oracle and only the sensitivity evidence
needed for this claim: statement/axiom audit for Lean; applicable mutation,
counterexample, or exact limit for custom symbolic work; solver status, error
model, refinement, and one independent method or soluble limit for numerics.
For a numerical oracle, link the analytic-closure receipt and confirm that the
computed predicate is exactly its named remainder. Reuse an unchanged receipt.
Do not create a verifier for this record or expand adjacent observations into
the transaction.

## Impact-Bounded Dependency Replay
List only consumers that can change under the declared claim or implementation
delta and their replay command. Preflight relevant compatibility hazards before
assigning a scientific verdict. Do not replay unrelated sectors merely because
they share a corpus, registry, or historical campaign.

| Consumer | Why affected | Command or proof | Result | Repair if needed |
| --- | --- | --- | --- | --- |
|  |  |  |  |  |

## Foundational Revision Gate
Leave empty unless independent evidence shows accepted foundations are inconsistent. If opened, link a separate proposal containing the pre-existing inconsistency, two or more repair candidates, minimum-change decision, migration map, independent review, and full replay. This is an active advancement route when the evidence qualifies; never use it merely to rescue a favored concept.

## Debt Ledger
Track hidden assumptions, unresolved promises, unexplained fitted parameters,
broken affected consumers, or contradictory accepted narration created inside
the promoted scope. Declared hypotheses, honest exclusions, the still-open
parent objective, adjacent observations, and unselected routes are not debt.

| Debt | Source | Effect | Discharge | Status |
| --- | --- | --- | --- | --- |

## Independent Claim Review
Link one `claim-review.md` instance per claim proposed for acceptance or changed
accepted statement. Use `evidence-attachment-review.md` for changed attachment
roles; one record may cover a coherent set of entrypoints. Review once, then
perform one correction check limited to requested changes.

## Results and Continuation
Lead with the strongest meaningful positive result and its reproduction command.
Then record checkpointed progress, the active obligation, complete route
frontier, and findings as
`blocking in-boundary`, `minimum correction`, or `follow-up`. Reserve `refuted`
for an explicit contradiction or counterexample; missing support remains
unverified or qualified. If one route becomes repetitive, stop repeating that
route and continue through a method repair, representation change, or materially
different candidate. Do not hand off or open a PR at the milestone.

## Promotion and Materialization
Record extracted APIs/tests, accepted registry entries, immutable campaign location, release id, generated docs command, accepted-memory synchronization, and proposal/attempt memory separation.

## Scientific Exhaustion Certificate
Leave incomplete while any plausible route remains. Define the candidate
universe; inventory historical, external, preregistered, and failure-generated
routes; partition equivalent variants; attach each route verdict and
continuation-ladder record; include an independent adversarial generation pass;
and provide coverage/no-go evidence for infinite classes. Exhaustion requires
an empty route frontier and review of this certificate.

| Candidate class | Tried routes | Equivalent variants | Failure-derived routes | Coverage evidence | Verdict |
| --- | --- | --- | --- | --- | --- |

- Independent candidate-generation artifact:
- Routes remaining:
- Exhaustion review:

## Done Gate
Check applicable success conditions from `AGENTS.md` once at the unchanged
promotion boundary. Record one receipt containing the base/head or tree hash,
impact surface, exact command, and result. Do not rerun it for review prose,
counts, evidence-role corrections, or generated summaries that cannot affect
the oracle. Open one campaign PR only after the complete positive contract or
the scientific-exhaustion certificate in `AGENTS.md` passes.

## Cross-References
Link the proposal, source artifacts, claims, reviewers, campaign, release, generated outputs, and related memory entries.
```
