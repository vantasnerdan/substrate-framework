# Campaign Proposal Template

Instantiate before a campaign computes or inspects comparator values. Use `theorem-synthesis.md` instead when the target is a fixed higher theorem composed from accepted claims. Store the prose contract in memory and create a matching `proposals/<id>/proposal.yaml` manifest. Run `PYTHONPATH=src .venv/bin/python scripts/validate_repository.py`, validate repository-local memory with `memory validate --base "$PWD" "$PWD/memory"`, and preserve any schema or memory failure before opening the source body or comparator values; a prose contract alone is not the freeze gate. The memory path is a required positional target; `--base` alone does not select it.

Begin every section with a plain-prose sentence. Inline code, a table, or a
list does not satisfy the memory index's first-content disclosure contract.

```md
---
description: <positive campaign objective>
author: <agent-id>
created: '<ISO-8601>'
updated: '<ISO-8601>'
tags:
- substrate-framework
- campaign-proposal
category: proposals
confidence: exploratory
status: active
---

## Question and Positive Deliverable
State the exact question and object to derive. A no-go, failed concept, residual, or honest account of an obstruction does not complete this campaign.

## Obligation Graph and Closure Map
Decompose the parent objective into dependency obligations, not milestone-sized
deliverables. For each node record its positive intent, prerequisites,
`pass_licenses`, `does_not_license`, maximum verdict, failure scope, and the next
nodes it unlocks. A successful node activates its next unsatisfied dependency; a
failed route leaves the node active. No node, subclaim, utility, or clean result
opens a PR or terminates the campaign.

| Node | Positive intent | Requires | Pass licenses | Does not license | Maximum verdict | Failure scope | Unlocks | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |

## Base Release and Provenance
Record the accepted release and commit. Resolve inventory paths against the pinned source root, record both locations, and verify the source checkout commit and file hash before execution; a queue path need not be relative to the framework working directory. List source claims and modules actually read. For predecessor work, name each hash-pinned `migration/source-claims.yaml` unit and its current disposition; its bridge is the candidate unit while linked dossiers, formalizations, and legacy rungs are evidence rather than extra claims. Newer directories and working-tree prose are not authority.

## Source Inventory and Access Gate
For campaigns whose objective is external literature (a paper, a theory, a dataset): enumerate every load-bearing external source with verified access status BEFORE preregistration — in hand (local path), open (URL), paywalled, or missing — and the exact claims and page/equation numbers extracted from each. An inaccessible primary source blocks the campaign objective as written: escalate to the requester with options (supply the document, restate the objective against the accessible corpus, substitute) before any computation. Auditing secondary literature about an unchecked primary is a defect (skill `quantitative-verification`, AP-14), not a partial result.

| Source | Access status | Extracted claims (with page/eq) |
| --- | --- | --- |

## Invariants, Conventions, and Allowed Imports
Freeze what the campaign must preserve and every input it may use. Accepted canon governs the base release but remains challengeable. Record evidence that would distinguish a candidate defect from an independent canonical inconsistency; anything added later becomes explicit debt and requires proposal revision.

## Frozen Candidate Universe
Derive the in-scope route families from the user's original objective,
invariants, source inventory, and permitted imports before attempts begin. This
boundary expands append-only when historical, external, or failure-generated
concepts are discovered. It may not shrink during an exhaustion audit without
explicit user approval.

| Route family | Why in scope | Known concepts | Coverage strategy |
| --- | --- | --- | --- |

## Candidate Preregistration
Register at least two plausible concepts when selecting among scientific mechanisms unless uniqueness is proved. If the statement is a fixed theorem, declare that target kind and use one complete proof route rather than inventing a competitor. Do not retrofit the framework after selecting a mechanism.

| Candidate | Description | Assumptions | Parameters | Framework-fit prediction | Decisive test |
| --- | --- | --- | --- | --- | --- |
| A |  |  |  |  |  |
| B |  |  |  |  |  |

## Selection Criteria and Blinding
Order the structural selection criteria and state the comparator-blinding point. Numerical agreement is evaluated only after structural selection freezes.

## Proposed Claim Delta
List claims proposed or challenged, their dependencies, evidence plan, and consumers. Before assigning a claim identifier, search the registry, campaigns, and durable memory; rejected or provisional identifiers remain reserved even when absent from the accepted registry. Do not use `supersedes` before acceptance.

## Analytic Specification and Numerical Licenses
For every obligation separately, freeze `object -> symmetry/conservation ->
ensemble and exact variational functional -> admissible function space and
representation -> analytic scale/asymptotic structure -> observable ->
irreducible numerical remainder -> numerical approximation -> permitted
verdict`. Record the configuration space and all admissible variations;
equivalence or gauge, topology and bundle charts; group-action domain/period and
generator normalization; invariance and conservation identity; complete
Euler--Lagrange, Legendre, and second-variation objects required by the claim;
dimensions, non-dimensional groups, dominant balances and scaling laws;
available identities, bounds, coercivity/monotonicity results, virial identities,
limits, perturbative reductions, asymptotic principal symbol and essential
spectrum; and the strongest conclusion already established by algebra, calculus,
or an applicable theorem.

Only after that analytic-closure receipt passes may the obligation name one
residual proposition for production numerics. State why the current analytic
ladder does not decide it, without pretending to prove that no unknown analytic
method exists; then freeze background residual and forward error, branch
identity, tangent space and gauge, kinetic-metric rank/sign, representation
coverage and excluded sectors, observable and contamination exclusions, imposed
support/mask/taper, chart or ansatz, box and boundary treatment, discretization
family, fit form, tolerance, and maximum verdict. Give each license a durable
identifier and evidence-backed status. A downstream `requires` entry resolves
to a license identifier; an absent or unearned license blocks production
numerics rather than becoming a negative scientific result. Pre-gate numerical
sampling is `exploratory_only`: it may generate hypotheses or debug code but may
not select the candidate, set gates, support a claim or anti-claim, or later be
relabelled as production evidence.

## Implementation and Oracle Plan
Name importable APIs, claim-appropriate exact/numeric/formal oracles, applicable mutations and counterexamples, refinements, independent routes, and impact-bounded replay commands. A kernel-checked Lean proof calls for statement, import, proof-escape, axiom-footprint, and physical-encoding audits rather than a ceremonial mutation of the kernel. Campaign verifiers run directly with `PYTHONPATH=src`; import reusable package APIs rather than repository scripts, which remain CLI adapters. Pin the campaign's own source, claim, and release evidence, but never make future valid work fail by asserting unrelated queue units stay pending or mutable `current` remains the historical release; replay old campaigns through durable snapshots or their canonical modules/tests. For each replay inventory, record lexical check-call sites, runtime check executions, and assertion nodes separately; loops and dynamic dispatch can make the runtime tally differ legitimately, so equality is not an oracle. Predeclare a compatibility preflight: canonical integration uses `trapezoid_integral`, mutable current-environment scripts use `np.trapezoid`, and executable syntax is checked for direct, imported, and dynamic legacy access. An eager fallback such as `getattr(np, "trapezoid", getattr(np, "trapz"))` is legacy access because the default is evaluated first. Repair mutable code to the current name or a safe two-step fallback; give immutable source an alias-only recorded replay before scientific adjudication. Do not count that native compatibility abort as candidate rejection. State why SymPy, Lean, or a particular SciPy method fits each obligation. Do not plan a numerical rerun as independent evidence when an exact result already fixes its right-hand side or output; classify it as regression coverage and prefer exact sensitivity or Taylor separation for tractable counterexamples. Before labeling a downstream tail, dispersion, normalization, or consistency route independent, eliminate shared intermediate variables and compare the resulting equations or positive solution sets. For cross-sector matches, freeze distinct field types, kinetic metrics, action measures, and coefficient conversions; equal symbols, shapes, or dimensions are not maps. Structural oracles must evaluate the claimed object rather than a literal boolean, stand-in constant, copied period, or unrelated bounded sample. For differential forms, predeclare the full graded Leibniz/cyclic expansion and keep nonvanishing, closedness, global non-exactness, period normalization, filling dependence, and gauge descent as separate gates. For genuinely unresolved ODE/BVP/PDE or quadrature work, specify precision, equations, domain, initial/boundary data, discretization, mesh/time/sample refinement, tolerances, error norm, invariants or controlled dissipation, solver-status gate, and method cross-check. For FFT differentiation or spectral line claims, freeze the active frequencies and window, require commensurability or measured endpoint closure, distinguish an identity on one FFT coefficient from independent evidence, and predeclare the claimed line's minimum norm or power fraction. Express near-zero and agreement thresholds in a declared dimensional or scale-relative error model, and keep exact analytic nulls separate from numerical roundoff regressions.
Authoring practices from reviewed campaign defects: before implementing
any numerical scheme, search installed skills for its regime (soft modes,
small ratios, stiff-plus-soft optimization) and apply those prescriptions;
write units, signs, and geometric factors as a comment block before code
and compute test expectations independently of the implementation; capture
verifier stdout into `attempts/000N/` on first execution rather than
rerunning completed runs to materialize records; and when a symbolic check
repeatedly misbehaves, preserve the attempt and rederive it from source with a
different exact identity, manual calculation, independent CAS route, or formal
encoding. Numerical samples may debug that derivation but cannot settle or
refute the exact proposition.
When patching existing files, anchor by content rather than remembered
line numbers: search for a unique pattern at the edit site first, prefer
AST-aware rewrites for nested or multi-site changes, and switch to one
full rewrite from a full read after repeated failed patch rounds — stale
anchors were the largest single source of P242 rework.

## Attempts and Continuation
Append every route with its active obligation, verdict, diagnosed layer,
licenses earned/missing, method repair, representation change,
failure-generated alternative, and routes considered/tried/remaining. Stop
repeating a dead route but continue the obligation. An ill-fitting concept is
rejected or reformulated; unrelated earlier work is not rewritten to save it.
If a conflict survives independently of the candidate, open a separate
`challenges` or foundational-revision proposal rather than treating canon as
irrevisable or silently changing it. Execution runs in declared waves whose
inputs are explicit; research/grounding is wave 0 and a dependent wave
(implementation, verification, report) opens only after every input its wave
declares has settled—a research subagent's output is a hard prerequisite, never
a race (AP-15).

Every symbolic or numerical result also records `computed_predicate`, the exact
`proposition` it implies, licenses consumed and earned, `parent_effect`,
`cannot_decide`, and its maximum verdict. Agreement on an unlicensed object or
background cannot supply the missing bridge.

| Attempt | Obligation | Route | Verdict and failure scope | Licenses earned/missing | Method repair | Representation change | Alternative concept | Routes remaining |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |

## Debt Ledger
Track hidden assumptions, unexplained fitted parameters, unsupported promises,
convention conflicts, and broken affected consumers inside the proposed claim.
Declared hypotheses, honest exclusions, open candidate routes, and adjacent
repository observations are frontier rather than debt.

## Review and Promotion Plan
Name one claim-level review per proposed or changed claim, the frozen transaction,
package extraction, release update, generated documentation, and accepted-memory
synchronization. Classify changed attachments by evidence role; do not turn them
into additional claim reviews. Preserve the strongest meaningful positive
statement and use the minimum correction before rejecting it. Reserve
`refuted` for a contradiction or counterexample. For predecessor migration,
edit `migration/dispositions.yaml` and regenerate `migration/source-claims.yaml`;
never hand-edit the queue. Record one impact-selected validation receipt and
reuse it. After corrections, check only changed statements and affected edges.

## Scientific Exhaustion Certificate
Leave this section active and incomplete unless positive success passes first.
To claim exhaustion, define the full in-scope candidate universe; inventory
preregistered, historical, external, and failure-generated concepts; partition
equivalent variants; give every route a verdict and continuation-ladder record;
attach the independent adversarial candidate-generation result; supply coverage
or no-go arguments for infinite classes; and show `routes_remaining: []`. A
failed optimizer, soft signal, nonconvergence, elapsed effort, or finite route
count cannot satisfy this section.

| Candidate class | Routes and equivalence partition | Method repair | Representation change | Alternative concepts | Coverage evidence | Terminal verdict |
| --- | --- | --- | --- | --- | --- | --- |

- Independent candidate-generation artifact:
- Infinite-class coverage/no-go artifacts:
- Routes remaining:
- Exhaustion review:

## Done Gate
The campaign PR gate opens only on complete positive success or the scientific
exhaustion certificate in `AGENTS.md`. Until then, checkpoint with commits and
continue; do not open a rung, milestone, subclaim, utility, or progress PR. Each
next step states its positive contribution to resolution: name the object it
constructs, question it closes, or distinction it establishes. A refutation
counts only after the object, ensemble, observable, admissibility, and
representation-coverage prerequisites pass, and only at its declared failure
scope. Avoidance and risk reduction do not describe a next route.
```

Matching manifest:

```yaml
id: P000
schema_version: 2
base_release: <release-id>
source_baseline: <repository and immutable commit or release>
question: <exact question>
candidate_universe:
  scope: <original user objective and in-scope class>
  frozen_from: [<objective>, <invariants>, <source inventory>, <permitted imports>]
  route_families: [<family A>, <family B>]
  append_only_expansions: []
obligation_graph:
  nodes:
    - id: <node-id>
      positive_intent: <object or proposition that must become established>
      requires: []
      pass_licenses: []
      does_not_license: []
      maximum_verdict: <typed verdict>
      failure_scope: <candidate class touched by failure>
      unlocks: []
      status: <pending | active | established | exhausted | refuted>
      license_chain:
        object: <exact mathematical/physical object>
        symmetry_or_conservation: <identity, license id, or not applicable>
        ensemble: <fixed charge | fixed frequency | statics | mode | dynamics>
        variational_functional: <complete functional varied>
        admissible_space: <smoothness, constraints, gauge, topology, boundaries>
        representation_coverage: <included and excluded sectors>
        observable: <definition and contamination exclusions>
        numerical_representation: <blocked until analytic closure; then background, branch, operator, error budget>
        permitted_verdict: <maximum typed verdict>
license_registry:
  - id: <license-id>
    proposition: <exact proposition supplied>
    status: unearned
    evidence: null
    earned_by: null
source_inventory:
  - source: <paper/dataset, author title year>
    access: <in-hand path | open URL | paywalled | missing>
    extracted: <claims/pages actually read, or escalation note>
invariants:
  - <accepted invariant>
allowed_imports:
  - <claim id or external source>
candidates:
  - id: A
    description: <candidate concept>
  - id: B
    description: <candidate concept>
uniqueness_evidence: null
selection_criteria:
  - framework invariant compatibility
  - assumption and parameter economy
  - correct limits and cross-sector composition
claims_proposed:
  - <claim id>
comparators_blinded_until: <artifact or gate>
route_frontier:
  active_obligation: <node-id>
  considered: []
  tried: []
  failure_generated: []
  remaining: []
execution_state: active
objective_state: active
exhaustion_certificate:
  historical_routes: []
  external_routes: []
  failure_generated_routes: []
  equivalence_partition: []
  adversarial_generation_artifact: null
  adversarial_reviewer: null
  adversarial_reviewer_role: null
  infinite_class_coverage: []
  route_verdicts:
    - route: <exact route id from route_frontier.considered>
      verdict: <established | refuted | blocked, with route scope>
      evidence: <attempt artifact>
      continuation: <method repair, representation change, and alternative generation record>
  routes_remaining: null
  review: null
status: draft
```
