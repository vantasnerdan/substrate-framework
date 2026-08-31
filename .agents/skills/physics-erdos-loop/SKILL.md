---
name: physics-erdos-loop
description: Run persistent, verifier-backed physics research and framework reconciliation from candidate generation through claim-level promotion. Use for physics derivations, equations, ODE/PDE work, symbolic or numeric checks, simulations, Lean formalization, campaign design, claim migration, framework-wide consistency audits, or any proposal that might change accepted scientific claims. Enforces typed problem deconstruction, campaign-terminal PRs, exhaustive route continuation, natural framework fit, append-only attempts, claim-appropriate verification, impact-bounded replay, and generated canonical records.
---

# Physics Erdős Loop

Use this loop to produce a positive, framework-consistent, verified result—not merely an articulate account of why one attempt failed. It adapts the Erdős-style persistent proof loop to physics while adding claim governance and software reuse.

## Non-negotiable outcome

Honesty is required but is not itself success. Classify a failed candidate, no-go, obstruction, residual, bound, or inconclusive computation as attempt evidence and continue. Do not close the effort on it. Equally, do not turn honest review into a ratchet toward trivial claims: preserve the strongest meaningful positive result, make the minimum truthful scope repair, and leave the larger objective open when it is not yet met.
Weakening is subject to the same persistence discipline as failure-chasing. Each
attempt must change something material — candidate concept, numerical method, oracle
class, formalism, or target statement — rather than re-examine a known route under
tighter standards; repeating a route with only stricter evidence thresholds is drift,
not diligence. A new acceptance criterion invented mid-campaign must name the observed
failure it prevents, be recorded once, and bind future attempts rather than
retroactively invalidating finished ones. Close each attempt with two typed
fields: `route_verdict` is established, refuted with mechanism, or blocked with
that route's missing construction; `evidence_scope` states the strongest earned
scientific meaning, such as `NUMERICALLY_UNRESOLVED` or
`REPRESENTATION_SCOPED`. Neither field propagates automatically to the
obligation or campaign. The next attempt inherits a direction, not a smaller
question.

Success requires the requested object plus all of these gates:

- accepted dependency closure and declared imports;
- natural compatibility with framework invariants, or a separately accepted minimal foundational revision;
- strongest practical oracle and claim-appropriate sensitivity evidence;
- preregistered comparison of plausible candidates when a mechanism is being selected;
- impact-bounded downstream replay;
- reusable importable implementation;
- claim-by-claim review and release promotion;
- synchronized generated docs and memory;
- empty debt ledger.

Read [governance.md](references/governance.md) before changing a claim, convention, invariant, or canonical API. Read [oracles.md](references/oracles.md) when selecting or auditing verification.
Read [problem-deconstruction.md](references/problem-deconstruction.md) before
freezing a multi-rung campaign or any numerical representation of a physical
object, ensemble, localization, mode, stability, force, or interaction claim.

## Supporting skills

Load these at the moments they change decisions, not after results:

- `small-ratio-numerics/SKILL.md` — required before freezing the verifier
  design of any gate whose quantity is a soft Hessian eigenvalue, a
  stability-window edge, a Morse index, a force, or an energy difference
  within roughly three orders of magnitude of the discretization,
  quadrature, or roundoff floor — including reproduced gates and gates that
  reuse committed machinery, where inherited floors otherwise go unmeasured.
  The attempt manifest records which prescriptions bind (error budget,
  zero-mode gauge, eigenpair residuals, lambda_min/lambda_2, observed-order
  extrapolation, jitter sign test) and how each was satisfied: by
  measurement, a stated bound, or explicit scoping of the claim below the
  floor.
- `theorem-synthesis/SKILL.md` — load when composing accepted claims into a
  higher theorem; reground on it periodically during long synthesis
  campaigns to keep the composition logic sharp.

## Core distinctions

- A campaign is an immutable research event, not canonical truth.
- A proposal can challenge a claim; only an accepted claim can supersede one.
- Accepted canon governs release and promotion decisions, but remains falsifiable and reviewable; a conflict is a diagnosis to investigate, not a reason to erase evidence or halt conditional artifact work.
- A verifier passing is necessary, but proves only its asserted predicate.
- A support gap is not a refutation. Reserve “false” for a contradiction or
  counterexample under the stated hypotheses.
- Evidence can be an exact proof, corroborating subclaim, regression,
  applicability test, or provenance record. Only the first is expected to prove
  the whole encoded statement; label the others instead of discarding them.
- Exact algebra, calculus, scaling, asymptotics, and applicable theorems define
  the claim before production numerics. A numerical verifier starts from the
  named remainder of that analytic specification; it does not supply missing
  semantics or choose its own continuum question.
- Exact proof and empirical applicability are separate obligations: measurement can test nature without becoming the proof of a symbolic or formal implication.
- A campaign completion decision, scientific claim decision, and PR merge decision are independent. An active campaign nevertheless opens no PR until positive completion or certified exhaustion; use `research-pr-harvest` only for an externally supplied PR or after that terminal gate opens.
- Numeric agreement is a comparator, never a concept-selection mechanism or hidden derivation input.
- A failed concept is evidence about that concept. It is not permission to retrofit the framework around it.
- An attempt is bounded and append-only. The effort continues until the success contract is met or the user changes the objective.

### License the object before measuring it

Every computational claim follows the typed chain `mathematical object ->
symmetry or conservation license -> ensemble and exact variational functional ->
admissible function space and representation -> analytic scale/asymptotic
structure -> observable -> irreducible numerical remainder -> numerical
approximation -> permitted verdict`. Before a downstream node can earn
scientific evidence, its manifest records the upstream licenses it consumes and
an analytic-closure receipt: exact equations and variations, symmetries and
constraints, scaling and non-dimensional groups, analytic bounds and limits,
asymptotic operator or continuum threshold, strongest non-numerical conclusion,
and the one proposition left for computation. A computable tangent, finite
quadratic norm, fitted frequency, sampled spectrum, or stable discretization
does not create a symmetry, conserved charge, cyclic coordinate, physical
ensemble, observable, or continuum statement.

For each obligation node record `requires`, `pass_licenses`,
`does_not_license`, `maximum_verdict`, `failure_scope`, and `unlocks`. Passing or
failing a node changes only the propositions named there. After success, advance
to the next unsatisfied parent obligation. After route failure, keep the node
active and continue through method repair, representation change, and
materially different candidate generation. No fixed number of routes proves
exhaustion.

## Phase 0 — establish authority and recall

Before deriving anything:

1. Read `AGENTS.md`, `governance/releases/current.yaml`, and `governance/claims.yaml`.
2. Inspect git status and history. Distinguish accepted release, committed provenance, uncommitted proposal work, and generated files.
3. Search repository memory with `memory search ... --base "$PWD/memory"` and `memory grep ... --base "$PWD/memory"`; treat hits as pointers and verify facts at source. Validate repo-local memory with an explicit repository base and absolute target, `memory validate --base "$PWD" "$PWD/memory"`, so host configuration cannot redirect relative paths.
4. Search importable source, campaign artifacts, tests, and dependency consumers. For predecessor migration, start from the hash-pinned unit in `migration/source-claims.yaml`; do not double-count its dossiers, frozen rungs, formalizations, or memory entries as independent claims.
5. Record the last accepted boundary and the genuine unresolved objective in a contract from `memory-templates/`.

Run `.agents/skills/physics-erdos-loop/scripts/preflight.sh` to check the local tools and governance surfaces.

## Phase 1 — write the success contract

Instantiate `memory-templates/research-arc.md` for physics work or `campaign-proposal.md` for a campaign. State:

- the exact positive deliverable;
- the parent-to-obligation closure graph and the next dependency each node unlocks;
- the accepted base release;
- definitions, variables, units, domains, quantifiers, and conventions;
- invariants that must survive;
- permitted imports and assumptions;
- the claim delta and downstream consumers;
- each proposed claim identifier and a repository-wide registry, campaign, and durable-memory collision search, because rejected provisional identifiers remain reserved;
- what each oracle must establish;
- the analytic-closure receipt for every planned numerical obligation, including
  the unresolved proposition, why the current analytic ladder does not decide
  it, frozen numerical design freedoms, and maximum numerical verdict;
- the empty-debt and canonicalization gates.

Do not include failure, no-go, residual, or “best effort” as an accepted outcome.

Validate the matching proposal manifest with `PYTHONPATH=src .venv/bin/python scripts/validate_repository.py` before opening a predecessor source body or comparator values. The frozen prose and YAML must agree, and a schema failure is an append-only attempt rather than permission to proceed informally.

## Phase 2 — preregister the actual choice

Register at least two plausible candidates when the scientific mechanism is genuinely being selected, unless uniqueness is already proved. For a fixed theorem statement, declare `target_kind: fixed_theorem` and register one complete proof route; compare proof routes only when a second route materially reduces uncertainty. Never invent a rival mechanism to satisfy a form. For each applicable candidate or proof route, record its new objects, assumptions, parameters, expected limits, affected claims, and likely consumers.

Freeze selection criteria before inspecting comparison values:

1. consistency with accepted invariants;
2. explanatory and predictive reach;
3. fewer new assumptions, imports, and parameters;
4. correct symmetries, dimensions, topology, and limits;
5. compatibility with other accepted sectors;
6. numerical robustness and implementability.

Keep empirical comparators blinded until equations, conventions, criteria, and structural tests are frozen when practical. Record and justify any exception.

## Phase 3 — close the analytic specification, then derive through importable APIs

Build the smallest dependency-first claim ladder. Implement reusable equations, constants, units, solvers, and transformations under `src/substrate_framework/`. Keep proposal scripts thin: import canonical functions and evaluate a candidate.

Before implementing a production numerical verifier, work the analytic ladder as
far as the claim permits: exact identities and elimination; complete first and
second variation; conservation and symmetry reduction; non-dimensionalization
and scaling; bounds, coercivity, monotonicity and convexity; virial identities;
limits, asymptotics and perturbation theory; and applicable existence,
uniqueness, compactness, or no-go theorems. Record the strongest result obtained
and one irreducible remainder. Absence of a convenient closed form is not by
itself a numerical-necessity argument, and the receipt need not prove that no
unknown analytic method could ever exist.

Match implementation to that frozen claim. Use exact symbolic algebra for
tractable identities and reductions, formal proof when the encoded theorem is
the real obligation, and SciPy only for the named root, spectrum, integral,
optimization, ODE, BVP, or PDE remainder. Exploratory numerical sampling before
the receipt is labeled hypothesis-generation or debugging; it cannot choose a
candidate, set gates, establish or refute an obligation, enter claim review as
scientific evidence, or later be relabeled as production evidence. Reuse
`substrate_framework.numerics` for common IVP, method-of-lines, BVP, refinement
evidence, and sampled trapezoidal integration; canonical modules call its
`trapezoid_integral` compatibility API, while mutable standalone scripts
targeting the current environment call `np.trapezoid`, never removed
`np.trapz`. Preflight executable syntax for direct, imported, and dynamic legacy
access. Never use an eager nested fallback such as `getattr(np, "trapezoid",
getattr(np, "trapz"))`; use the canonical helper or a two-step `None` fallback.
The claim implementation must still own and expose its equation, discretization,
boundary data, error norm, and physical acceptance thresholds.

Do not:

- execute simulations at import time;
- duplicate `check()` helpers, profiles, constants, or convention conversions;
- encode the expected answer as an input;
- reinterpret earlier variables to make a new concept fit;
- edit generated documentation;
- patch existing files from remembered line numbers — re-anchor by content
  search at the edit site, prefer AST-aware rewrites (ast-grep) for nested
  or multi-site changes and unique-pattern substitution for single lines,
  and switch to one full rewrite from a full read after repeated failed
  patch rounds; stale anchors were the largest source of P242 rework;
- implement a numerical scheme without first searching installed skills for
  its regime — soft modes, small ratios, stiff-plus-soft optimization are
  covered by `small-ratio-numerics`, whose prescriptions (non-dimensionalize
  the dominant balance, second-variation operators plus pseudoinverse rather
  than gradient methods on soft manifolds, mesh extrapolation, pinned BLAS
  threads) preempt whole debug cycles;

Use [verify_claim.py](assets/verify_claim.py) for exact or general claim checks and [verify_pde.py](assets/verify_pde.py) for a SciPy method-of-lines and mesh-refinement pattern. Both use shared framework APIs rather than campaign-local solver or tally copies.

## Phase 4 — run append-only attempts

Create `attempts/0001/`, `0002/`, and so on. Preserve candidate source, command, environment, stdout, stderr, elapsed time, verdict, and diagnosed mechanism. Never overwrite an attempt. Capture stdout into the attempt directory on first execution — records generated by rerunning already-completed verifiers re-spend exactly the compute the records exist to justify.

Each attempt records the active obligation, route verdict, failure scope,
licenses earned and still missing, routes considered/tried/remaining, and the
next materially different route. A numerical route may earn `refuted` only when
its object, ensemble, observable, admissibility, and representation coverage are
licensed and a converged error enclosure excludes the candidate's success set.
Nonconvergence, branch loss, continuum drift, unresolved state error, singular
normalization, imposed boundary support, or an incomplete perturbation space
sets `evidence_scope: NUMERICALLY_UNRESOLVED` or
`evidence_scope: REPRESENTATION_SCOPED`; it does not propagate a scientific
refutation to the obligation.

After failure, choose the next action from the diagnosis:

- legacy-library alias only → detect direct `np.trapz`, imported `trapz`, dynamic `getattr(np, "trapz")`, and eager nested-default access; repair the mutable script to `np.trapezoid` or a safe two-step fallback and rerun the same scientific route; for immutable hash-pinned source, preserve the native abort and run an explicit alias-only compatibility replay, without counting the environment abort as a rejected scientific candidate;
- implementation defect → repair and rerun;
- unstable numerics → change discretization, solver, precision, or oracle;
- bad representation → change variables, gauge, basis, coordinates, or formalism;
- concept conflicts with the framework → determine whether the defect belongs to the candidate or is independent evidence against accepted structure; reject or reformulate a defective candidate, otherwise open a `challenges` or foundational-revision proposal and keep the frontier active;
- target was misstated → correct the claim while preserving the user's objective;
- accepted foundation appears inconsistent → open a separate foundational-revision proposal.

Before inventing a replacement, inspect the nearest accepted campaign or
canonical module with a related obligation. Verify its source, extract the
construction and selection logic that actually worked, and record the
transferable assumptions and the present mismatch. External research may
supplement this step but does not replace framework-context reconciliation.

Persistence applies to the scientific objective, not repeated work on a dead
route. Commit clean milestones on the campaign branch and continue immediately;
do not open, land, or hand off a rung-level PR. Stop repeating a route when its
mechanism is known, then change method, representation, or concept until the
obligation is satisfied or the complete solution space meets the exhaustion
contract in `AGENTS.md`.

## Phase 5 — audit the verifier once at the frozen boundary

Choose the strongest practical oracle using [oracles.md](references/oracles.md), then audit it against the named claim:

Audit the analytic-closure receipt and objective bridge before auditing
predicates. Record the computed
predicate, the mathematical proposition it implies, every upstream license, the
maximum verdict it can earn, the parent obligations it advances, and the
questions it cannot decide. A correctly computed predicate with an absent
bridge or an unclosed analytic specification is exploratory or provenance
evidence only.

The audit asks whether the verifier establishes the named claim and where its
positive support ends. Do not create a second meta-verifier, add checks for
review prose or reviewer identity, or reopen unrelated accepted claims. A
discovered adjacent concern is a follow-up unless it directly falsifies the
proposed statement, removes a proposition used from a declared dependency, or
breaks an affected consumer. Record it once and return to the frozen boundary.

Do not count a weaker oracle as independent evidence when a stronger result already fixes its input. In particular, after exact algebra removes a parameter from an ODE right-hand side, local uniqueness proves same-initial-data trajectory independence; integrating that identical right-hand side twice is only regression coverage. Likewise, eliminate shared intermediate variables before calling a downstream tail, dispersion, or normalization check independent: if it yields the same equation or positive solution set, record it as a dependent regression. Cross-sector matching additionally requires explicit field, kinetic-metric, action-measure, and coefficient maps; equal names, shapes, or dimensions do not supply them. Use exact sensitivity or initial Taylor coefficients for analytically accessible counterexamples, and reserve simulation for behavior the exact result does not decide.

- confirm process status zero and terminal tally independently; record lexical check-call sites, runtime check executions, and assertion nodes as distinct inventories, without demanding equality when loops or dynamic dispatch multiply executions; with `CheckLedger`, preserve the formatted tally while returning its status-zero success token rather than a positive count as an OS exit code;
- pin a verifier's own source, claim, and release evidence, but do not freeze unrelated future queue dispositions or require mutable `current` to remain a historical release; historical replay uses durable snapshots when available and otherwise targets canonical modules/tests rather than rewriting an adjudicated campaign;
- mutate each load-bearing input and require a relevant check to fail for custom symbolic/numeric verifiers and translation layers. For a Lean theorem checked by the kernel, audit the exact statement, imports, proof escapes, axiom footprint, and physics encoding rather than performing a ceremonial mutation of the kernel;
- test wrong signs, normalizations, conventions, and counterexamples;
- check dimensions, symmetries, conserved quantities, and known limits;
- run resolution, timestep, domain, and tolerance refinement for numerics;
- compare against an independently implemented or analytically solvable case;
- inspect the exact statement of formal theorems and their axioms.
- reject structural predicates implemented as literal booleans, stand-in constants, copied periods, or samples that do not evaluate the defining object;
- for differential forms, expand all graded Leibniz and cyclic terms and audit nonvanishing, closedness, global non-exactness, periods, extension dependence, and gauge descent independently.

For SciPy work, record the routine and algorithm (`solve_ivp`, `solve_bvp`, sparse eigensolver, optimizer, quadrature, or another justified method), floating-point precision, mesh and domain, initial/boundary data, tolerances, stopping status, and error norm. Tie near-zero and agreement thresholds to a declared dimensional or scale-relative error model. When an absolute threshold fails at roundoff scale, preserve the attempt, show refinement or conditioning evidence, and repair the oracle with a justified scale-sensitive bound; do not blur a separately exact null into a floating-point claim. Treat solver success as a prerequisite, not the verdict. A PDE claim additionally needs spatial and temporal refinement, stability evidence, conservation or controlled-dissipation checks, and a method cross-check or soluble limit appropriate to the equation.

A large pass tally with insensitive predicates does not promote a claim.

## Phase 6 — assess framework fit before data fit

Compare candidates using the preregistered criteria. Structural fit precedes empirical closeness. If the favored candidate requires reinterpretation of unrelated claims, convention mixing, compensating imports, or narrative edits merely to preserve it, reject it and continue the search. If the mismatch is reproduced independently of that candidate, treat it as evidence about canon and route it through a separate challenge rather than assuming acceptance status resolves the science.

Do not revise foundations merely to save a candidate. A foundational revision is nevertheless a legitimate advancement route when it demonstrates an independent pre-existing inconsistency, compares at least two repairs, selects the minimum coherent change, enumerates the migration, and passes global replay. Until adjudication, the conflict blocks promotion and accepted downstream use, not truthful conditional APIs or continued investigation.

Only after the structural choice is frozen should you open the comparator gate and report predictive agreement or disagreement.

## Phase 7 — replay the dependency graph

Before review:

1. Enumerate direct and indirect consumers of every changed claim and canonical symbol.
2. Re-run targeted unit, symbolic, numeric, simulation, and formal checks. If a mutable consumer aborts on direct or dynamic access to removed `np.trapz`—including an eagerly evaluated nested `getattr` default—repair it to `np.trapezoid` or a safe two-step fallback and rerun before classifying the consumer or campaign; use an alias-only recorded replay for immutable source.
3. Re-check units, conventions, signs, limits, free-symbol sets, imported constants, and parameter counts.
4. Compare generated outputs and narrative consumers.
5. Record and discharge defects or hidden promises created inside the proposed
   claim delta. Keep unrelated observations and unfinished parent work on the
   follow-up frontier rather than in the promotion debt ledger.

Local success with broken downstream consumers is failure.

## Phase 8 — independent claim review

Use `memory-templates/claim-review.md`. Review each claim proposed for acceptance
or changed accepted statement once, not every sentence, theorem entrypoint, or
evidence attachment in the campaign. The reviewer must have the raw artifacts
and acceptance criteria, not the proposing agent's preferred conclusion.

Pin the review boundary to the proposed claim delta and its declared
dependencies. Review an evidence attachment at the exact scope it claims; do
not require it to rederive the parent claim or share every parent object when
it is explicitly labeled as a narrower corollary, regression, applicability
test, or provenance record. Use
`memory-templates/evidence-attachment-review.md` when attachment roles need an
audit; one record may cover a coherent group of entrypoints. After requested
corrections, perform one correction check limited to the changed statements and
directly affected dependency edges. Do not start a second substantive pass.

For a scope defect, report in this order: strongest meaningful supported
statement, exact unsupported extension, minimum repair, and evidence that would
restore the stronger version. Prefer a correct quantifier, explicit hypothesis,
layer split, or evidence-role change over deletion. Reject the whole statement
only when no useful positive core survives.
Scope reduction is legitimate only when it changes what downstream work can conclude;
reduction that leaves every consumer equally informed is drift, and the reviewer states
alongside any narrowing what evidence would restore the stronger form, so the next
agent inherits a path forward. When successive reviews of one claim have each reduced
its scope without new contrary evidence, stop narrowing: name the shrink history in
the review record and escalate to the user, because the review chain — not the claim —
is failing.

Assign each claim independent verification, review, compatibility, and epistemic statuses. Unaccepted work stays under `proposals/`. Use `challenges` until a replacement claim is accepted; only then add `supersedes`.

## Phase 9 — promote and materialize

For accepted claims:

1. Extract reusable logic into `src/substrate_framework/` and tests.
2. Update `governance/claims.yaml` and a pinned release manifest. Compute release closure from accepted registry membership (`review: accepted` / non-null `accepted_in`), not from `epistemic: active`; accepted qualified claims remain release members. Use the governance validator rather than a verifier-local status filter. If predecessor units were consumed, edit only `migration/dispositions.yaml`, preserve any unmigrated subclaims explicitly, and regenerate `migration/source-claims.yaml`; never hand-edit that generated queue. Materialize every evidence path before registering it. A final attempt that summarizes the promotion gate may begin with an explicit in-progress status, then be finalized after the gate and checked with only record-sensitive repository/generation validation.
   Terminal `qualified`, `refuted`, `duplicate_evidence`, and `out_of_scope` dispositions require their structured reason and durable evidence paths; use `qualified` when a mixed unit also maps accepted claims.
3. Move the adjudicated campaign record into the immutable `campaigns/` log.
4. Run `scripts/render_docs.py`; never hand-edit `docs/generated/`.
5. Generate or synchronize accepted claim/release memory. Keep proposal and attempt memory separate.
6. Run targeted scientific checks, the scope selected by `scripts/validate_changed.py`, and `git diff --check` once at the final frozen boundary. An append-only synthesized leaf theorem may remain scoped when it changes no existing claim or contract and its exact proof, registry/rendering, direct consumers, and formal surface replay; shared machinery, altered existing claims, foundational revisions, or uncertain impact require `scripts/validate.sh --full`. Record one content-addressed receipt and reuse it; review prose, evidence counts, and generated summaries do not stale a scientific gate. Run validation and commit in separate process invocations so an unguarded shell cannot continue past a failed gate and mask it with a later successful command.

## Phase 10 — done gate

Declare the campaign positively complete only when every applicable item in the
success contract passes and the in-boundary debt ledger is empty. Otherwise
continue until the scientific-exhaustion certificate in `AGENTS.md` passes. The
campaign branch is the checkpoint surface; no rung, subclaim, utility, partial
goal, diminishing-return point, or clean milestone opens the PR gate or ends the
run. After positive completion or certified exhaustion, freeze and submit one
terminal campaign PR. An exhaustion PR preserves the strongest supported
results without claiming the positive objective and uses `Advances`, not
`Fixes`, unless success also passed. A pause caused by user authority, runtime,
or an external dependency preserves the active contract and exact next
executable route; it is not scientific exhaustion. When a repeated workflow
defect or tooling gotcha is discovered, correct and consolidate the relevant
instruction in `AGENTS.md`, this skill, and the applicable memory template; do
not merely append another overlapping rule.

## Working with delegated agents

When delegation is authorized, give each worker one child contract and a disjoint write surface. Use one reviewer for one frozen claim transaction unless the promotion contract or user explicitly requires more. Fresh reviewers receive sourced inputs, the claim, and criteria—not the parent agent's interpretation or expected answer. They may report direct blockers and one compact follow-up list, but may not recursively commission broader audits. Reconcile blockers once and rerun only checks invalidated by the edits.
