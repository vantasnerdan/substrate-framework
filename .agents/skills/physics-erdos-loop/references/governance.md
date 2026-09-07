# Claim governance and framework fit

Use this reference whenever work may change a scientific claim, convention, invariant, dependency, release, or canonical module.

## Authority model

Authority flows from a pinned accepted release through accepted claims to their evidence. Campaign order records chronology only. Commit history establishes provenance only. Neither makes a claim true.

Every claim carries four independent status axes:

| Axis | States |
| --- | --- |
| Verification | `unverified`, `symbolic_verified`, `formal_verified`, `numeric_evidence`, `simulation_evidence` |
| Review | `unaudited`, `audited`, `accepted`, `rejected` |
| Compatibility | `unassessed`, `native`, `compatible_extension`, `conflict` |
| Epistemic | `proposed`, `active`, `qualified`, `superseded`, `refuted` |

Earn each axis independently: commits supply provenance, review supplies its own evidence, compatibility supplies framework fit, and derivation remains distinct from empirical agreement.

## Artifact roles

- `proposals/`: mutable candidate work before adjudication.
- `attempts/`: append-only records inside a proposal; never canonical.
- `campaigns/`: immutable adjudicated research events.
- `governance/claims.yaml`: machine-readable claim graph.
- `governance/releases/`: reproducible accepted claim sets.
- `src/substrate_framework/`: reusable implementation of accepted definitions and derivations.
- `docs/generated/`: canonical views rendered from accepted state.
- memory `efforts`/`proposals`/`attempts`: work state and recall pointers.
- memory `claims`/`releases`: generated or synchronized accepted-state summaries.

## Achievement — a sourced proposal

Before calculating, earn a frozen proposal by recording:

- base release and source commit;
- exact question and positive completion object;
- accepted invariants and conventions;
- permitted imports and assumptions;
- at least two candidate concepts when a scientific mechanism is being selected, unless uniqueness is proved; a fixed theorem target may register one complete proof route;
- selection criteria fixed before comparator values are used;
- proposed claim delta and anticipated consumers;
- comparator-blinding point;
- validation and impact-bounded replay plan.

Preserve the freeze through an explicit proposal revision when these fields
change. Replay only checks whose inputs or proposition changed.

For synthesis campaigns, set `campaign_type: synthesis` and
`target_kind: fixed_theorem`, name one higher claim and the structural gap, and
list at least two distinct accepted composition dependencies. The accepted
atoms are not re-reviewed. The new claim records `category: synthesized`, its
SymPy or Lean glue proof, assumptions, and exclusions. An optional
`layer: interpretive` theorem must state hypothesis H explicitly and cannot feed
the core layer.

## Achievement — natural framework fit

A candidate fits naturally when it:

- reuses accepted primitives without redefining them;
- preserves symmetries, topology, units, conventions, and known limits;
- reduces or leaves unchanged the assumption and parameter ledger;
- composes with other accepted sectors through explicit APIs;
- produces consequences without importing their desired values;
- requires no unrelated narrative rewrite to appear compatible.

Use a concrete mismatch to reformulate the candidate or generate another concept; the reward is a better fit, not rescuing the favorite.

## Foundational revisions

Earn a separate foundation proposal with evidence of a pre-existing inconsistency that does not depend on the proposed replacement. Compare at least two repairs, identify the smallest coherent change, list every affected claim and consumer, obtain independent review, and replay the whole dependency closure.

Keep this proposal separate from the favored candidate so independent evidence decides whether the foundation itself needs revision.

## Promotion transaction

Promote each proposed/changed claim individually while retaining unchanged
acceptance and accurate attachment roles:

1. Freeze the proposal and its attempt history.
2. Audit the exact claim and verifier sensitivity.
3. Assign all four statuses.
4. Validate accepted dependency closure.
5. Extract reusable implementation and tests.
6. Replay consumers.
7. Add accepted claim and graph edges.
8. Pin a release.
9. Generate docs and accepted memory.
10. Commit the complete transaction together.

Use one substantive review and one correction check. Evidence attachments carry
an explicit role—exact proof, corroborating subclaim, regression, applicability,
or provenance—in the review record and registry `scope` text; no new schema key
is implied. Related entrypoints may be grouped under one proposition.
Reserve `refuted` for an explicit contradiction or counterexample. When scope
needs repair, preserve the strongest meaningful positive statement through the
minimum quantifier, hypothesis, layer, or evidence-role correction.

Before acceptance, a proposal may record `challenges`. Only an accepted claim may record `supersedes`.

## Success and continuation

Bank a failed route as honest `failed_attempt` evidence, extract its reusable mechanism, and pursue a repaired or different candidate. The elite achievement is resolving the positive objective through better constructions; truthful contradictions remain visible. Honor user pauses under AGENTS.md.
