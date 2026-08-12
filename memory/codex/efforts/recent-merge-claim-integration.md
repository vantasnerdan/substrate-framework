---
description: Integrate P227 and P229 primitives with accepted gauge conventions and prepare exact claim-review packets
author: codex
created: '2026-08-12T11:27:58+02:00'
updated: '2026-08-12T11:36:38+02:00'
tags:
- substrate-framework
- effort
- claim-integration
category: efforts
confidence: working
status: active
---

## Goal and Success Contract
This effort delivers exact typed bridges from the recent P227 and P229 merge
artifacts into existing framework conventions, narrows their provisional claims
to the propositions actually established by their strongest oracles, and leaves
one complete review packet per claim for a distinct reviewer. It is complete as
an author handoff only when the bridge APIs and mutations pass, every affected
consumer replays, the proposal and review-request records agree, one final
repository validation passes, the debt ledger is empty, and a PR advancing
issue #59 is ready for a distinct reviewer. Scientific promotion remains
incomplete until that reviewer adjudicates each claim and an accepted registry,
release, generated-document, campaign, and memory transaction lands.

## Accepted Baseline
The accepted authority boundary is release `v0.159.0` at
`substrate@6d1f4e0`. The implementation base is
`substrate-framework@640c47f`. Accepted claims `C-LIE-001` and `C-RGE-005`
supply the canonical SU(3) generators and declared gauge-beta convention.
Proposal P227 supplies the unpromoted worldline and Lorentz atoms; proposal
P229 supplies the unpromoted chromomagnetic atoms. Merge history, proposal
memory, and the delivered P229 report establish provenance only.

## Constraints and Invariants
The work preserves P227 mostly-plus worldline conventions, positive einbein
and massive/massless branch separation, P229 background-Feynman-gauge and
`gB` conventions, the `T_a=lambda_a/2` SU(3) normalization of `C-LIE-001`,
and the `mu*dg/dmu=b*g^3/(16*pi^2)+...` beta convention of `C-RGE-005`.
No field, action, measure, or coefficient is identified merely by symbol or
dimension. The SU(2) bridge is only the pure-gauge specialization
`C_background=-b_SU2/(32*pi^2)` with `C_A=2`; it does not rederive the general
group or matter weights of `C-RGE-005`. The author may implement and propose
but may not independently accept or merge any claim.

## Candidate and Selection Record
The candidates were fixed before implementation. Candidate A derives typed
bridges through the canonical gauge-beta and SU(3) generator APIs and supplies
claim-specific review packets. Candidate B adds documentation-only links while
retaining literal coefficients and generators. Candidate C promotes the
existing composite P227/P229 narratives unchanged. Selection criteria, in
order, are accepted-dependency closure, convention typing, exact derivation,
mutation sensitivity, narrow claim semantics, consumer reuse, and maintenance
cost. Candidate A is selected structurally; B leaves the work as an island and
C exceeds the verified scope. There is no empirical comparator and therefore
no comparator gate.

## Claim Delta
The reserved P227 identifiers remain `C-WLN-001`, `C-WLN-002`, `C-WLN-003`,
`C-LOR-001`, and `C-LOR-002`. Their review packet must state the exact
massive/massless domains, local constant-e gauge condition, nontrivial massive
Weyl obstruction, H2/H3 scope, and matrix-commutator-only little-group scope.
The reserved P229 identifiers `C-CMV-001` and `C-CMV-003` are narrowed to the
conditional SU(2) spectrum, exact heat-kernel/log coefficient, declared-scheme
minimum, imaginary-part magnitude, exact sector cancellation, and canonical
SU(3) T3 root charges. `C-CMV-002` and `C-CMV-004` stay proposed and unpromoted
because the independent two-loop master-formula derivation and the actual
inhomogeneous-background fluctuation determinant do not exist in the current
evidence.

## Decomposition
Work proceeds through these dependency-ordered steps and continues after failed attempts.

1. [x] Recall and source verification.
2. [x] Candidate and selection-criteria preregistration.
3. [x] Importable bridge implementation.
4. [x] Verifier and sensitivity audit.
5. [x] Framework-fit and downstream replay.
6. [ ] Independent claim review and any later promotion transaction.

## Attempts
Attempts are append-only and individually reproducible. A failed row names the
diagnosed mechanism and the next materially different action.

| Attempt | Candidate or repair | Artifact and command | Verdict | Mechanism | Next attempt |
| --- | --- | --- | --- | --- | --- |
| 0001 | Candidate A typed bridges and narrow review packet | `src/substrate_framework/chromomagnetic_background.py`, `chromomagnetic_sectors.py`, seven claim-review requests, and two proposal verifiers | Strong milestone | The exact SU(2) beta/background bridge closes, canonical SU(3) roots replace a literal tuple, and each candidate has an explicit review boundary | Replay all affected consumers and the final repository gate |
| 0002 | Repair exact spectrum encoding | `pytest tests/test_chromomagnetic_background.py::test_algebraic_spectrum_matches_encoder` and the P229 narrow verifier | Passed | The new verifier exposed `2.0` and a `0.0` default that contaminated symbolic input with a SymPy Float; integer literals restore exact `-gB` while preserving numeric behavior | Retain the exact symbolic regression and proceed to full replay |

## Validation
Validation covers the exact bridge equations, wrong-Casimir and wrong-Cartan
mutations, P227 exact verifier, P229 consumers, proposal schema, public exports,
and final generated and repository state.

- Targeted scientific command: P227 verifier passed 14 exact checks; P229 narrow verifier passed 11 exact checks.
- Mutation and counterexample command: wrong mass sign, Weyl/reparametrization weights, little-group sign, Casimir, scale weight, ghost weight, and Cartan generator were rejected.
- Dependency replay: GitNexus reported low risk, 24 mapped changed symbols,
  no affected execution process, and the full repository suite passed all
  2,193 tests.
- Targeted tests: 54 focused worldline/chromomagnetic tests passed; one exact spectrum node was rerun after repairing the detected Float contamination and passed.
- `scripts/validate.sh --full`: passed all fixed workflow checks and 2,193
  tests in 313.90 seconds at the final scientific/code boundary.
- `git diff --check`: passed in a separate invocation; rerun after this
  evidence-only status update.

## Debt Ledger
Every defect or unsupported promise inside the proposed merge unit must be
discharged before handoff. The independent review and excluded P229 capstones
are campaign frontier, not merge debt.

| Debt | Introduced by | Why it is real | Discharge artifact | Status |
| --- | --- | --- | --- | --- |

## Results
The recent work is no longer an island at the proposal boundary. The new
`SU2ChromomagneticBetaBridge` derives `C_A=2` from the accepted Pauli-half
matrices, consumes C-RGE-005's exact beta ledger, and closes independently on
P229's `11/(48*pi^2)` background coefficient. The SU(3) T3 charged roots now
come from C-LIE-001's canonical generator API rather than a literal tuple.
P227's local constant-e and nontrivial massive-Weyl scopes are pinned by the
proposal verifier. Seven claim-specific review requests give a distinct
reviewer exact statements, dependencies, mutations, exclusions, and promotion
boundaries.

The audit also found and repaired exactness contamination in the charged-mode
encoder: symbolic input had inherited `0.0`/`2.0` literals and produced a
SymPy Float. The exact verifier now requires `-gB`. C-CMV-002 and C-CMV-004
remain research frontier for the independent two-loop derivation and actual
inhomogeneous-background requantization described in the P229 frontier record.

## Canonicalization
This author PR may change canonical source APIs, tests, P227/P229 proposal
records, and reviewer-request artifacts. It does not edit accepted claims,
release manifests, immutable campaigns, generated documentation, or accepted
claim memory before the independent verdict.

## Done Gate
The author implementation, review packets, dependency replay, and final full
validation are complete with an empty merge-scope debt ledger. The handoff
remains active only until its commit, PR, and issue comment exist. The
scientific promotion objective on issue #59 remains open until independent
claim decisions and any accepted promotion transaction land.

## Cross-References
Canonical issue: https://github.com/vantasnerdan/substrate-framework/issues/59.
Source proposals: `proposals/P227-einbein-worldline-harvest/` and
`proposals/P229-preparata-qcd-vacuum-audit/`. Source APIs are under
`src/substrate_framework/{relativistic_particle,lorentz_orbits,lorentz_little_groups,chromomagnetic_background,chromomagnetic_sectors}.py`.
