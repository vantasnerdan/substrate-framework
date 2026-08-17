---
description: Extract exact einbein worldline and Lorentz geometry atoms into claim-ready framework APIs
author: codex
created: '2026-08-11T14:14:03+02:00'
updated: '2026-08-17T20:45:00Z'
tags:
- substrate-framework
- campaign-proposal
- einbein
- worldline
category: proposals
confidence: established
status: archived
---

## Question and Positive Deliverable
This campaign asks whether the exact reusable worldline and Lorentz-geometry atoms in merged PRs #38 and #39 can become pure, dimension-aware, mutation-sensitive framework APIs. The positive deliverable is a focused harvest PR with importable implementations, independent exact tests, replayed consumers, and individually reviewable proposed claims; artifact merge does not itself promote a claim.

## Base Release and Provenance
The accepted base is release `v0.159.0`; the implementation branch starts at `substrate-framework@0b79c54443ed79e51b112dab9c26d10c2a102d49`. Source objects are the merged 2+1D and 3+1D tutorial Markdown and test suites from PRs #38 and #39, `src/substrate_framework/lorentz_orbits.py`, and the geometry functions in `src/substrate_framework/covariant_sine_gordon_action.py`. Related accepted claims C-SG-008, C-GW-002, C-MAX-001, C-GOR-001, and C-KIN-001 are authority and scope cross-checks, not conclusions imported into the proposed claims.

## Invariants, Conventions, and Allowed Imports
The campaign preserves the mostly-plus signature, Lorentz matrices acting on contravariant vectors, canonical momentum as a covector, positive `e` and `c0`, positive `m` only on the massive-elimination branch, and a separate `m=0` specialization. SymPy exact algebra and repository exact-input helpers are permitted. The worldline Weyl result is restricted to its displayed action transformation and cannot be called a Maxwell, Yang-Mills, photon, graviton, parity, helicity, anyon, continuous-spin, vDVZ, or substrate derivation.

## Candidate Preregistration
The candidates distinguish reusable representations rather than alternative physics mechanisms.

| Candidate | Description | Assumptions | Parameters | Framework-fit prediction | Decisive test |
| --- | --- | --- | --- | --- | --- |
| A | Compose a scalar-contraction worldline ledger with generic exact local-metric and Lorentz-algebra APIs. | The local velocity norm and metric data are typed inputs. | No new physical parameter. | Small APIs serve both dimensions and preserve explicit branches. | Both tutorial suites and independent exact mutations consume the same canonical objects without convention loss. |
| B | Build a fully coordinate-indexed SymPy action using applied coordinate functions and derive every object through `euler_equations`. | Coordinate charts and applied functions are supplied. | No new physical parameter, but a larger symbolic surface. | Strong direct derivation but higher cost and poorer composability for simple algebraic identities. | Exact derivation remains tractable for generic 3x3 and 4x4 metrics and replay time stays acceptable. |
| C | Keep the document-local formulas and tests as the implementation boundary. | Each document owns duplicate helpers. | None. | Lowest edit cost but fails importability and consumer integration. | GitNexus and import scans remain disconnected, so the issue success gate fails. |

## Selection Criteria and Blinding
Selection is ordered by accepted-convention compatibility, exact massive/massless branch separation, assumption and parameter economy, reuse across both tutorial dimensions and existing consumers, mutation sensitivity, independent rederivability, and maintenance cost. There is no empirical comparator. Candidate A is structurally preferred before implementation because the core auxiliary-field algebra depends only on the contracted norm while metric and Lorentz operations remain separately typed; Candidate B remains the independent derivation route, and Candidate C is retained as the disconnected baseline.

## Proposed Claim Delta
The provisional keys were collision-searched across governance, releases, campaigns, proposals, migration, memory, documentation, source, and tests before allocation. `C-WLN-001` proposes the massive auxiliary-field equation, positive elimination, square-root recovery, covector momentum, and pure-constraint Hamiltonian. `C-WLN-002` proposes the zero-mass constraint and non-affine-to-affine null-geodesic specialization. `C-WLN-003` proposes exact worldline reparametrization invariance and the narrow massless Weyl transformation with the massive obstruction explicit. `C-LOR-001` proposes only the exact induced `H^2` and `H^3` future-unit-timelike-vector orbit metrics. `C-LOR-002` proposes only the displayed 2+1D and 3+1D massive/null stabilizer Lie algebras. All have no accepted-claim dependency; the named accepted claims are consumers or scope ceilings, not premises.

## Implementation and Oracle Plan
Reusable APIs will live in `relativistic_particle.py`, `pseudo_riemannian.py`, `lorentz_orbits.py`, and `lorentz_little_groups.py`, with explicit public exports and no import-time calculation. The canonical route uses exact SymPy expressions and matrices. Independent tests reconstruct Euler derivatives, Legendre transforms, Christoffel contractions, orbit pullbacks, and commutators without calling the canonical helper that produces the compared quantity. Mutations cover the mass-term sign, reparametrization weight, covector/vector placement, zero-mass use of the positive massive root, Weyl weight, and little-group generator signs. Numeric tutorial checks remain regression evidence only where exact algebra is already available. Direct, imported, and dynamic legacy `np.trapz` access is preflighted even though no quadrature is planned. Consumer replay includes the two tutorial suites, covariant sine-Gordon geometry tests, Lorentz-orbit tests, related canonical import tests, and the full repository validator because public exports and shared geometry are affected.

## Attempts and Continuation
Attempts are append-only; the first implementation attempt begins only after this prose and its matching proposal manifest pass the schema gate.

| Attempt | Candidate or repair | Artifact and command | Verdict | Mechanism | Next attempt |
| --- | --- | --- | --- | --- | --- |
| 0001 | Candidate A canonical extraction | `src/substrate_framework/{relativistic_particle,pseudo_riemannian,lorentz_little_groups}.py`; generalized `lorentz_orbits.py`; exact targeted replay and `campaigns/P227-einbein-worldline-harvest/verify.py` | Strong milestone | The scalar-contraction ledger, generic geometry layer, and narrowly typed stabilizer ledgers served both tutorial dimensions. Independent derivations and load-bearing mutations all disagreed when signs, weights, or index placement changed. | Submit the artifact boundary for distinct-agent review; do not promote claims in the author PR. |

## Unit-Level Harvest Disposition
The einbein action and positive massive elimination, massless specialization, reparametrization/constant-e gauge/Hamiltonian family, and narrow worldline Weyl identity are implemented as claim-ready exact ledgers under `relativistic_particle.py`. The H2/H3 future-unit-timelike-vector orbit geometry is integrated by dimensionally generalizing the existing unpromoted `lorentz_orbits.py`, while coordinate Christoffel and Ricci construction is refactored into `pseudo_riemannian.py` and replayed through the existing covariant sine-Gordon wrappers. The displayed massive and massless little-group matrices and commutators are harvested narrowly in `lorentz_little_groups.py`. Printed Christoffel tables remain tutorial regression/documentation rather than a new scientific claim. Representation-classification, helicity, parity, anyon, continuous-spin, photon, graviton, Maxwell/Yang-Mills conformal, and substrate interpretations remain only cited tutorial narrative and PR history; they are neither encoded nor proposed for promotion here.

## Verification Results
The proposal verifier originally executed 13 exact checks, including explicit rejection of wrong mass-term and little-group signs. The affected targeted replay passed 89 tests across the new canonical modules, both tutorial consumers, Lorentz orbits, and the covariant sine-Gordon geometry consumer. GitNexus reported low integration risk, no affected execution processes, and only the expected touched geometry/tutorial symbols; the new files were covered directly because the pre-change index could not map previously nonexistent symbols. `scripts/validate.sh --full` passed all fixed workflow checks and all 2,126 repository tests in 230.50 seconds. A final review then repaired concrete-numeric einbein elimination without changing the symbolic result; that boundary passed all fixed checks plus the 60 directly affected worldline and tutorial tests. The issue #59 pass adds a fourteenth exact check for the positive local constant-e gauge rate and narrows the massive Weyl check to distinguish `Omega=1` from a nontrivial rescaling.

## Debt Ledger
This ledger tracks new imports, convention conflicts, unresolved verifier sensitivity, broken consumers, and unsupported interpretation inside the proposed harvest. It remains empty after targeted implementation replay; the proposal deliberately retains claim adjudication as campaign frontier rather than misclassifying it as artifact debt.

| Debt | Introduced by | Why it is real | Discharge artifact | Status |
| --- | --- | --- | --- | --- |

## Review and Promotion Plan
A distinct agent will review the PR claim by claim using the raw APIs, exact tests, mutations, and affected-consumer replay. The author PR leaves `governance/claims.yaml`, releases, generated docs, accepted memory, and immutable campaigns unchanged. If the reviewer accepts any proposed claim, promotion must be a separate complete governance transaction with adjudication evidence, registry/release closure, generated documentation, accepted-memory synchronization, and full validation. Rejected or narrowed claims retain their provisional identifiers and proposal history.

## Issue 59 Review Integration
The 2026-08-12 integration pass created one author-side review request for each
of C-WLN-001, C-WLN-002, C-WLN-003, C-LOR-001, and C-LOR-002. The packets
state exact domains and exclusions rather than accepting the claims. In
particular, the constant-e construction is local for smooth positive einbein
unless global interval coverage is separately shown, and the massive Weyl
term is an obstruction only for nonidentity positive rescalings. The P227
verifier now checks those distinctions directly. Issue #59 remains open for a
distinct review and any later governance transaction.

## Done Gate
The harvest PR is ready only when the importable objects exist, the exact and mutation oracles pass, affected consumers replay, the public API and convention maps are documented, repository validation succeeds, and debt inside the harvest boundary is empty. The scientific campaign remains active until an independent review either promotes or rejects each provisional claim; merge provenance alone is not completion.

## Adjudication Result — 2026-08-17

P227 moved unchanged into `campaigns/P227-einbein-worldline-harvest/` for
adjudication. The first prose-review draft failed audit because its C-WLN-002
Euler residual and C-LOR-002 generator products were false. Attempt 0001
preserves that failure. A fresh executable independent review repaired the
oracle, passed 51 exact mutation-sensitive checks, and agreed with the
canonical APIs. C-WLN-001 through C-WLN-003 and C-LOR-001 through C-LOR-002
were then promoted individually in release `v0.160.0` with no dependencies or
supersession edges. The generated registry, documentation, and accepted memory
are the current authority; this proposal record remains historical provenance.
