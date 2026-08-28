---
description: Earn an isolated fixed-J M5 clock and determine its two-body interaction law
author: Main
created: '2026-08-28T10:40:00+00:00'
updated: '2026-08-28T14:20:00+00:00'
tags:
- substrate-framework
- campaign-proposal
- m5
- isolated-clock
category: proposals
confidence: exploratory
status: active
---

## Question and Positive Deliverable
P247 asks whether the accepted spectral-Cartan M5 model class supplies a genuinely isolated fixed-angular-momentum clock — an asymptotically natural, finite-energy stationary field at fixed internal angular momentum `J` with finite positive inertia, nonzero frequency, and stability against the full admissible perturbation class — and, if it does, what leading two-body interaction law two well-separated relaxed clocks obey. The positive deliverable is at least one such one-clock construction, a nonsingular independent two-clock collective-coordinate sector, and a derived-and-measured leading interaction (channel, sign, falloff) — or, for each admissible candidate that fails, a mechanism-named obstruction that redirects the search. This is issue #178's milestone ladder; no boundary condition, ansatz, topology, action extension, or numerical method is prescribed.

## Base Release and Provenance
The accepted base is release v0.168.0 at `substrate-framework@552f2f6`; development runs from current `main` (`6afc49e`). Accepted inputs read for this contract: C-M5S-001 through C-M5S-014 in `governance/claims.yaml`, the canonical modules `src/substrate_framework/m5_covariant_action.py`, `m5_stationary_fields.py`, `m5_fluctuation_spectrum.py`, `m5_kinetic_axis.py`, and the merged P240 attempt records 0041-0043 under `proposals/P240-m5-kinetic-axis/attempts/`. The reported attempt 0044 de-boxing verdict (PR #153 comment 2026-08-22 and issue #151 comment 2026-08-22) is report prose only: its artifacts were appended to `research/p240-phase1-stability` after merge, the branch was auto-deleted, and neither `proposals/P240-m5-kinetic-axis/attempts/0044/` nor the referenced memory entry `memory/giuliano/efforts/p240-deboxing-verdict.md` exists in the repository. Milestone 0 of this campaign reproduces and materializes that verdict from merged machinery before any use.

## Source Inventory and Access Gate
All load-bearing sources are in hand in the repository.

| Source | Access status | Extracted claims (with page/eq) |
| --- | --- | --- |
| governance/claims.yaml C-M5S-001..014 | in-hand repository source | census, cutoff, composition, consumer, pair-interaction, spectrum, zero-point, constraint, stress, initial-data scopes |
| proposals/P240-m5-kinetic-axis/attempts/0041 (cpu_energy.py, solve_radial_1d.py, radius_stability_scan.py, certify_morse_index.py, result.yaml) | in-hand repository source | committed torch CPU functional: biaxial hedgehog basis with structural wall pin, potential `-Tr(X^2)/2 - Tr(X^3) + (Tr X^2)^2 + 1/2`, clock generator `[n-hat]_x`, inertia density `4*sum_d ||[GX-XG, d_d X]||^2`, `omega = 1/(2I)`, `E_J = E_stat + 1/(4I)` |
| proposals/P240-m5-kinetic-axis/attempts/0042 (phase1_ladder.py, xspace_energy.py, largeR-roots.json, result.yaml) | in-hand repository source | R-ladder root production and independent x-space energy route |
| proposals/P240-m5-kinetic-axis/attempts/0043 (pair_oracle.py, result.yaml) | in-hand repository source | boxed pair oracle behind accepted C-M5S-008 |
| PR #153 comment 2026-08-22T00:13:32Z and issue #151 comment 2026-08-22T00:13:56Z | open URL (GitHub), quoted in issue #178 | reported 0044 verdict: symbolic V2 = diag(5,3,3,0,0,6)/2 local frame; linear E(R) growth ~0.33/unit on R=12..18; lambda_min(A) < 0 flattening at -(5+-2)e-7, Morse index 1, purely split-channel; stability window R in [~8,~34]; free-wall R=8: E=50.44 vs pinned 60.17, wall value 0.99625, I 0.67638 vs 0.67662, omega 0.73923 vs 0.73929; aliasing gate (quadrature >= 4x basis order, cross-quadrature energy reproduction to 1e-6) |

## Invariants, Conventions, and Allowed Imports
The campaign preserves the accepted P239 spectral-Cartan action and its pinned spectrum potential, the committed biaxial hedgehog conventions (mostly-plus signature, `omega = 1/(2I)`, `E_J = E_stat + J^2/(4I)`, the structural wall pin at the rank-1 projector `n n^T` wherever the pinned machinery is run), the accepted confined-clock claims C-M5S-001..014 in their stated confined and gravitational scopes, and the rule that no accepted confined-clock claim is reinterpreted as an isolated-particle claim. Comparator blinding is impossible for milestone 0 (the reported 0044 numbers are published in the comments); those numbers are therefore declared transfer references, and all reproduction gates below were frozen before this campaign computed any new number. P244 attempt 0010's de-boxing ladder for the zero-point shift and its preregistered R=16/18 extension are adjacent evidence about a different quantity and are not substituted for the milestone-0 reconstruction. The stationary axisymmetric Einstein-M5 co-solve is owned by #177/#174 and stays outside this campaign.

## Candidate Preregistration
Three candidate families are frozen before any pair-force comparator is inspected, per issue #178; the starting families are routes, not commitments.

| Candidate | Description | Assumptions | Parameters | Framework-fit prediction | Decisive test |
| --- | --- | --- | --- | --- | --- |
| A | Current-action, broader-sector route: keep the accepted positive fixed-J construction and search beyond the tested radial realization — the two exactly potential-flat shear channels (np, n-a) reported after P240 as the only long-range multipole candidates, plus full angular/topological sectors and asymptotic configurations | none beyond accepted canon | none | shear channels carry gradient-only stiffness, so they are the only candidate long-range sector inside the accepted action | de-boxed shear-channel far-field construction with finite energy and a signed falloff |
| B | Minimal positive extension route: add the smallest local Lorentz-covariant structure that vanishes on the complete arbitrary static 3x3 sector yet changes isolated asymptotics or split-channel stability | one new local term, no fitted frequency or force | minimal, declared | extension must recover the static sector exactly and repair the split-channel sign or asymptotics | exact static-sector recovery plus de-boxed one-clock existence and stability |
| C | Alternative completion/topology route: a different Lorentz-covariant 4x4 completion or topological clock sector under the same criteria | alternative action or sector, no fitted answer | minimal, declared | must match the static 3x3 recovery, boundedness, symmetry, and parameter economy of the accepted baseline | same gates as A/B under the frozen criteria |

## Selection Criteria and Blinding
Structural selection criteria, frozen before any pair-force comparator is opened: bounded Hamiltonian and regular charge-frequency map; local Lorentz covariance; exact or explicitly challenged 3x3/static recovery; natural finite-energy asymptotics; assumption and parameter economy; topology and symmetry fit; implementability; predictive reach. No pair-force number exists yet in this campaign; when one is produced, the milestone-4 field observable (surface/multipole/stress or cross-inertia) is derived first and the fixed-J Legendre algebra comparison happens only after that observable is established.

## Proposed Claim Delta
P247 reserves C-M5S-015 for the milestone-0 adjudication: in the committed P240 model class at the accepted base, the de-boxed limit of the confined radial fixed-J clock branch is not a finite-energy stationary isolated clock — the reproduction states the exact scope (pinned and free-wall radial branches, the tested R ranges, the pencil/stability verdicts, and the two exactly flat shear channels) with reproduced numeric evidence. Collision search performed 2026-08-28: `governance/claims.yaml` ends at C-M5S-014, and registry, proposals, campaigns, and durable memory contain no P247 or C-M5S-015 identifier. Later claim identifiers (isolated existence, stability, two-clock sector, interaction law) are named only in manifest revisions once their constructions are concrete.

## Implementation and Oracle Plan
Attempt 0001 (milestone 0) reconstructs the reported 0044 computation from merged machinery: exact SymPy second variation of the committed potential about the rank-1 projector background in the local frame (G1); the clean R-ladder on the committed torch functional with the aliasing gate — solving quadrature at least 4x basis order and cross-quadrature energy reproduction to 1e-6 — measuring E(R) growth, constrained-second-variation lambda_min, Morse index, and channel fractions (G2); the pencil stability window (G3); and the free-wall variation with a relaxed wall amplitude (G4). All gates were frozen in this contract before execution; stdout is captured into `proposals/P247-isolated-clock-two-body-law/attempts/0001/` on first execution. Later attempts follow the issue ladder: isolated one-clock existence (finite-energy asymptotics from the candidate action on unbounded/matched/converged domains), constrained second variation with deflated moduli and symmetry-breaking perturbations, the independent two-clock 2x2 inertia sector, and the conditioned far-field interaction observable. Reusable asymptotic, solver, and stability machinery is extracted into `src/substrate_framework/` at promotion. Campaign verifiers run with `PYTHONPATH=src` against canonical modules; committed attempt infrastructure (0041/0042 helper APIs and roots) is imported as source-audited infrastructure following the P244 precedent.

## Attempts and Continuation
Attempts are append-only under `proposals/P247-isolated-clock-two-body-law/attempts/`. Each attempt ends in a decisive verdict — established, refuted with mechanism, or blocked with the one missing construction named. A no-go for the tested radial branch advances the search and routes the next attempt to the surviving candidate families; it does not close the issue. The aliasing gate adopted from the reported 0044 work (spurious fixed-order large-R stationary points with divergent independent-quadrature energies) binds every later root acceptance in this campaign.

Attempt 0001 (milestone 0, 2026-08-28) reached a MIXED verdict with the full adjudication materialized in `proposals/P247-isolated-clock-two-body-law/attempts/0001/result.yaml`. Established exactly: the potential second variation about the rank-1 projector background is diag(5,3,3,0,0,6)/2 with exactly flat np/na shear channels; the de-boxed energy grows monotonically without plateau along the clean aliasing-gated ladder R=12..18 (55.1042 to 59.2462, fitted 0.685 per unit R, R^2=0.9984), so the branch has no finite-energy isolated limit; the free-wall R=8 endpoint reproduces (E=50.5233, wall value 0.99631736, index 0, clock sector invariant). Not reproduced, each with a named mechanism: the reported persistent negative split-channel mode (actual: lambda_min(A) = -1.8751e-06 at R=12 and -7.8309e-07 at R=13 with q-channel-dominated negative directions, then the positive +-1.7e-8 resolution floor with index 0 from R=14); the bounded stability window (actual: lower edges rise with background radius, upper edges unbounded for R>=14); the 0.33 slope (actual 0.685); and the pinned 60.17 reference (actual 50.52966035 from every committed seed path; pinning costs 0.012 percent, not 16 percent).

Attempt 0002 (milestone 1, 2026-08-28, branch campaign/p247-attempt-0002-two-body-law,
no PR per owner direction) verdict DECISIVE_NO_GO_WITH_STABILITY_RESTORATION.
Milestone 1 closed on the frozen PERSISTENT_GROWTH band: C[c] 461 -> 1197
through R=30 without resolved saturation (saturating fit loses by dAICc 7.6),
E(R) slope steepens to 0.855 - scoped isolation no-go for the radial de-boxed
family, route to candidate B. The preregistered basis-order cross-check (16/20/24,
Richardson with budgets, jitter sign tests) VOIDED attempt 0001's order-16 floor
claims at every rung: the refined lambda_min(A) is negative above budget everywhere
(R=14 -8.6e-7, R=16 -3.4e-7, R=18 -1.5e-7, R=12 -8.7e-6 with index 2), restoring
the reported 0044 magnitude (-5+-2)e-7 while refuting its split-channel attribution
(the negative mode is q-channel at fraction >= 0.99; the split channel is the
near-null positive direction converging to +0). G3 window reopens (order-24 pencil
recomputation = attempt-0003 gate W1); candidate-B term enumeration on paper is
gate W2. Shear probe (direct projector construction; the projector is not
representable in the modal ansatz): G1 potential-flatness and the
gradient-only-stiffness premise confirmed at the projector (static a2 +11.9/+8.9
np/na, potential ~0); on the branch the shear channels stay statically stable
(+9.2 to +3.8) while softening toward marginality. C-M5S-015 promotion postponed
to issue completion per owner direction; its flatness-to-floor draft wording is
superseded by attempt 0002 (draft never accepted, no registry change).

Strongest supported statement (after attempt 0002): the radial de-boxed branch
of the committed model class fails both isolation (curvature class grows without
saturating through R=30) and stability (a genuine q-channel negative mode at
every rung, magnitude flattening near -(5e-7..1e-6) at refined basis order);
the exactly flat np/na shear channels are potential-flat with gradient-only
stiffness, statically stable on the branch and softening toward the projector
marginality. Candidate A's radial-family route is exhausted; the search routes
to candidate B (minimal extension) with candidate A's shear sector as the
measured long-range phenomenology to preserve. The reported 0044 artifacts must
still not be cited as accepted evidence; this materialized reproduction remains
the citable record, and the 0044 adjudication must be corrected at the next
materialization (magnitude restored, channel attribution refuted).

## Debt Ledger
The milestone-0 reproduction inherits the reported verdict's numbers only as transfer references; any gate that fails is recorded as MIXED with the named mechanism rather than silently narrowed. Open frontier (not debt): the full-action scope beyond the tested radial branch, the isolated-existence question for candidates A-C, the two-clock sector, and the interaction law. The absent 0044 source artifacts are external provenance facts, recorded here and in the attempt record, not hidden assumptions inside any proposed claim.

## Review and Promotion Plan
C-M5S-015 receives one independent claim review on the frozen transaction (claim text, reproduced evidence, declared dependencies C-M5S-001..014 scopes, affected consumers) plus at most one correction check. Promotion extracts reusable verification machinery into `src/substrate_framework/`, updates `governance/claims.yaml` and the release manifest, runs `scripts/render_docs.py`, synchronizes accepted-claim memory, and records one content-addressed validation receipt at the frozen boundary.

## Done Gate
P247 closes only when one candidate supplies the isolated finite-energy fixed-J clock, the full declared stability result, a nonsingular independent two-clock sector, and a resolved leading interaction law with claim-level review and promotion. A reproduced de-boxing no-go, a confined replacement, or a conditional force formula is progress evidence that keeps the issue open, and each failing milestone names the construction that unblocks the next attempt.
