---
description: 'P251: derive Cosserat (micropolar) equations from coarse-grained Euler with vortex-structure rotational degrees of freedom (issue 198, Comparsi item via Dan)'
author: giuliano
created: '2026-09-03T14:40:00+00:00'
updated: '2026-09-05T06:00:39+02:00'
tags:
- substrate-framework
- campaign-proposal
- p251
- cosserat
- euler
- vortex
- issue-198
category: proposals
confidence: exploratory
status: active
---

## Question and Positive Deliverable

The campaign answers whether coarse-graining a declared ensemble of vortex structures — whose rotational degrees of freedom are carried by exact stationary knotted vortex tubes of the Enciso–Peralta-Salas type — under constant-density incompressible Euler yields an exact conditional derivation of the linear Cosserat (micropolar) equations. The positive deliverable is the claim ladder C-CST-001..007 with computed effective moduli, identified balance laws, two-branch dispersion, the no-Cosserat contrast side, and the EPS conditional bridge; a no-go, residual, or honest obstruction does not complete this campaign.

## Obligation Graph and Closure Map

The objective decomposes into seven dependency-ordered obligations, tracked in the manifest's `obligation_graph`. N1 (frame-transport kinematics) and N2 (Biot–Savart micro-moduli) are the roots; N3 (ensemble Cauchy–Born with triads) consumes N2; N4 (balance-law identification) consumes N3; N5 (dispersion) and N6 (no-Cosserat side) consume N4/N3; N7 (EPS conditional bridge with a realized Beltrami example) is last and consumes N2 and N5. Each node records positive intent, pass licenses, failure scope, and its unlocked successors; no node, subclaim, or utility opens a PR or terminates the campaign.

## Base Release and Provenance

The accepted base release is v0.171.0 pinned at substrate-framework@e24497b0272f8582a71d4189581db818a0e276aa; the canonical goal issue is #198 (vantasnerdan/substrate-framework), opened 2026-09-03. Modules actually read: `src/substrate_framework/homogenization.py` (declared affine filament homogenization, sphere moments, coefficient matching — conditional unpromoted P242 infrastructure), `src/substrate_framework/numerics.py` and `verification.py` (CheckLedger/IVP evidence APIs), `proposals/P242-elastodynamics-from-euler/proposal.yaml` and its verifiers (ladder and mutation patterns mirrored here). P242's claims C-ELS-001..007 are proposed-but-unpromoted and are imported only as conditional machinery.

## Source Inventory and Access Gate

All load-bearing external sources are in hand, verified 2026-09-03 (PDF saved, text layer extracted to sidecar, md5 recorded in `sources/MD5SUMS`, provenance in `sources/PROVENANCE.md`). The four load-bearing papers are committed: arXiv:1210.6271 (knotted vortex tubes, existence backbone), arXiv:1003.3122 (knotted vortex lines), arXiv:1505.01605 (Beltrami eigenfields — realized-example pool for N7), arXiv:2103.14458 (in-hand statement of the target linear micropolar system plus the coarse-grain-rotations methodology precedent). arXiv:2401.09506 is context-only (37 MB PDF hash-pinned locally, not committed). MIT 2.25 §9 is the declared classical input, re-verified open 2026-09-03. Eringen 1966 and Cosserat 1909 are cited provenance, not load-bearing. No source is missing; nothing blocks.

## Invariants, Conventions, and Allowed Imports

Constant-density incompressible Euler is the sole microscopic dynamics; every emergent statement is conditional on declared, independently falsifiable ensemble premises with no hidden fitted constant — bending/twist moduli and microinertia are derived from the Biot–Savart self-energy, never postulated. Material frame indifference, isotropy, and SO(3) frame covariance hold; the fluid limit recovers Euler and the structure-free limit recovers the P242 Navier–Cauchy sector. Exact claims use SymPy; dynamical claims use SciPy with declared refinement and mutation probes; incompressible and compressible branches are never mixed. The incompressible constraint is declared by the ensemble itself (P242 precedent). A candidate defect versus a canonical inconsistency is distinguished by the N4 identification test: failure confined to the candidate's moduli or frame choice is a candidate defect; a failure that survives candidate substitution independently opens a separate foundational-revision proposal.

## Frozen Candidate Universe

The in-scope route families derive from the original objective plus the invariants and permitted imports: (family 1) filament-triad affine transport, where rotational DOF are centerline-frame bend/twist; (family 2) core-polarization patches, where rotational DOF are cross-section vorticity-polarization rotations; (family 3) ergodic no-topology closure, the no-Cosserat contrast side. The universe expands append-only with failure-generated or historical concepts and does not shrink without user approval.

## Candidate Preregistration

Two candidates are registered for mechanism selection. Candidate A (Bishop-triad affine transport) carries rotational DOF on segment centerline frames, derives moduli from the N2 self-energy via affine transport of curved/twisted segments, and predicts couple stress proportional to the twist-modulus moment closure. Candidate B (core-polarization patch) carries rotational DOF in the tube cross-section's vorticity-polarization profile, with torsion stiffness from profile-distortion energy under the same Biot–Savart functional. A's assumptions are thinner at the segment level (no cross-section profile model); B introduces a cross-section polarization ansatz but may realize the microrotation variable more directly. Decisive test for both: exact coefficient matching (N3) plus balance-law identification (N4) with mutation probes, then dispersion (N5).

## Selection Criteria and Blinding

Frozen ordering: (1) framework invariant compatibility (Euler microdynamics, Kelvin frozen-in, frame indifference); (2) assumption and parameter economy (no postulated modulus; fewer declared premises); (3) correct limits and cross-sector composition (fluid limit Euler; structure-free C-ELS Navier–Cauchy; rigid rotation classical spin); (4) symmetry and dimensional consistency of every effective coefficient; (5) implementability with claim-appropriate oracles. No empirical comparator exists; the structural target form is identified from in-hand 2103.14458 before attempts, and the blinding point is the structural freeze of the N3/N4 moduli derivations, after which no numerical value may alter any candidate or criterion.

## Proposed Claim Delta

Claims C-CST-001..007 are proposed (collision search over registry, campaigns, proposals, and durable memory on 2026-09-03 found no C-CST use; identifiers free). C-CST-001 frame-transport kinematics; C-CST-002 Biot–Savart micro-moduli; C-CST-003 ensemble moduli; C-CST-004 balance-law identification; C-CST-005 dispersion; C-CST-006 no-Cosserat side; C-CST-007 EPS conditional bridge. Dependencies: 001/002 are roots; 003 needs 002; 004 needs 003; 005 needs 004; 006 needs 003; 007 needs 002 and 005. Consumers: future substrate sectors wanting rotational continua; no accepted claim depends on these. None uses `supersedes`.

## Analytic Specification and Numerical Licenses

Per-obligation analytic specifications live in the manifest `license_chain` entries (object, symmetry, ensemble, variational functional, admissible space, representation coverage, observable, blocked numerical representation, permitted verdict). N1, N2, N3, N4 are exact claims whose oracles are SymPy identities with mutation and counterexample probes — no production numerics is licensed for them. N5 earns L4 (symbolic two-branch dispersion) before its IVP runs; the IVP uses `substrate_framework.numerics.solve_ivp` with declared tolerances, a single harmonic mode (nonperiodic-window contamination excluded by construction), refinement-stability evidence, and a scale-relative error budget. N6's Monte-Carlo average runs only after the exact isotropy argument, with declared seed, sample count, and scale-relative tolerance. N7's spectral sampling runs only after N2–N5 analytic closure (L5 gating), on a declared resolution ladder. Pre-gate sampling anywhere is `exploratory_only` and cannot select candidates or support claims.

## Implementation and Oracle Plan

Reusable pieces extend `src/substrate_framework/homogenization.py` patterns: triad/tangent joint moment closures, micropolar Cauchy–Born energy, coefficient-matched moduli, and balance-law residuals are new importable functions with tests, extracted only at harvest; campaign verifiers stay thin under `PYTHONPATH=src` importing `CheckLedger`, `numerics`, and homogenization APIs. Every exact claim is checked two-route (direct construction + independent symbolic route), mutated (each load-bearing input flipped → relevant check must fail), and probed for wrong conventions (sign of torsion term, skew-vs-symmetric stress split, log-cutoff placement). Units convention: ρ density [M/L³], Γ circulation [L²/T], a core radius, R outer cutoff, L_v tube length per volume [1/L²], moduli in micropolar units derived consistently; the convention block is written before each verifier and expectations are computed independently of the code under test. Verifier stdout is captured into `attempts/000N/` on first execution.

## Attempts and Continuation

Attempts are append-only under `proposals/P251-cosserat-from-vortex-euler/attempts/000N/` with candidate source, command, environment, stdout, verdict, and diagnosed mechanism; the manifest `route_frontier` tracks considered/tried/failure-generated/remaining routes. A failed route names its layer (implementation, numerics, representation, concept, target statement, accepted foundation) and triggers the next materially different attempt per the classification ladder; a candidate that conflicts with the framework is rejected and the next candidate runs. Execution runs in declared waves: wave 0 grounding (complete 2026-09-03), wave 1 exact ladder (N1–N4), wave 2 dynamics and contrast (N5, N6), wave 3 bridge (N7); a wave starts only when every declared input has settled.

## Debt Ledger

Currently empty: no hidden assumption, unexplained fitted parameter, unresolved promise, convention conflict, or broken consumer inside the proposed claims. Declared hypotheses (scale separation a ≪ window ≪ field scale, affine transport, inviscid persistence, isotropy with triad locking) are premises of the conditional claims, not debt. The uncommitted 37 MB context PDF is recorded provenance, not debt.

## Review and Promotion Plan

Each claim C-CST-001..007 gets one claim-level independent review with raw artifacts and acceptance criteria (reviewer must rederive the load-bearing step); the frozen transaction is the claim delta plus its evidence records and extracted APIs. Evidence roles: `exact_proof` for SymPy two-route derivations, `regression` for SciPy refinement runs, `applicability` for the N7 realized example, `provenance_only` for the EPS import. Promotion (if review accepts) extracts reusable APIs into `src/substrate_framework/` with tests, updates `governance/claims.yaml` and a pinned release manifest, moves the adjudicated record into `campaigns/`, runs `scripts/render_docs.py`, and synchronizes accepted-claim memory. The PR to main goes through fork `giuliano-vantasner/substrate-framework` with a distinct merger; `Advances #198` until the full objective completes.

## Scientific Exhaustion Certificate

Active and incomplete: no exhaustion is claimed at preregistration. The section fills only if positive success passes first, per AGENTS.md.

## Done Gate

The campaign PR gate opens only on complete positive success (all success-contract gates passed, debt ledger empty) or the scientific exhaustion certificate; until then the campaign checkpoints with commits on `research/p251-cosserat-from-vortex-euler` and does not open rung, milestone, or utility PRs. Federico Comparsi's sketch (requested; reply pending) lands as intake context before the structural freeze; a late arrival after freeze is recorded as intake evidence and may spawn an append-only candidate-universe expansion, never a silent criterion change.

## Intake Amendment 2026-09-03

Federico Comparsi (@Splinterrone) posted the item's derivation sketch on #198 (comment 5526924549) at 14:00 UTC, before any attempt ran. The sketch pins the structural target: his equations (1) and (2) are the standard linear micropolar system with coupling modulus alpha, couple moduli beta/gamma/epsilon, and microinertia J_i; N4 identification now targets that form, with Rudge 2103.14458 as the consistent secondary statement. His micro-decomposition v = V_macro + Phi x r + v' weighs toward the core-polarization family but does not change the frozen candidate set or criteria. The manifest source_inventory carries the full entry; notation register: our derivation writes skew(sigma) and epsilon_c to avoid the Levi-Civita collision present in the sketch's prose.

## Independent Review Correction 2026-09-04

PR #199 remains an active campaign rather than a terminal-success record. N1
is established, while N2 is the active obligation and licenses L1--L5 remain
unearned. The exact N3/N4/N5 algebra is retained as a useful conditional model,
but the Euler-to-Cosserat dependency edge is still open: exact relative
Green--Lagrange line stretch is invariant under a free frame rotation, so the
claimed `alpha = L_v T/6` term cannot be obtained by truncating the relative
rotation at first order and then retaining its quadratic norm. The next
construction is an Euler-derived core-polarization or inter-tube Biot--Savart
interaction whose exact second variation locks macro- and micro-rotation.

Bounded defects have been repaired in attempt 0025. The transverse determinant
now follows directly from the N4 operator and contains `4 alpha^2 k^2`; its
SciPy DOP853 refinement is executed. The structure-free limit removes the spin
coordinate when all of its action weights vanish, and the longitudinal branch
is labelled as a formal unrestricted micropolar extension rather than a
physical incompressible-Euler mode. N6 now establishes only zero coherent
signed response under its explicit ergodic closure, while retaining nonzero
quadratic frame fluctuations. N7 imports the hash-pinned EPS vortex-tube
existence statements and uses the elementary periodic Beltrami field only as a
steady-Euler equation regression, not as a localized or knotted tube.

The earlier empty-debt and terminal-completion statements are superseded by
this correction. The campaign frontier is a scientific construction, not a
refutation: close N2 and derive the missing frame-locking interaction, then
replay the conditional N3--N7 ladder and individual claim reviews.

## Completion continuation 0029–0034

The new exact Cartesian audit (0029) corrects the earlier review's reliance on
0019: original N2 velocity formulas have the correct Coriolis signs; 0019's
replacement operator reverses them. Exact differentiation of K1 in 0031 gives
the bending-branch constant 1/4-EulerGamma and preserves the m=2 coefficient
-1/6; signed Bessel roots replace the previous coarse residual scan.

Attempt 0028's rate coefficient has inertia units and is not an elastic alpha.
The repaired N3 verifier preserves the conditional angle-energy and moment
identities. Genuine angle dynamics now exists in 0030 (elliptic Euler patch
with prescribed ambient flow) and 0032 (closed six-vortex collective action).
0033 corrects the energy/helicity confusion in the crossing route, and 0034
derives a consistent collective field map plus a conditional long-wave
continuum correspondence with cell density L_cell=L_v/6.

The original exact smooth-Euler/EPS continuum objective remains active. Its
next construction is the finite-core, reaction-bearing affine-cell action
whose translation, orientation and spatial coefficients share the derived
Euler symplectic structure. Neither a green local oscillator nor an unbounded
ambient-flow model closes that objective; the new constructions are preserved
as executable progress on the campaign branch.

## Reviewed conditional construction checkpoint 0108--0113

Codex continuation constructs the full conditional finite-core alternative:
compact Euler angular pairs inside actual stationary EPS tubes, one canonical
material phase action, common-mean momentum reduction, nonzero physical
translation/core-angle coupling, positive shear/curvatures and explicit
population limits. Independent0108 accepts new C-CST-008..010 exactly as
frozen in0105;0113 records their release175 promotion. Original001..007
identifiers remain bound to the earlier statements, and the old five-scalar
straight-tube formula is not revived. The alternative role closure is in
0105/route-closure-map.md and0108/goal-closure-assessment.md.

This does not close the stronger unrestricted same-cell Euler realization.
The exact complement/current remains in0095,0109 gives controlled actual
parcel spin and translation, and0112 continues on the actual EPS periodic
core. The campaign branch remains active beyond the promotion checkpoint;
no terminal PR or full-unrestricted-objective verdict is asserted.

## Actual same-core continuation 0112--0117

0112/0114 construct actual Euler amplitude transport on the EPS periodic
core, a positive action in a fixed zero-winding geometric frame, and a
nonzero actual material-spin change over one core period.0115 supplies
the exact moving symplectic/complement bookkeeping used in that action.
0116 continues with identical periodic knotted cells and complete Bloch
stress transfer.0117 retains the actual spin AND symmetric shape-rate
dipoles in the centroid-plus-ambient Fourier momentum. With the declared
isotropic full response this yields a first-gradient physical centroid
response to the changing spin, even when the resolved point dipole cancels.
The stronger autonomous spin/action coefficient and full spatial-dynamics
join remain active; none of these observations is renamed C-CST-009's
constant constitutive law. Accepted008..010 are unchanged.
