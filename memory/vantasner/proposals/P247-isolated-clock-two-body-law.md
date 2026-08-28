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

Attempt 0003 (W1/W2/W3, 2026-08-28, same branch, no PR) verdict
DECISIVE_PER_GATE. W3 SIGN_CHANGE: at order 24 the q-channel static instability
is confined - lambda_min(A) = -2.08e-9 at R=20 (q-dominated), zero crossing in
(20,22), then +5e-10 split-channel floor through R=30 with Morse index 0; the
de-boxed family is statically stable for R >= 20 at refined resolution. W1
REFUTED_ORDER24 for the universal band: order-24 windows are background-relative
(lower edge 0.77-0.80 x R_own, finite upper edges for R <= 18; every background
stable at its own radius), self-consistent band [19.20, 29.96] fails the frozen
[8,34]+-20 percent gate at the lower edge while the reported upper edge is
rehabilitated in-gate. Refined roots carry order-improvement 0.7-12.7 percent
growing with R (order-16 roots have ~1e-2 stationary error); E(R) stays monotone
so the isolation no-go survives. W2 SELECTED candidate B: the boost-alignment
mass s_m = (m^2/2)((u^T eta xi)^2 - 1), m^2 = Lambda^2 (C-M5S-002 scale, no
fitted parameter) - masses the C-M5S-001 massless boost orbit (the named
isolation-failure mechanism), vanishes exactly on the static 3x3 sector,
isorotation-invariant (clock map unchanged); declared cost: explicit internal
boost-symmetry breaking. Execution honesty: a spurious divergent R=24 basin
(E=6.9e516) passed the relative alias gate at 3e-8; two preregistered
amendments (root-continuity gate 0.25 after measuring legitimate
order-improvement; R-continuation reseeding) plus one in-place band correction
(tag leak) are recorded in manifest and result. Attempt-0004 gates C1-C4
preregistered: reduced functional derivation, crossing bisection, one-clock
existence with Yukawa tail, pencil stability.

Strongest supported statement (after attempt 0003): the radial de-boxed branch
of the committed model class fails isolation (curvature class grows without
saturating; E monotone at refined roots through R=30) while its static
instability is confined to R <~ 20-21 and heals at refined resolution;
stability windows are background-relative, not universal. The accepted action's
massless boost orbit (C-M5S-001) is the named isolation mechanism, and the
selected candidate-B repair is the boost-alignment mass with m^2 = Lambda^2 -
zero new continuous parameters, exact static-sector survival, isorotation
invariance. The reported 0044 artifacts must still not be cited as accepted
evidence; the materialized reproduction remains the citable record (0044
magnitude and bounded-window-per-background claims rehabilitated, universal
band and split-channel attribution refuted).

Attempt 0004 (candidate-B construction, 2026-08-28, same branch, no PR)
verdict DECISIVE_PER_GATE. C1 IDENTIFIED: the boost sector is structurally
wired to the clock - the clock velocity sources the boost field (F_{01} mixed
term), the chi=0 regression against the committed functional is exact (diffs
0.0), and C1b pinned the normalization with zero new parameters (boost
projector-current metric entry exactly 1/16, accepted census eigenvalues
reproduced; mass m^2 = Lambda^2). C2 MEASURED: crossing bracket [20.0, 20.25]
- the static instability is confined to R <= 20.0(1). C3 NO_LOCALIZED_ROOT by
co-aligned source cancellation: the director-aligned boost field is a
spectator for the radial clock - grad E_J at chi=0 is exactly zero (boost
axis = director = rotation axis kills the contraction; the C1-r3 structural
note promoted to the full contraction). C4 not applicable (no new stationary
point). Successors registered from the mechanism: D_tilted_boost (fixed-axis
boost; breaks isorotation covariance - clock map rebuilt) and
E_localized_core (localized wall-core object; no new breaking, heavier
existence problem). Attempt-0005 gate D1: sympy full-contraction survival of
the linear source under non-co-aligned axes decides the routing.

Attempt 0005 (D1 routing, 2026-08-28, same branch, no PR) verdict
DECISIVE_D1_VANISHES: the linear boost source vanishes EXACTLY for both
co-aligned and fixed-axis boosts on the symmetric one-clock root (machine
precision; graph connectivity and numeric central derivative cross-checked;
finite-chi graph check confirms the functional is live). Mechanism: the
one-clock source is parity-odd on mu-symmetric roots, so the chi EOM about
chi = 0 is homogeneous - the boost-alignment mass cannot localize the
symmetric fixed-J clock at any order. One-clock no-go for the B-family
reduction, mechanism named. Reduction defect recorded: the truncated sigma
Hamiltonian evaluates negative at finite chi (not the positive Legendre
density). Redirect: F_pair_exchange (a parity-broken CLOCK PAIR sources the
gapped boost field - the interaction-law milestone has a live Yukawa route)
and G_fixed_axis_rotor (one-clock symmetry breaking by a fixed-axis rotor;
needs a non-spherical solver). Attempt-0006 gates F1 (pairing derivation) /
F2 (numeric check) / G1 (rotor specification) preregistered.

Attempt 0006 partial (D1', 2026-08-28): the one-clock no-go is COMPLETE -
the linear boost source vanishes exactly in EVERY angular channel (odd mu
1/3, even 0/2, combined) on the symmetric fixed-J root, not only the radial
channel probed in 0005. The boost-alignment mass cannot localize the
one-clock in the single-centered reduced class. Remaining extended-class
routes: F_pair_exchange (two-centered machinery) and route G (non-spherical
rotor solver, specified in attempts/0006/f1-g1-specification.md). The
extended functional gained a chi-angular-degrees parameter; the C3-0 exact
regression re-verified after the change.

Attempt 0006 CLOSED (2026-08-28): Yukawa/pair-exchange route REFUTED by two
independent mechanisms. (1) Exact Z2 evenness of the extended functional -
E[+chi] = E[-chi] at machine precision on symmetric AND randomized
backgrounds, all angular channels - so chi = 0 is the exact classical
minimum for EVERY configuration: no static pair exchanges through the boost
channel at any order (F_pair_exchange dead; route G's boost-tail dead). (2)
Box-growth channel-attributed to the S-potential sector (accelerating:
increments 0.29 -> 0.73 per DeltaR = 2 over R = 20..30; curvature
saturating; boost sector exactly 0 at chi = 0) - candidate B masses a
classically inert channel and cannot cure the static isolation no-go. (3)
No S row fits a massive tail in the accessible box (exp_resid 0.35-0.64) -
the Yukawa pairing premise fails; the committed-class two-clock law is the
slowly-decaying S-overlap (confined anchor +456.6 d^-1.696, C-M5S-008).
Remaining routes: de-boxed S-overlap law exponent (two-centered solver),
S-potential-localizing isolation construction, route E localized core,
route G as one-clock existence only. Artifacts: attempts/0006/
f1-routing-z2.json, f1-law-measurement.json, f1-derivation.md, result.yaml.

Attempt 0007 CLOSED (2026-08-28): de-boxed two-clock interaction law
MEASURED. Machinery: exact V(S) = -0.5 trS^2 - trS^3 + (trS^2)^2 + 0.5
extracted from the functional, pointwise-validated to 8.5e-13 (an initial
0.5-quartic transcription error was caught by the validation gate before
any law was fitted). V(P) = 0 exactly on the radial-projector family.
Superposing deviations phi_i = S_i - P_i on the common background and
telescoping gives an attractive saturating law Delta E(d) = -A exp(-gamma d),
A ~ 6.3 +- 0.2, gamma ~ 0.25 +- 0.03 (box/refinement stable), no repulsive
core down to d = 6; the confined C-M5S-008 power law is a different regime.
Projector Hessian: NO moduli - physical masses m^2 = 5, 3, 3, 6; gauge tilts
exactly flat in V only. Isolation mechanism sharpened: box-growth = branch
spreading (V-support widens linearly with R on the V zero set; core density
falls 0.0615 -> 0.0313 over R = 20 -> 30) - the committed-class branch has
no fixed-width member; obstruction is branch/EOM structure, not a flat
potential. next_loop 0008: full quadratic-operator spectrum (identify the
gamma mode), fixed-width branch probe, route E. Artifacts:
attempts/0007/v-law-measurement.json, v_law.py, manifest, result.yaml.

Attempt 0008 CLOSED (2026-08-28): GHOST DIAGNOSIS. The full fixed-R
Hessian of the extended fixed-J energy at the R=24 root has 21 negative
directions of 96; lambda_min = -1.9471e4 (FD-verified at three step
sizes) is a PURE chi-row direction: the sigma-truncated boost kinetic
term has the wrong sign in the reduced functional - the attempt-0005
reduction defect is the measured leading instability. Block spectra:
A_pot has exact zero modes (gauge orbits + the marginal W3 structure),
D_curv positive, D_fixedj min -185. Method lesson: zero-VALUED terms at
chi=0 carry nonzero cross-Hessians - block-sum Hessians that omit them
are wrong (caught by FD cross-checks); W3's A-only stability ladder is
an incomplete criterion for the extended class. The candidate-B class
cannot deliver stable one-clock or pair-law claims until the chi sector
is reconstructed from the accepted positive Legendre density.
next_loop 0009: chi-sector kinetic reconstruction (claim-level repair),
S-sector gamma-mode identification, fixed-width probe, route E.
Artifacts: attempts/0008/ghost-diagnosis.json, width-mode.json,
H_*.npy, width_mode.py, result.yaml.

## Debt Ledger
The milestone-0 reproduction inherits the reported verdict's numbers only as transfer references; any gate that fails is recorded as MIXED with the named mechanism rather than silently narrowed. Open frontier (not debt): the full-action scope beyond the tested radial branch, the isolated-existence question for candidates A-C, the two-clock sector, and the interaction law. The absent 0044 source artifacts are external provenance facts, recorded here and in the attempt record, not hidden assumptions inside any proposed claim.

## Review and Promotion Plan
C-M5S-015 receives one independent claim review on the frozen transaction (claim text, reproduced evidence, declared dependencies C-M5S-001..014 scopes, affected consumers) plus at most one correction check. Promotion extracts reusable verification machinery into `src/substrate_framework/`, updates `governance/claims.yaml` and the release manifest, runs `scripts/render_docs.py`, synchronizes accepted-claim memory, and records one content-addressed validation receipt at the frozen boundary.

## Done Gate
P247 closes only when one candidate supplies the isolated finite-energy fixed-J clock, the full declared stability result, a nonsingular independent two-clock sector, and a resolved leading interaction law with claim-level review and promotion. A reproduced de-boxing no-go, a confined replacement, or a conditional force formula is progress evidence that keeps the issue open, and each failing milestone names the construction that unblocks the next attempt.
