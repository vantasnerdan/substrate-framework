# P250 terminal campaign: shell, wall, and bubble sectors of the exterior-degenerate M5 clock (C-M5W-001..008, release v0.173.0)

## Objective and issue

Terminal campaign PR for issue #195: the shell, wall, and bubble sectors of
the exterior-degenerate M5 clock, constructed and promoted on the accepted
C-M5C-001..004 completion. The campaign answers the issue's four questions
with eight accepted claims, exact SymPy oracles, a rigorous interval
certification of the degenerate-crossing frequency, an exact phase
decoupling theorem, the constructed wall/bag family with its thin-wall law,
and the thin-wall stability split. The exact global closure of
omega_c^2 = omega_*^2 (G3 case C) is parked as an independently runnable
future test by owner directive; everything else the issue asks for is
established and promoted.

Canonical issue: #195

Issue relationship: Fixes #195

PR eligibility: terminal campaign success (positive completion receipt
below), extended in place by owner-directed continuation (G1-G4 gap plan).

Campaign terminal evidence: all five obligation nodes S1-S5 of
`proposals/P250-shell-bubble-clock/proposal.yaml` established; the
owner-directed G1 (wall tension) and G2 (bag family) gaps closed and
promoted as C-M5W-006/007; G4 (stability) established analytically at the
reduced-radial scope and promoted as C-M5W-008; G3 advanced to a certified
partial (exact KKT cases A and B closed at omega*^2, minimal polynomial and
isolating enclosure certified; case-C counting parked — see
`attempts/0004/derivation.md`); claims C-M5W-001..008 individually reviewed
and accepted into `governance/claims.yaml`; release `v0.173.0` pinned;
generated docs and accepted-claim memory rendered and valid; in-boundary
debt ledger empty (parked G3 case C is named frontier, not debt).

Active campaign relationship: closes campaign P250.

Authoring agent: Main (executor); v0.172.0 review by distinct agent
ClaimReview (record:
`proposals/P250-shell-bubble-clock/reviews/C-M5W-001-005-review.md`,
correction checks MC-1..MC-4 all closed); v0.173.0 review by three distinct
reviewer agents, one per claim, each returning request_changes with concrete
findings, followed by applied minimum repairs and a one-pass correction
check (Accept, 0.99) — record:
`proposals/P250-shell-bubble-clock/reviews/C-M5W-006-008-review.md`.

Intended merger: repository owner or another distinct merger (no self-merge
direction recorded; the authoring agent leaves the validated PR ready).

- [x] The canonical issue existed before this PR was submitted.
- [x] The merger is distinct, or explicit owner/user self-merge direction is
  recorded. (Distinct merger by default; not self-merged.)

Coordination note: issue #195's owner comment asks contributors not to
comment on issues/discussions until PR review; the executor-protocol comment
that issue #195 requests is therefore recorded here instead: slice = the
whole S1-S5 campaign plus the owner-directed G1-G4 continuation; proposal
P250; claim namespace C-M5W (collision-searched free); branch
`research/P250-shell-bubble-clock`; write surfaces =
proposals/P250-shell-bubble-clock/, src/substrate_framework/m5_wall_clock.py,
tests/test_m5_wall_clock.py, governance (append-only claims + release),
docs/generated, memory, AGENTS.md + two skills (self-improvement).

## Change classification

- [x] Reusable implementation or tooling
- [x] Scientific proposal or campaign evidence
- [x] Claim promotion transaction
- [x] Terminal scientific campaign: positive completion

## Authority and scope

Base release and commit: v0.171.0 at 8b74d3a (exactly the hash frozen in
issue #195); campaign head bd4ace3. Intermediate releases: v0.172.0
(C-M5W-001..005), v0.173.0 (C-M5W-006..008, this extension).

Accepted claims and canonical sources used: C-M5C-001..004 in their exact
accepted scopes; C-M5S scopes not consumed.

Declared imports, assumptions, and conventions: standard linear algebra,
variational calculus, Noether theory, ODE existence/monotonicity arguments,
real-algebraic root isolation (Sturm), declared SymPy methods, mpmath
interval certification, SciPy BVP/collocation with itemized error budgets.
Every threshold derives from the canonical action; the accepted witness
ratio 45/16 is used only as a provenance cross-check and is proven a strict
upper bound (omega_c^2 < 5/3 < 45/16) by exact rational witnesses.

Headline results, one per accepted claim:

- C-M5W-001: exact wall-slice reduction with kinetic metric
  diag(1/4,1/2,1/2,1/2), EL system, and the mechanical first integral
  T - V_omega = const; tension as exact variational functional.
- C-M5W-002: static shells cannot exist (unique vacuum at omega = 0; exact
  Derrick T + 3U = 0 radially) — the stationary shell is a rotating-frame
  object.
- C-M5W-003: the degenerate locus is fixed by the whole S1 orbit; phase
  jumps across it cost exactly zero; two clock regions of arbitrary
  frequencies across one shell have mismatch coefficient of
  (omega_1 - omega_2)^2 exactly zero — area price only (the static shell
  tension), no volume term.
- C-M5W-004: the Maxwell (degenerate-depth) wall exists; deep branch on
  m = 0 with an exact rational Maxwell system; the crossing
  omega_*^2 = 1.663945700059150298856193... is rigorously enclosed to width
  ~1.2e-43 by a two-step Krawczyk iteration with a positive-definite
  fixed-omega Hessian (Gershgorin margins 8.17/4.88/13.13); exact rational
  witnesses give omega_c^2 < 5/3 < 45/16 with omega_c^2 <= omega_*^2
  definitionally. The global identification omega_c^2 = omega_*^2 is the
  named SOS closing route and is not asserted.
- C-M5W-005: fixed-charge bag selection R = 2 sigma/p exactly in the
  thin-wall family, exact envelope closure dE/dQ = omega, interior inertia
  envelope iota_int = -2 dV_min/d(omega^2); the wall-to-infinity limit is
  exactly the P249 exterior-degenerate picture and the Discussion #186
  isofrequency-wall picture is the same wall sector arranged between two
  ticking regions: one mechanism, two arrangements; no M5.32 equivalence.
- C-M5W-006: at the certified Maxwell frequency the slice wall BVP admits
  the monotone kink profile with tension
  sigma_0 = 0.72929841786(58) — the bracket derived from the certified
  budget (total 5.836e-10, b1-dominated), route spread 4.9e-13 across the
  three mechanical routes plus Gauss-Kronrod quadrature, corrected
  boundary-treatment agreement 3.7e-13, headline anchored to a status-0
  solve. (The v0.172.0-era headline 0.7292984178707203 was superseded by
  review: its solver row terminated at status 1, its B-end asymptotic basis
  was evaluated at the wrong bulk point, and its displayed bracket did not
  carry the budget — all repaired with committed evidence; see
  attempts/0002/CORRECTIONS-2026-09-02.md.)
- C-M5W-007: seven stationary S1 wall-bag solutions for
  omega^2 = omega*^2 + delta, delta = 0.001..0.007, mesh-stable radii
  following the exact selection law with chi = 0.999562 -> 0.987601
  (|chi - 1| ~ delta^1.72), envelope dE/dQ = omega to <= 1.9e-4 at the
  physical charge normalization Q = omega I, rung-1 energy matching
  E_crit = (16 pi/3) sigma_0^3/p^2 to 0.16%, and dQ/d omega < 0 across the
  family (Q_true 2.3351e10 -> 6.5152e7). (The archived charge column was
  omega^2-normalized; corrected columns and the retraction of an
  unsupported per-rung acceptance conjunct are recorded in
  attempts/0003/CORRECTIONS-2026-09-02.md and
  attempts/0005/bag_results_corrected.json.)
- C-M5W-008: thin-wall stability split — F''_omega(R_c) = -8 pi sigma < 0
  exactly, so the fixed-omega bag is a Morse-index>=1 saddle of the reduced
  radial problem (critical nucleation bubble, not a minimum), closing the
  radial part of the stability frontier C-M5W-005 named; plus criterion
  satisfaction dQ/d omega < 0 across the certified family and the exact
  envelope identity. No constrained-minimum conclusion is asserted (no
  accepted dependency supplies it); the full spectrum, shape modes, l-lift,
  and thick-wall regime are named frontier (descoped by owner direction).

Owner-directed parking note: the G3 exact-closure certificate
(`attempts/0004/certificate.py`) is committed but not run end-to-end — its
exact number-field machinery certified the KKT case A enumeration, the
degree-32 minimal polynomial of omega*^2, the isolating enclosure, and the
exact emptiness of case B; the case-C root counting, Newton/Krawczyk
solve-back, and the final OMEGAC_CLOSURE verdict remain as a future
independently runnable test. Optimizing that script is not needed for this
PR or its review.

Self-improvement carried in this PR (owner-directed): the campaign's
validation defects — kernel-only derivations, uncommitted results, a
transcribed enclosure denominator, an omega-vs-omega^2 normalization slip,
a headline accepted from a failed solver state — are consolidated into
positive reward language at the decision points: AGENTS.md authoring
discipline (the committed script-and-output pair as the unit of proof;
recompute, don't transcribe), physics-erdos-loop Phase 4 (the replayable
trail), and small-ratio-numerics reproducibility (convergence-status
acceptance, producer/consumer normalization typing, budget-derived
uncertainty brackets).

Files or concerns intentionally out of scope (frontier, not debt): full bag
second-variation spectrum (shape modes, l-lift, thick-wall regime);
fixed-charge constrained-radial second-variation bridge; G3 case-C
certificate completion (parked, independently runnable); planar
vacuum-to-vacuum static loop; non-spherical and multi-shell geometry
corrections; dynamical/radiating shells; two-clock force law; particle
identification; gravity; the P249 audit repair (separate bounded
transaction on the #190 thread).

New or materially changed public interfaces:

| Symbol or module | Authority status | Accepted claim IDs | Owning issue/PR | Conditional boundary |
| --- | --- | --- | --- | --- |
| src/substrate_framework/m5_wall_clock.py | accepted claim | C-M5W-001..008 | #195 | symbolic_verified + numeric_evidence scopes; no particle/two-clock/gravity interpretation |
| exterior_channel_pencil / charged_continuum_edge_squared (pre-existing, untouched) | accepted claim data | C-M5C-002 | #190 audit thread | literal-returning defect class not inherited by P250 oracles; repair owned by the separate audit transaction |

## Useful-unit disposition

| Unit | Local claim or purpose | Dependencies | Evidence and tests | Proposed disposition |
| --- | --- | --- | --- | --- |
| Wall-sector module + tests | exact slice, EL, first integral, locus, slip, Maxwell system, bag laws, resultants (collision-fixed) | C-M5C-001..004 | tests/test_m5_wall_clock.py (16), consumer replay tests/test_m5_exterior_clock.py (9) | merge |
| Claim transaction C-M5W-001..008 | promoted shell/wall/bubble claims | C-M5C scopes | review records + validation receipts | merge |
| Certification artifacts (attempts/0001) | Krawczyk + Gershgorin certificates, route history | C-M5W-004 | certify.py reproduces certificate | merge |
| G1/G2/G3/G4 attempt artifacts (0002-0005) | wall tension, bag family, exact KKT partial, stability split | C-M5W-006..008 | committed scripts + captured outputs; correction notes | merge |

- Artifact merge: yes — novel reusable exact wall-sector API with derived,
  mutation-sensitive oracles.
- Claim promotion: C-M5W-001..008.
- Campaign terminal state: positive completion (G3 case C parked as named
  frontier per owner directive).
- Goal completion: yes — issue #195's S1-S5 nodes established, ten
  achievement gates recorded in the proposal and reviews.

## Validation receipt

```text
# PYTHONPATH=src .venv/bin/python scripts/validate_repository.py
# WORKFLOW VALID: 260 claims, 260 accepted, 12 proposals, 4 skills;
#   MIGRATION QUEUE: 218 units, 0 pending, 0 partial
# bash scripts/validate.sh --pytest-scope tests/test_m5_wall_clock.py
# 16 passed; ALL FIXED REPOSITORY CHECKS AND REQUESTED PYTEST SCOPE PASS
#   (wall-clock suite 15 prior + resultants regression test; consumer
#    replay 9 exterior-clock tests green in the v0.172.0 receipt and
#    untouched by this extension's canonical surface)
# PYTHONPATH=src .venv/bin/python attempts/0005/stability_check.py
#   ALL 18 CHECKS PASS (exact stability split + corrected charge columns)
# PYTHONPATH=src .venv/bin/python attempts/0005/evidence_repair.py
#   ALL 10 CHECKS PASS (status-0 headline, corrected masses/BC/charges)
# proposals/P250-shell-bubble-clock/attempts/0004/certificate.py
#   PARKED (owner directive): partial run certified grad V == system,
#   coercivity, mu deg-32 irreducible, enclosure isolation, case B empty;
#   case C counting + OMEGAC_CLOSURE verdict remain
# git diff --check
# clean
# receipt: proposals/P250-shell-bubble-clock/validation-receipt.yaml
```
