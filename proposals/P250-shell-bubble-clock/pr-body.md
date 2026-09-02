## Objective and issue

Terminal campaign PR for issue #195: the shell, wall, and bubble sectors of
the exterior-degenerate M5 clock, constructed and promoted on the accepted
C-M5C-001..004 completion. The campaign answers the issue's four questions
with five accepted claims, exact SymPy oracles, a rigorous interval
certification of the degenerate-crossing frequency, and an exact phase
decoupling theorem.

Canonical issue: #195

Issue relationship: Fixes #195

PR eligibility: terminal campaign success (positive completion receipt below).

Campaign terminal evidence: all five obligation nodes S1-S5 of
`proposals/P250-shell-bubble-clock/proposal.yaml` established; claims
C-M5W-001..005 individually reviewed and accepted into
`governance/claims.yaml`; release `v0.172.0` pinned; generated docs and
accepted-claim memory rendered and valid; in-boundary debt ledger empty.

Active campaign relationship: closes campaign P250.

Authoring agent: Main (executor), review by distinct agent ClaimReview
(record: `proposals/P250-shell-bubble-clock/reviews/C-M5W-001-005-review.md`,
correction checks MC-1..MC-4 all closed).

Intended merger: repository owner or another distinct merger (no self-merge
direction recorded; the authoring agent leaves the validated PR ready).

- [x] The canonical issue existed before this PR was submitted.
- [x] The merger is distinct, or explicit owner/user self-merge direction is
  recorded. (Distinct merger by default; not self-merged.)

Coordination note: issue #195's owner comment asks contributors not to comment
on issues/discussions until PR review; the executor-protocol comment that
issue #195 requests is therefore recorded here instead: slice = the whole
S1-S5 campaign; proposal P250; claim namespace C-M5W (collision-searched
free); branch `research/P250-shell-bubble-clock`; write surfaces =
proposals/P250-shell-bubble-clock/, src/substrate_framework/m5_wall_clock.py,
tests/test_m5_wall_clock.py, governance (append-only claims + release),
docs/generated, memory.

## Change classification

- [x] Reusable implementation or tooling
- [x] Scientific proposal or campaign evidence
- [x] Claim promotion transaction
- [x] Terminal scientific campaign: positive completion

## Authority and scope

Base release and commit: v0.171.0 at 8b74d3a (exactly the hash frozen in
issue #195); campaign head 0455f88.

Accepted claims and canonical sources used: C-M5C-001..004 in their exact
accepted scopes; C-M5S scopes not consumed.

Declared imports, assumptions, and conventions: standard linear algebra,
variational calculus, Noether theory, ODE monotonicity arguments, SymPy
exact arithmetic, mpmath interval certification. Every threshold derives
from the canonical action; the accepted witness ratio 45/16 is used only as
a provenance cross-check and is proven a strict upper bound
(omega_c^2 < 5/3 < 45/16) by exact rational witnesses.

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

Files or concerns intentionally out of scope (frontier, not debt): bag
linear stability; thick-wall regime near the charged edge; non-spherical
and multi-shell geometry corrections; dynamical/radiating shells; global
SOS certificate for omega_c^2 = omega_*^2; two-clock force law; particle
identification; gravity; the P249 audit repair (separate bounded
transaction on the #190 thread).

New or materially changed public interfaces:

| Symbol or module | Authority status | Accepted claim IDs | Owning issue/PR | Conditional boundary |
| --- | --- | --- | --- | --- |
| src/substrate_framework/m5_wall_clock.py | accepted claim | C-M5W-001..005 | #195 | symbolic_verified scope; no particle/two-clock/gravity interpretation |
| exterior_channel_pencil / charged_continuum_edge_squared (pre-existing, untouched) | accepted claim data | C-M5C-002 | #190 audit thread | literal-returning defect class not inherited by P250 oracles; repair owned by the separate audit transaction |

## Useful-unit disposition

| Unit | Local claim or purpose | Dependencies | Evidence and tests | Proposed disposition |
| --- | --- | --- | --- | --- |
| Wall-sector module + tests | exact slice, EL, first integral, locus, slip, Maxwell system, bag laws | C-M5C-001..004 | tests/test_m5_wall_clock.py (15), consumer replay tests/test_m5_exterior_clock.py (9) | merge |
| Claim transaction C-M5W-001..005 | promoted shell/wall/bubble claims | C-M5C scopes | review record + validation receipt | merge |
| Certification artifacts (attempts/0001) | Krawczyk + Gershgorin certificates, route history | C-M5W-004 | certify.py reproduces certificate | merge |

- Artifact merge: yes — novel reusable exact wall-sector API with derived,
  mutation-sensitive oracles.
- Claim promotion: C-M5W-001, C-M5W-002, C-M5W-003, C-M5W-004, C-M5W-005.
- Campaign terminal state: positive completion.
- Goal completion: yes — issue #195's S1-S5 nodes established, ten
  achievement gates recorded in the proposal and review.

## Validation receipt

```text
# PYTHONPATH=src .venv/bin/python scripts/validate_repository.py
# WORKFLOW VALID: 257 claims, 257 accepted, 12 proposals, 4 skills;
#   MIGRATION QUEUE: 218 units, 0 pending, 0 partial
# bash scripts/validate.sh --pytest-scope tests/test_m5_wall_clock.py tests/test_m5_exterior_clock.py
# 24 passed; ALL FIXED REPOSITORY CHECKS AND REQUESTED PYTEST SCOPE PASS
#   (run at the MC-3 boundary 99ea88f and re-verified after MC-4: the MC-4
#    diff is claim-statement wording + generated renders, which does not
#    stale a scientific oracle; pytest scope and registry validation were
#    re-run after MC-4)
# proposals/P250-shell-bubble-clock/attempts/0001/certify.py
#   Krawczyk strict inclusion (two steps); Gershgorin 8.17/4.88/13.13;
#   omega_*^2 enclosure width ~1.2e-43
# git diff --check
# clean
# receipt: proposals/P250-shell-bubble-clock/validation-receipt.yaml
```
