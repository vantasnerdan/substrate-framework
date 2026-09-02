---
description: Constructive review of C-M5W-001..005
author: ClaimReview
created: '2026-09-02'
updated: '2026-09-02'
tags:
- substrate-framework
- claim-review
category: decisions
confidence: working
status: active
---

## Claim and Positive Role

This record reviews the five staged P250 wall-sector claims C-M5W-001..005 as
one claim delta. The claims answer, on the accepted C-M5C-001 canonical action
only: whether the planar stationary wall sector reduces exactly to an aligned
real-psi slice with an exact mechanical first integral (C-M5W-001), why a
static shell cannot exist and why the rotating frame is necessary (C-M5W-002),
whether two clock regions can be decoupled through an orbit-fixed degenerate
locus at exactly zero bulk mismatch (C-M5W-003), whether a stationary wall
exists at an exactly certified degenerate-depth frequency with exact rational
witness bounds (C-M5W-004), and what the thin-wall fixed-charge bag selection,
charge closure, and P249/#186 comparison laws are within the reduced family
(C-M5W-005). None of the five restates a definition or a vacuous special case:
each carries a new exact construction (slice reduction, nonexistence verdict,
zero-mismatch bookkeeping, certified Maxwell point, Laplace/envelope laws) with
stated hypotheses, quantifiers, and exclusions.

## Frozen Transaction

The reviewed transaction is substrate-framework base 8b74d3a, head 5d79e85
(branch research/P250-shell-bubble-clock), two commits by the authoring agent
(96caed1 proposal + attempt 0001 derivation; 5d79e85 wall-sector module, tests,
certification artifacts, proposal update). The claim delta is
C-M5W-001..005 staged in
`proposals/P250-shell-bubble-clock/attempts/0001/claims_block.yaml` (working
tree at review start; not git-tracked at head; the authoring agent confirmed by
IRC that the statements appended to promotion staging are content-identical to
the reviewed snapshot, up to leading-indent normalization). Changed
implementation and evidence records: new
`src/substrate_framework/m5_wall_clock.py` (289 lines, sha256 prefix 5d0ee495a2afb9ed)
and `tests/test_m5_wall_clock.py` (15 tests, sha256 prefix 4a2bf6c7f0f8a05e);
new attempt artifacts `certify.py`, `maxwell_certificate.json`,
`maxwell_point.json`, `derivation.md`, and route scripts `eliminate*.py`
(eliminate.py in 96caed1; eliminate2/3 and the JSON certificates in 5d79e85);
`proposal.yaml` updated in 5d79e85. `src/substrate_framework/m5_exterior_clock.py`
and `tests/test_m5_exterior_clock.py` are untouched by the delta (empty diff
8b74d3a..5d79e85). Accepted dependency propositions actually consumed:
C-M5C-001 (canonical action, orbit, inertia, Legendre envelope), C-M5C-002
(V_M5 >= 0 with equality iff S = N N^T, charged-edge provenance), C-M5C-003
(hylomorphic family for the C-M5W-005 limit statement; 45/16 as provenance
only). Post-head working-tree changes (deletion of the superseded
eliminate*.py scripts and claims_block.yaml, a route-history section added to
derivation.md, release draft governance/releases/v0.172.0.yaml) are outside
this transaction; the route-history verdicts it records were independently
corroborated (see Oracle Audit). Existing validation receipt at review time:
`PYTHONPATH=src .venv/bin/python -m pytest tests/test_m5_wall_clock.py
tests/test_m5_exterior_clock.py -q` — 24 passed (15 wall + 9 exterior) in
4.95 s, re-run by this review.

## Strongest Supported Positive Statement

The strongest result the sourced artifacts support is the five-claim chain at
its stated scopes, with one quantifier reattribution on C-M5W-004. C-M5W-001,
C-M5W-002, C-M5W-003, and C-M5W-005 are supported exactly as written: every
symbolic identity they assert was re-derived independently in this review (see
Evidence Map and Oracle Audit), and their exclusions (no tension value, no
dynamical stability, no thick-wall or stability assertion, no particle or
gravity content) match the evidence. C-M5W-004 is supported at the certified
scope once its certified numeric interval is attributed to the deep-branch
degenerate-crossing frequency omega_*^2 rather than to the global
degenerate-depth omega_c^2 = 2 inf V0/iota: the Krawczyk certificate (re-run
byte-identically) isolates exactly one deep-branch Maxwell solution and
certifies the fixed-omega Hessian positive definite there, the K-box enclosure
omega_*^2 in [1.66394570005915029884, 1.66394570005915029887] is contained in
the claimed bracket [1.6639457000591502, 1.6639457000591504], and the exact
rational witnesses V_w(0, 31/100, 13/20, 41/50) = -13739/18750000 < 0 at
omega^2 = 5/3 and -33854777/25000000 < 0 at omega^2 = 45/16 (both re-verified
exactly here, the second not covered by any test) bound the GLOBAL omega_c^2
strictly: omega_c^2 <= omega_*^2 < 5/3 < 45/16. The identification
omega_c^2 = omega_*^2 (and the parenthetical attainment of the global inf) is
the named global-SOS closing route, exactly as the claim's own first exclusion
states; the statement text should carry that scoping in the interval sentence
too (minimum correction MC-1 below). No claim required narrowing beyond this
reattribution, and no evidence supports a weaker chain.

## Evidence Map

The evidence groups below cover all five claims; each row states what the
group proves, its evidence role, the bridge to the claims, and its limit.

| Evidence | Proposition established | Role | Bridge to claim | Limit |
| --- | --- | --- | --- | --- |
| `m5_wall_clock.py` slice primitives (`clock_slice_tensor`, `wall_slice_potential`, `wall_slice_inertia`, `wall_slice_gradient_density`) | V_w is composed from canonical primitives; reduction identity V_slice == hand form on f >= 0 re-verified exactly in this review against `full_clock_potential`/`clock_inertia_density` | Exact proof (independent rederivation) | C-M5W-001 slice reduction; base for all five claims | Slice sector only; Abs-branch fixed to f >= 0 by assumption, as declared |
| `first_integral_conservation_identity` + tests 95-96 | d/dx(T - V_w) == -sum EL_i u_i', hence T - V_w const on stationary profiles; EL form 2K u'' = dV/du re-verified with K = diag(1/4,1/2,1/2,1/2) | Exact proof | C-M5W-001 first integral and tension chain sigma = 2 int V_w dx = 2 inf int sqrt(V_w) ds_g | Identities exact; the tension VALUE is carried, not evaluated, per exclusion |
| `locus_orbit_fixation`, `orbit_invariance_identity` + tests 99-113 | R(q) diag(m0,c0,c0) R(q)^T fixed for all q (general tangent-isotropic point, not only the vacuum); potential, inertia, velocity-energy, static-gradient densities exactly orbit-invariant | Exact proof | C-M5W-003 locus fixation, zero inertia at locus, single-region exactness | Identities verified on the aligned slice chart; orbit is the declared diagonal S1 |
| `phase_slip_energy_bound` + test 116-121 | Uniform phase-slip cost L^2 dtheta^2 eps / 4 -> 0; formula matches the pinching construction | Exact proof (bookkeeping bound) | C-M5W-003 zero-cost phase jump in H1 | Bound only; no dynamical assertion, per exclusion |
| `static_shell_derrick_residual` + tests 124-130; C-M5C-002 positivity; V_M5(-NN^T) = 2 (re-verified: wall_slice_potential(-1,0,0,0,0) == 2) | Static radial shell excluded: scaling residual T + 3U = 0 with T, U >= 0 forces T = U = 0; no clock-active interior state at omega = 0 | Exact proof (Derrick assembly on canonical normalization) | C-M5W-002 nonexistence and rotating necessity | Radial static shells and the omega = 0 interior requirement only, per exclusion; planar static loop left open as frontier |
| `maxwell_system` + independent derivation (this review): df V_w = 0 gives omega^2 = 32f^2-12f+6-24b exactly; pA, pB, pC == dV/dc, dV/db, 2V_w after substitution; dm V_w = m*(...) factorization | The exact rational deep-branch Maxwell system, derived by differentiate-then-substitute (order matters and is respected) | Exact proof (independent rederivation) | C-M5W-004 system statement | Deep branch m = 0; global branch classification not earned (see route history) |
| `certify.py` + `maxwell_certificate.json`; re-run by this review, byte-identical artifact, strict Krawczyk inclusion True, Gershgorin lower bounds 8.1699/4.8790/13.1268 | Exactly one real solution of the Maxwell system in the rational box (halfwidth 1e-12); fixed-omega Hessian strictly PD over the box; K-box omega enclosure inside the claimed bracket | Exact proof (interval arithmetic, 40 dps) | C-M5W-004 local certification and nondegeneracy | Local to the box; global inf identification NOT certified, per exclusion |
| Exact rational witnesses (tests 189-198, 217-223; this review additionally verified the 45/16 witness -33854777/25000000, which no test covers) | omega_c^2 < 5/3 < 45/16 exactly; the accepted C-M5C-003 ratio is a strict upper bound | Exact proof (global; independent of branch identification) | C-M5W-004 witness bounds | Bounds only; no equality claim |
| `envelope_identity`, `inertia_envelope_identity`, `thin_wall_bag_radius` + tests 226-233 | dE~/dQ = omega for any smooth stationary family (general symbolic proof verified); chain-rule decomposition total = partial + chain verified on a generic instance; R = 2 sigma/p | Exact proof (envelope) + corroborating bookkeeping (inertia chain instance) | C-M5W-005 closure and selection law | Chain-rule instance does not itself exercise stationarity-kills-chain (one-line generic calculus); thin-wall family scope per exclusions |
| `derivation.md` route history (post-head edit, corroborated by this review: MW = Res_f(Res_c(pA,pB), Res_c(pA,pC)) has degree 108 in omega^2, reproduced in 29.4 s) | Route A Groebner blocked, Route B computed-but-prohibitive (global enumeration NOT EARNED), Route C PSLQ falsified, Route D adopted | Provenance / route honesty | C-M5W-004 exclusion 1 credibility | No claim rests on global root enumeration; artifact `elimination.json` was never produced |
| `tests/test_m5_wall_clock.py` mutation gates (140-176) + module docstring invariant | Lock-strength mutation (6 -> 6+eps) moves the derived pB away from the lock manifold while invisible on it; scalar mutation W + eps f^2 shifts the f-relation by exactly 2 eps, invisible at f = 0 | Regression / false-green exclusion (#190 class) | All claims' oracle integrity | Two channels covered (lock, scalar); structural identities end in symbolic == 0 and are not literal-returning |
| `tests/test_m5_exterior_clock.py` (9 tests, untouched, green) | Consumer replay: canonical primitives unchanged and consistent | Regression | Compatibility axis for all claims | Replay only; no new consumer obligations |

## Oracle Audit

The strongest practical oracle for the load-bearing proposition (C-M5W-004's
certified Maxwell point) is `certify.py`, and this review re-ran it rather than
trusting the artifact. The analytic-closure receipt precedes any numerics: the
Maxwell system is exact rational polynomial algebra derived by differentiation
from the canonical potential (differentiate first, substitute omega^2 after —
the module comment at `m5_wall_clock.py:201-203` correctly identifies the
ordering hazard, and the independent rederivation confirms the result), and the
only numerical step is rigorous interval arithmetic (mpmath iv, 40 dps). The
Krawczyk formula in `certify.py:53-65` is the standard operator
K = y - Y f(y) + (I - Y J(X))(X - y) with Y = J(y)^{-1} the midpoint inverse,
and strict inclusion K ⊂ int X implies exactly one root of the system in X —
the standard Krawczyk theorem, correctly applied with J(X) evaluated as an
interval matrix over the whole box. The computed predicate answers exactly the
stated remainder (local isolation + nondegeneracy of the deep-branch
degenerate-depth point), not a different question. The Gershgorin step is valid
for the FIXED-omega Hessian: H = d²V_omega/dc²,db²,df² with V_omega =
V0 - (omega²/2)(f²+4b²) has entries affine in omega² (verified structurally in
this review), and the omega² enclosure W = [32 f_lo² - 12 f_lo + 6 - 24 b_hi,
32 f_hi² - 12 f_hi + 6 - 24 b_lo] is a valid monotonicity-based enclosure over
the box (64f - 12 > 0 and -24 < 0 on the box), so interval evaluation of each
entry is a valid enclosure and row-wise min-diagonal minus max-offdiagonal-sum
> 0 certifies every matrix in the enclosure is positive definite. Re-run
result: strict inclusion True, bounds 8.1699/4.8790/13.1268, artifact
reproduced byte-identically (diff clean), runtime 0.5 s. Mutation coverage
verifies the gates are not literal-returning: both channel-local mutations
(lock strength, scalar W) demonstrably move the derived system while staying
invisible at the witness class, and no test compares an oracle output to a
hardcoded literal (the boolean identity gates end in symbolic simplification to
zero; the artifact-checking tests have regression role only, with the live
oracle re-runnable and re-run). Computed predicates elsewhere answer their
stated remainders: the Derrick residual solves the static-shell scaling
condition T + 3U = 0 (correct 3D expansion scaling, re-checked), the phase-slip
bound formula matches the pinching construction exactly
(2 eps · (L eps)² (dtheta/2eps)²/2 = L² dtheta² eps/4), and the ψ_I slice
invariance claimed in C-M5W-001 was verified here directly (EL residual for the
imaginary scalar part vanishes identically at p ≡ 0; potential stiffness in p
equals 2(12 b f + 16 f³ - 6 f² + 3) > 0 for f > 0, minimum 2.875 at f = 1/4
before the nonnegative 12 b f term). The exterior-vacuum nondegeneracy consumed
from C-M5C-002 was spot-verified: all eight generalized vacuum mass squares are
strictly positive (raw-entry chart {6,6,8,8,10,10,22,44}; the accepted chart
multiset is {4,4,6,6,10,10,22,22} — see finding F4 on the derivation receipt's
normalization wording; positivity holds in every chart). No pre-gate
exploratory sample entered any claim: the numeric landscape scan in
derivation.md is marked exploratory, and no claim text cites it.

## Findings

No current blocker exists: no counterexample or contradiction was found under
the stated hypotheses, no load-bearing step is absent or circular, every
declared dependency supplies what is used, and the consumer replay passes. The
findings below are one minimum correction (boundary of the certified scope) and
follow-ups.

| Finding | Direct evidence | Blocking in boundary / minimum correction / follow-up | What would resolve or overturn it |
| --- | --- | --- | --- |
| F1: C-M5W-004 attributes the certified interval to omega_c^2, but the certified object is the deep-branch degenerate-crossing frequency omega_*^2; omega_c^2 <= omega_*^2 with equality = the global-SOS closing route the claim itself declines to assert. The parenthetical "(inf attained, by coercivity of the ratio and the exact positivity of the excluded directions)" is part of the same unasserted global identification. | claims_block C-M5W-004 statement ("giving omega_c^2 in [1.6639457000591502, 1.6639457000591504]") vs its own exclusion 1 ("the identification of that branch's crossing as the global inf ... not asserted"); certify.py K boxes; witnesses (global) | Minimum correction (wording/quantifier; not refutation — the local certificate, nondegeneracy, and all witness inequalities hold under either reading) | Rewrite the interval sentence with omega_*^2 (deep-branch crossing), state omega_c^2 <= omega_*^2 and omega_c^2 < 5/3 < 45/16 (the latter global and exact via the witnesses), and move "inf attained" into the named-SOS closing route. Upgrade to the original wording: complete the global SOS certificate (or an equivalent exact global comparison of all real branches), which would also earn the route-B/B-side global enumeration currently recorded as not earned |
| F2: the recorded certificate backs a wider enclosure than the claim quotes: `maxwell_certificate.json` records w(X) = [1.6639456999950312, 1.6639457001232694], while the narrow claimed bracket follows only from the K boxes that certify.py prints but does not record. Separately, `maxwell_point.json` records 20-digit coordinates while its residual_max 4.81e-50 refers to the untruncated 50-dps solve; the recorded coordinates give residual ~5.0e-19 (verified). | `certify.py:64-67` (K boxes printed, not persisted); `maxwell_certificate.json` w_interval; `maxwell_point.json`; this review's re-run output and 60-dps residual check | Follow-up (artifact completeness; should land with the promotion since the JSON is release evidence) | Record the K-box (or K-box-derived omega) enclosure in maxwell_certificate.json; record >= 50-digit coordinates in maxwell_point.json or restate the residual at the recorded precision. Overturn check: none needed — the narrow interval is a true consequence of the certified K boxes (verified: K-box omega enclosure = [1.66394570005915029884, 1.66394570005915029887] ⊂ claimed bracket) |
| F3: C-M5W-001's leading quantifier ("every planar stationary wall ... is carried, up to the exact orbit gauge, by the aligned real-psi slice") rests on exact verified ingredients (ψ_I-equation invariance verified here; shear stiffness 8x²-3x+1 strictly positive per variable, discriminant -23; orientation stiffness 48 b f² > 0 off the locus; orbit invariance) plus the standard elliptic-uniqueness assembly (homogeneous equations with positive stiffness and slice/vacuum boundary data), which derivation.md states but does not discharge as an explicit exact step. | `m5_wall_clock.py` identities; tests 99-121; derivation.md items 1, 5; this review's ψ_I and stiffness checks | Evidence-role qualification (accepted with the role noted: exact ingredients + standard assembly) | Record the shear/orientation/ψ_I homogeneous-ODE uniqueness argument as an explicit exact derivation entry (one subsection); that upgrades the quantifier to fully discharged. Overturn would require a planar stationary wall escaping the slice despite positive channel stiffness — excluded by the verified stiffness signs |
| F4: derivation.md item 4 says the re-derived exterior pencil "matches the accepted multiset", listing {10,22,10,4,4,6}; C-M5C-002's accepted multiset is {4,4,6,6,10,10,22,22} (its chart) and the raw-entry chart gives {6,6,8,8,10,10,22,44}. The listing matches neither as printed. Only positivity is consumed by this delta and holds in all charts. | derivation.md item 4; this review's generalized-eigenvalue computation; governance/claims.yaml C-M5C-002 | Follow-up (receipt wording; C-M5C-002 itself is an unchanged accepted dependency, outside this review) | State the chart normalization or drop the multiset-equality sentence in favor of citing C-M5C-002 directly. No upgrade beyond wording; no claim consumes the multiset |
| F5: `inertia_envelope_identity` verifies the chain-rule decomposition (total = partial + chain) on a generic instance with arbitrary u(omega²); the stationarity-kills-chain step is one line of generic calculus not separately exercised by the instance. The companion `envelope_identity` (dE~/dQ = omega) IS a general symbolic proof. | `m5_wall_clock.py:265-289`; tests 226-228 | Evidence-role note (corroborating bookkeeping, not a gap — the claim text does not overstate it) | None required; optionally add a stationary-u instance so the test exercises the full iota_int = -2 dV_min/domega² chain |
| F6: orientation-Hessian convention: the claim states 48 b f²; a half-angle parameterization gives 12 b f² for the same physical quantity. | claims_block C-M5W-001; post-head derivation.md annotation | No defect (already annotated post-head) | None |

## Compatibility and Consumers

The transaction is a pure compatible extension: it adds one module and one test
file, touches no existing source, and the only affected consumers are the
exterior-clock tests, which replay green unchanged (9/9 within the 24-pass
run). Declared hypotheses bind as stated: the constrained auxiliary frame and
lock strengths are C-M5C-001 model data; fields are H1 where the phase-ramp
family needs it; the |psi| branch is fixed to the nonnegative slice; W(f) and
the kinetic normalization are canonical. Units and conventions are inherited
unchanged from C-M5C-001..004; no new constants are imported from any other
model (spot-checked: every constant in the claims traces to
`wall_slice_potential` composition — witnesses, omega relations, and stiffness
quadratics all re-derived from the canonical potential in this review).
Dependency closure verified proposition by proposition: C-M5W-001 uses only
C-M5C-001 (action, orbit, inertia density, envelope); C-M5W-002 uses
C-M5W-001 (slice V0 and gradient normalization) plus C-M5C-002 (V_M5 >= 0 with
equality iff N N^T, unique vacuum) — both supplied; C-M5W-003 uses C-M5W-001
(density layers, orbit invariance, inertia-free locus) plus C-M5W-004
(sigma_0 as the decoupling price) — supplied; C-M5W-004 uses C-M5W-001 (V_omega,
first integral, excluded-direction stiffness), C-M5W-002 (static nonexistence
motivating the rotating frame), and C-M5C-002 (positivity; edge 4 as provenance
only) — supplied, and the C-M5C-003 witness is re-derived fresh rather than
imported (no undeclared dependency); C-M5W-005 uses C-M5W-001 (V_omega,
inertia), C-M5W-003 (locus mechanism and the isofrequency reading),
C-M5W-004 (sigma(omega), omega_c, p(omega)), and C-M5C-003 (hylomorphic family
for the wall-to-infinity realization; dE/dQ = omega reaches it transitively via
C-M5W-001's C-M5C-001 dependency) — supplied. Declared hypotheses and the
open parent work (global SOS identification, thick-wall regime, bag stability,
planar static loop, curvature corrections, bundle charts) are frontier, not
debt: each is named in a claim exclusion or the proposal route frontier, and
none is consumed inside this transaction. No defect was introduced inside the
delta.

## Four-Axis Decision

Decisions per claim, without inflating one axis from another.

- C-M5W-001 (slice reduction and mechanical first integral):
  - Verification: verified — every asserted identity re-run or independently
    re-derived (reduction, EL residuals, first integral, orbit invariance,
    stiffness signs, ψ_I invariance, Derrick-supporting normalizations).
  - Review: accepted, with the F3 evidence-role note on the "every wall"
    quantifier (exact ingredients verified; elliptic-uniqueness assembly
    stated, not yet discharged as an explicit exact step).
  - Compatibility: compatible_extension — additive module/tests, consumers
    untouched and green.
  - Epistemic: active.
  - Relationship: depends on C-M5C-001 only; supplied.
  - Strongest accepted statement: as written, with the reduction quantifier
    carried by exact ingredients plus the standard assembly until F3's
    explicit entry lands.
- C-M5W-002 (static-shell nonexistence and rotating necessity):
  - Verification: verified — unique-zero decomposition from C-M5C-002 plus the
    canonical summands; V_M5(-N N^T) = 2 and inertia-free vacuum re-verified
    exactly; Derrick residual T + 3U = 0 correct for the 3D expansion scaling.
  - Review: accepted.
  - Compatibility: compatible_extension.
  - Epistemic: active.
  - Relationship: C-M5W-001 + C-M5C-002; both supplied.
  - Strongest accepted statement: as written, with its declared exclusions
    (radial static shells and the omega = 0 interior requirement only; planar
    static loop open frontier).
- C-M5W-003 (orbit-fixed locus and exact phase decoupling):
  - Verification: verified — locus fixation (general tangent-isotropic point),
    all four density layers orbit-invariant, inertia-free locus, uniform
    phase-slip bound, exactly zero (omega_1 - omega_2)² coefficient; the
    gluing statement accepted as the energy bookkeeping theorem its exclusion
    declares (no dynamical stability or wall-force assertion made).
  - Review: accepted.
  - Compatibility: compatible_extension.
  - Epistemic: active.
  - Relationship: C-M5W-001 + C-M5W-004; both supplied.
  - Strongest accepted statement: as written within the glued-configuration
    family scope.
- C-M5W-004 (Maxwell wall existence and frequency characterization):
  - Verification: verified at the certified scope — Krawczyk uniqueness and
    fixed-omega Gershgorin PD re-run byte-identically; the K-box enclosure
    lies inside the claimed bracket; both exact rational witnesses re-verified
    (including the test-uncovered 45/16 witness); the omega_c^2 interval
    attribution and the parenthetical global attainment are the one scope
    overreach (F1).
  - Review: accepted_with_minimum_corrections (MC-1 = F1 wording; MC-2 = F2
    artifact backfill).
  - Compatibility: compatible_extension.
  - Epistemic: active.
  - Relationship: C-M5W-001 + C-M5W-002 + C-M5C-002; all supplied; no
    undeclared import (witness re-derived fresh).
  - Strongest proposed statement: the corrected form — exactly one deep-branch
    Maxwell solution in the certified box, nondegenerate degenerate-depth
    crossing there with omega_*^2 in the claimed bracket, omega_c^2 <=
    omega_*^2, omega_c^2 < 5/3 < 45/16 exactly (global witnesses), global
    identification = named SOS closing route.
- C-M5W-005 (thin-wall bag selection, charge closure, comparison):
  - Verification: verified within the reduced thin-wall family — R = 2 sigma/p,
    dE~/dQ = omega (general symbolic proof), iota_int = -2 dV_min/domega²
    (chain rule with the F5 role note), limits to the P249 picture and the #186
    arrangement reading; thin-wall/thick-wall and stability scoping honest.
  - Review: accepted.
  - Compatibility: compatible_extension.
  - Epistemic: active.
  - Relationship: C-M5W-001 + C-M5W-003 + C-M5W-004 + C-M5C-003; all supplied.
  - Strongest accepted statement: as written, with "no M5.32 equivalence" and
    no particle/force/gravity content, per exclusions.

Merge readiness: READY TO MERGE conditional on applying MC-1 (C-M5W-004
wording reattribution, mandatory before accepted_in/review fields are set in
governance/claims.yaml and before governance/releases/v0.172.0.yaml is
finalized with its already-staged C-M5W-001..005 entries) and MC-2 (certificate
artifact backfill, to land with the promotion). No blocker exists; with MC-1/MC-2
applied, all five claims are acceptance-ready at verification symbolic_verified,
compatibility compatible_extension, epistemic active, review accepted (C-M5W-004
accepted after correction).

## Promotion Transaction

The promotion requires exactly these changes and no others: append the five
claims to governance/claims.yaml with the MC-1-corrected C-M5W-004 statement
and review/accepted_in fields set from this record's decisions; finalize
governance/releases/v0.172.0.yaml (the staged C-M5W-001..005 entries must
reference the corrected statements); land MC-2 in
`maxwell_certificate.json`/`maxwell_point.json` before release pinning since
they are release evidence; update the memory contract page
(memory/vantasner/proposals/P250-shell-bubble-clock.md) to the promoted state;
no implementation, test, or generated-doc changes are required beyond this —
module and tests are frozen and green, and the existing content-addressed
validation receipt (24-pass run, re-executed by this review) remains valid
because the code and inputs are unchanged. The license registry entries
L-S1-WALL through L-S4-COMPARISON in proposal.yaml remain "unearned" in the
frozen proposal; if the campaign convention is to mark them earned at promotion,
that flip belongs to the authoring agent's promotion commit, informed by this
record, and is not part of the claim review itself.

## Correction Check

Executed at head 84a8db4 (MC-1/MC-2 promotion commit), covering exactly the
corrected C-M5W-004 statement, the two JSON artifacts, certify.py, the registry
and release edges, and the unchanged sibling statements. Outcome: MC-1 landed
correctly in governance/claims.yaml (interval reattributed to omega_*^2,
omega_c^2 <= omega_*^2 definitional, witnesses global, inf-attainment inside
the named SOS route; C-M5W-005's two terminology adaptations are consistent
downstream renamings; C-M5W-001/002/003 content-faithful; 24/24 tests green at
the new head). MC-2 introduced a new defect, found by this check:

- MC-3 (required follow-up correction, non-blocking to claim content): the new
  `k_omega_hi` in certify.py pairs the wrong corners — it takes
  max(w(f_a,b_a), w(f_b,b_b)) instead of the max corner w(f_b,b_a) — so the
  recorded `krawczyk_omega_enclosure` upper endpoint 1.66394570005915029885691
  understates the true max of w over the persisted K boxes
  (1.6639457000591502988681629270496958798642644584868...; gap exactly
  24 x (b-box width) ~= 1.13e-20). The recorded enclosure therefore does not
  follow from the Krawczyk operator output it accompanies, and the C-M5W-004
  statement's bracket upper endpoint (1.6639457000591502988569) inherits it.
  The bracket remains numerically true of the crossing (true value
  ~1.663945700059150298856193..., inside the bracket), so nothing is refuted;
  the statement-artifact-certificate triangle is just inconsistent at the
  20th decimal. Fix: correct the corner pairing (use K[2].b with K[1].a),
  regenerate the certificate, and widen the statement bracket to the certified
  upper endpoint — or produce a rigorous tighter bracket first (one further
  Krawczyk iteration on the K box, cheap at 40 dps, would certify far tighter
  than either value). Also fix the related provenance slip: at the new
  maxwell_point.json's stored 55-digit coordinates the residual evaluates to
  2.912e-54 (this check, 80 dps), not the recorded 2.18e-60, so the note
  "residual evaluated at stored precision" is inaccurate; record the
  stored-precision residual or the full-precision coordinates, and align
  precision_dps (still 50) with the 55 stored digits.

No other corrected statement, evidence role, dependency edge, or consumer edge
required attention. Not needed for any other claim.

Closing verification of MC-3 at head 99ea88f (requested by the author; scope
limited to the corrected certification chain): the two-step Krawczyk iteration
is rigorous — the operator contains every root of its input box, strict
inclusion holds at both steps (K2 ⊂ K, reproduced), and the omega enclosure is
a direct interval evaluation of w over K2 (no corner pairing); the recorded
enclosure [1.66394570005915029885619300029616144415736458,
1.66394570005915029885619300029616144428226086] contains this review's
independent 120-dps monotone corner evaluation over the persisted K2 boxes and
the 55-digit point's crossing value; artifact reproduces byte-identically;
24/24 tests green; maxwell_point.json now records residual_max 2.912e-54
(matching this review's 80-dps evaluation of the stored coordinates exactly),
precision_dps 55, inaccurate note removed. MC-3 CLOSED. One residual
micro-defect introduced by the closure rounding:

- MC-4 (one-line micro-correction, same statement-not-implied-by-certificate
  class as MC-1/MC-3 at the 39th significant digit): the C-M5W-004 statement
  bracket endpoints are rounded INWARD relative to the recorded rigorous
  enclosure — statement lower 1.663945700059150298856193000296161444158 >
  recorded lower ...44415736458, and statement upper
  ...444282 < recorded upper ...44428226086. A valid bracket must round the
  recorded endpoints outward (e.g. lower ...444157, upper ...444283 at that
  digit count) or quote the enclosure endpoints in full. Magnitude ~2.6e-40,
  at the noise floor of the 40-dps interval arithmetic; claim content
  unaffected.

## Result and Frontier

The positive result retained is the full five-claim wall-sector chain at its
certified scopes: an exact slice reduction with mechanical first integral
(C-M5W-001), an exact static-nonexistence verdict forcing the rotating frame
(C-M5W-002), exactly zero bulk mismatch decoupling through the orbit-fixed
locus (C-M5W-003), a locally certified nondegenerate degenerate-crossing Maxwell
point with exact global witness bounds 5/3 < 45/16 and honest global-SOS
scoping (C-M5W-004, corrected wording), and exact thin-wall selection, charge
closure, and comparison laws with explicit frontier boundaries (C-M5W-005). The
one decisive open construction for the strongest form of C-M5W-004 is the
global SOS certificate identifying the deep-branch crossing with the global inf
of V0/iota (route history records Groebner blocked and degree-108 resultant
enumeration computationally prohibitive; the Krawczyk route is the adopted
local replacement). Adjacent items stay frontier and do not expand this
transaction: thick-wall regime and bag linear stability, planar static
homoclinic loops, geometry-ladder curvature corrections, bundle/director charts,
dynamical shells, and any particle, force, or gravitational content.

## Cross-References

Proposal: proposals/P250-shell-bubble-clock/proposal.yaml (frozen at 96caed1,
updated 5d79e85). Claims: attempts/0001/claims_block.yaml (working-tree
snapshot reviewed; content-identical promotion staging confirmed by author).
Dependencies: governance/claims.yaml C-M5C-001..004 (accepted v0.171.0, outside
this review). Evidence: attempts/0001/{derivation.md (head + post-head route
history), certify.py, maxwell_certificate.json, maxwell_point.json,
eliminate.py, eliminate2.py, eliminate3.py}. Implementation and tests:
src/substrate_framework/m5_wall_clock.py, tests/test_m5_wall_clock.py;
consumer replay tests/test_m5_exterior_clock.py. Parent arc: campaign
P249-exterior-degenerate-clock (adjudication, attempts 0001-0008, C-M5C
review), issues #195, #190, #186, release v0.171.0. Promotion: head 84a8db4
(MC-1/MC-2 corrections, registry acceptance, release v0.172.0). Reviewer
verification scripts were scratch-only (/tmp) and are not part of the record.

## Verdict

Per-claim four-axis decisions are recorded in the Four-Axis Decision section:
C-M5W-001 accepted (verified / accepted / compatible_extension / active, F3
role note), C-M5W-002 accepted (verified / accepted / compatible_extension /
active), C-M5W-003 accepted (verified / accepted / compatible_extension /
active), C-M5W-004 accepted with minimum corrections (verified at certified
scope / accepted_with_minimum_corrections / compatible_extension / active),
C-M5W-005 accepted (verified / accepted / compatible_extension / active).
Overall merge readiness: all requested corrections have landed and been
verified — MC-1 (statement reattribution) applied and verified at 84a8db4;
MC-3 (rigorous two-step Krawczyk enclosure replacing the defective corner
pairing) applied and independently verified, and closed, at 99ea88f, which
also reconciled the maxwell_point.json residual/precision metadata. The single
outstanding item is MC-4, a one-line endpoint-rounding micro-correction in the
C-M5W-004 statement bracket (round the quoted enclosure endpoints outward, or
quote them in full); claim content is unaffected. With MC-4 applied, the
merged state is final and acceptance-complete for all five claims.
Corrections ledger: MC-1 — statement reattribution — APPLIED AND VERIFIED
(upgrade to stronger wording: complete the global SOS certificate). MC-2 —
artifact backfill — APPLIED WITH DEFECT MC-3, SUPERSEDED. MC-3 — rigorous
two-step Krawczyk enclosure — APPLIED, INDEPENDENTLY VERIFIED, AND CLOSED at
99ea88f (residual micro-defect MC-4: statement bracket endpoints rounded
inward vs the recorded enclosure at the 39th significant digit; fix by
rounding the quoted endpoints outward or quoting them in full). Correction
checks: one substantive check at 84a8db4 plus the requested closing
verification of MC-3 at 99ea88f, both recorded above; this review pass is
complete.
