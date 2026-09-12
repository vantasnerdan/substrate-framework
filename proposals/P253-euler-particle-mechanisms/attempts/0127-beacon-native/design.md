# 0127-beacon-native - NATIVE back-reaction program, Phase 1 scoping (beacon)

Charter: shepherd NATIVE program (owner declined the R-EM2 import;
0112 draft stands declined-clean per its reversion clause §5: key-and-
lock (0111) stands conditional, charge work BLOCKED, Euler-native
undisturbed — no fault). Phase 1 = SCOPING ONLY. Builds charter
separately. No solves, no new members, no compute in this attempt.

## The wall (inputs, both measured)

- 0124 mechanism: s^6 skirt-curvature + src=0 free-boundary motion
  under ||du||_oo ~ 0.2 makes exact dQ ~ 9 ~ lam_min while linearized
  dQ ~ 4-5. Linear state-error transport is structurally defeated at
  the robust res ~1e-2 floor (3 meshes, P1/P2, 140+ iters, S2 clean).
- 0121 analysis: soft share drained 50% -> 0.2% by fitted mesh, yet
  dF persists 2.08 -> 1.84 — remainder is stiff-mode jf-amplified.
- 0125 asset: conditional lemma with checkable C(e) <= 0.9, trigger
  rho* ~ 0.03 (x7 below floor). UNSATISFIED.

Program thesis: stop pushing state error through the kink. Change
what is transported, or what carries the boundary.

## Route candidates (ranked by mechanism-match)

### R-A tube-uniformity (best match, no transport at all)

Witness lives on SETS: survey min-Q over the stall tube B(u_h, E)
directly with the tested pipeline (corners + axes + random interior,
n3=64/L3=8.0). If min-Q-over-tube >> 0, gap protection holds WITHOUT
any linearization — the mechanism (transport failure) is sidestepped,
not beaten. Cost: assemblies only (~seconds each).
- Threshold/trigger: min eig over >= 2^d corners + 32 interior samples
  >= 5.0 (x2 margin vs charter 0.9-scale; lam_h ~ 11) -> HOLD.
- Stop: min-Q-over-tube < 5.0 with sampling converged (doubling samples
  moves min < 10%) -> route dead, record min-Q map as the obstruction.

### R-B shape-calculus split (attacks the mechanism directly)

Split e = bulk + boundary-motion: explicit src=0 front parametrization
(Hadamard shape derivative for the kink velocity) + smooth-bulk
linear transport (bulk has no kink — 0124 shows bulk linearizes fine:
slope -> G1 clean once crossings stop). First-order exact in the
front displacement; bulk via sharp G1.
- Threshold/trigger: front-displacement norm from du (measured
  crossing count 1 node at full step) vs kink regularization width;
  shape-linear prediction of dQ within 20% of exact at half-scale
  (3.49 printed) -> HOLD the split.
- Stop: shape prediction misses half-scale exact by > 2x -> the
  skirt-curvature (not just front position) dominates; fall to R-A.

### R-C native stiffening (the physics bet — back-reaction proper)

Derive, from the Euler substrate action natively (no EM import), the
two-way coupled operator whose off-diagonal back-reaction blocks add
to Q. Frozen-flux-type argument: axisymmetric toroidal vorticity tied
to fluid motion -> boundary displacement induces a restoring term
from NATIVE dynamics. If the coupled Q_coupled = Q + B with B >= 0
explicit, the wall moves: protection from coercivity, not accuracy.
- Threshold/trigger: B explicit, symmetric, min-eig(B) >= 9.0 on the
  fitted state (covers the 9.1 exact-transfer scale) -> HOLD.
- Stop: no positive-semidefinite native B derivable without importing
  EM structure (owner line) -> route dead; document the obstruction
  as the precise missing construction.

### R-D adjoint DWR (diagnostic only, not a route)

Dual-weighted residual: dQ = <R, z> with adjoint z. Reuses residual
directly, but still linearizes the kink — same wall. Scope: use as a
CROSS-CHECK on R-A/R-B numbers only. Never a build.

## Order and gates

Phase 2 charter ask: R-A first (cheapest, kills-or-holds in one
survey) -> R-B (needs front machinery) -> R-C (needs derivation).
Each route carries its own stop above; any HOLD goes to drift
firewall before downstream use. C(e)/rho* instruments reused
unchanged as the per-state checkers.

## Phase-2 amendments (drift 0127 review, folded pre-charter)

- Expectation sequencing: R-A expected FAST KILL (banked G2' 9.1/3.5
  already spans the tube) — run first anyway (cheapest verdict +
  min-Q map); real action R-B/R-C. Kill-resistance R-B>R-C>R-A.
- R-A tube radius E: inherit G2's E = 2||du||_oo ≈ 0.43 (fitted 0.41,
  trust 0.43). Protection-cost note: smaller E weakens containment
  (u* may lie outside); any smaller-E survey states its weakened
  claim explicitly.
- R-A survey coordinates (F2): span{8 softest J-modes + mu/c bordered
  dirs} (d ≈ 10 collective coords — adversarial by 0121 evidence:
  error concentrates in the soft cluster; G2' corner shows bordered
  dirs matter) + 32 random full-space interior samples. A HOLD on
  anything less licenses nothing.
- R-B HOLD: pre-registered TWO points — half-scale 3.49 ±20% AND
  quarter-scale 1.51 ±20% (both printed). One more assembly defeats
  point-tuning.
- R-B fallback retargeted: miss -> skirt-obstruction record (R-C
  motive), NOT R-A.
- Lesson adopted (correction-review firewall pattern): every printed
  bound in Phase 2 gets a slope- or degeneracy-limit check.
## Non-scope (explicit)

- No member solves, mesh builds, or feed runs in Phase 1 (this file
  is the deliverable).
- No EM import in any form (owner line); R-C derives or dies native.
- No promotion, registry, or G-a2 unblock claimed here.
