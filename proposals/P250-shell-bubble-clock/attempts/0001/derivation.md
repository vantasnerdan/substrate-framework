# Attempt 0001 — wall-sector reduction, Maxwell system, decoupling and bag laws

Active obligation: S1 (feeds S2/S3/S4 directly).
Route: WALL_ROTATING_MAXWELL (primary), WALL_STATIC_PROFILE (competing, refuted), SLIP_LOCUS_DECOUPLING.

## Analytic specification receipt (exact layer, no numerics consumed)

Object: planar wall of the C-M5C-001 action, adapted real-psi aligned slice
S = diag(m, c+b, c-b), psi = f >= 0. Kinetic metric diag(1/4, 1/2, 1/2, 1/2)
in (m, c, b, f). Exact slice potential (verified == full_clock_potential on
f >= 0):

  V0 = V_M5(m,c,b) + 2c^2 + 2b^2 + 6(b - f^2)^2 + W(f),
  V_M5 = -r^2/2 - (m^3 + 2c^3 + 6cb^2) + r^4 + 1/2,  r^2 = m^2 + 2c^2 + 2b^2.

Rotating frame: V_w = V0 - (w^2/2) * iota, iota = f^2 + 4b^2 (canonical
inertia primitive, exact on the slice).

## Verified exact facts (SymPy, from canonical primitives)

1. Slice reduction: V_slice == hand form on f >= 0 (Abs(f) branch resolved);
   inertia == f^2 + 4b^2; canonical velocity energy Tr(Sdot^2)/4 + |psidot|^2/2.
2. EL system + mechanical first integral: d/dx(T - V_w) = -sum EL_i u_i'
   (identity verified); T - V_w = const on stationary profiles.
3. Static nonexistence: radial static shell killed by exact Derrick
   (T + 3U = 0 with T,U >= 0); interior state at w = 0 does not exist
   (V0 >= bound + axislock + phaselock + W >= 0, unique zero at S = NN^T,
   psi = 0; V_M5(-NN^T) = 2 > 0 exact). Planar static vacuum-to-vacuum
   homoclinic loops: not excluded by the first integral alone; recorded as
   representation frontier, not needed for the objective.
4. Exterior pencil re-derived from the canonical potential (closing the #190
   literal defect for P250 oracles): stiffness diag(5,22,10,4,4,6) in
   (a, split, trace, u, v, psi_Im), kinetic diag(1/2,1,1,1,1,1), masses
   {10,22,10,4,4,6}, scalars 6 — matches the accepted multiset.
5. Excluded directions along clock-active branches (exact):
   shear Hessians 2(8b^2-3b+8c^2-3c+4m^2-3m+1) (+ cross 8v^2/8u^2),
   orientation Hessian 48 b f^2 (true rotation angle; positive on b,f>0,
   vanishing exactly at the locus = microscopic perfect slip).
6. Orbit-fixed locus: R(q) diag(m0,c0,c0) R(q)^T = diag(m0,c0,c0) exactly;
   inertia density 0 there; potential density, inertia density, velocity
   energy, and static gradient density are all EXACTLY orbit-invariant
   (symbolic identities, trigsimp-reduced to 0).
7. Phase slip calculus: a phase jump of any size at the locus costs
   <= L_f^2 (delta theta)^2 eps / 4 -> 0 (f(x*) = 0); tensor orientation
   jump is exactly free by orbit-fixation. Hence S2 mismatch coefficient
   of (w1 - w2)^2 is EXACTLY zero, bulk energy densities are single-region
   exact (locality), no volume term.
8. Envelope identities (general symbolic proofs): fixed-Q family E~(Q) =
   E(w) + Q^2/(2I) satisfies dE~/dQ = w exactly; interior inertia
   iota_int = -2 dV_min/dw^2 exactly (chain rule + stationarity).
9. Thin-wall bag: dE_w/dR = 0 at fixed w gives the Laplace/volume-surface
   law R = 2 sigma/p exactly within the reduced family; d^2E_w/dR^2 < 0
   there (spinodal at fixed w); fixed-Q stability carried by the Q^2/2I
   term; dE/dQ = w by 8. Limits: p -> 0 as w -> w_c => R -> infinity
   (wall-to-infinity, P249 limit).

## Maxwell (degenerate-depth) frequency

Deep branch: m = 0 exactly (dm Vw = m(...) = 0 identically). Critical
system after substituting w^2 = 32f^2 - 12f + 6 - 24b (from d_f Vw = 0):

  pA(c,b) = 8c^3 - 3c^2 + (8b^2+1)c - 3b^2
  pB(c,b,f) = db Vw = 0
  pC(c,b,f) = Vw = 0

High-precision root (mpmath, 50 dps, residual 5e-50):
  c* = 0.30280764518677380369, b* = 0.65773437772324925193,
  f* = 0.81436149699856776719,
  w_c^2 = 1.663945700059150298856193...

Exact elimination to the minimal polynomial of w_c^2: resultants
R1 = Res_c(pA,pB), R2 = Res_c(pA,pC) computed exactly; final
Res_f/Res_b chain running (elimination.py). Numeric landscape scan
(exploratory only, cannot enter claims): interior minimum negative already
at w^2 ~ 1.7, consistent with w_c^2 < 45/16 (the accepted witness ratio is
a provenance cross-check, not an input).

## Route verdicts this attempt

- WALL_STATIC_PROFILE: refuted for the objective's purpose (no interior
  state at w = 0: exact; radial static shell: exact Derrick). Route-scoped;
  planar-loop existence recorded as frontier, not a resting verdict.
- WALL_ROTATING_MAXWELL: established at the structural level (system, first
  integral, kink existence argument via tension functional coercivity);
  exact w_c^2 characterization pending the elimination confirmation.
- SLIP_LOCUS_DECOUPLING: established exactly (items 6-7).

## Licenses earned / missing

Earned: slice-reduction license, first-integral license, orbit-invariance
license, zero-mismatch license, envelope licenses, Laplace-law license.
Missing: exact algebraic isolation certificate for w_c^2 (elimination run);
exact nondegeneracy of the Maxwell point Hessian (pending the algebraic
point); sigma_c value enclosures (variational functional exact; value not
needed by any claim).

## Next routes

- Finish exact elimination + Sturm isolation; classify ALL real branches of
  the Maxwell system (deep vs any shallow branch).
- Exact nondegeneracy (Hessian minors in the algebraic field).
- Formalize claims C-M5W-001..005 statements; module + tests; mutation gates.
