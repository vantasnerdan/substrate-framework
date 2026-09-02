# Attempt 0001 — wall-sector reduction, Maxwell system, decoupling and bag laws

Active obligation: S1 (feeds S2/S3/S4 directly).
Routes: WALL_ROTATING_MAXWELL (primary, established), WALL_STATIC_PROFILE
(competing, refuted), SLIP_LOCUS_DECOUPLING (established), BAG_THIN_WALL
(established).

## Analytic specification receipt (exact layer, no production numerics consumed)

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
   representation frontier, not needed by the objective.
4. Exterior pencil re-derived from the canonical potential (closing the #190
   literal defect for P250 oracles): derived stiffness diag(5, 22, 10, 4, 4,
   6) in this derivation's channel labels (a, split, trace, u, v, psi_Im);
   the accepted C-M5C-002 chart labels its middle entries (t, p) in the
   opposite order, giving the accepted multiset. Only positivity of all
   channels is consumed downstream, and the reviewer verified all vacuum
   generalized masses strictly positive in either chart.
5. Excluded directions along clock-active branches (exact):
   shear Hessians 2(8b^2-3b+8c^2-3c+4m^2-3m+1) (+ cross 8v^2/8u^2),
   orientation Hessian 48 b f^2 (true rotation angle; positive on b,f>0,
   vanishing exactly at the locus = microscopic perfect slip). A half-angle
   parameterization gives 12 b f^2 for the same physical quantity.
6. Orbit-fixed locus: R(q) diag(m0,c0,c0) R(q)^T = diag(m0,c0,c0) exactly;
   inertia density 0 there; potential density, inertia density, velocity
   energy, and static gradient density are all EXACTLY orbit-invariant
   (symbolic identities, trigsimp-reduced to 0).
7. Phase slip calculus: a phase jump of any size at the locus costs
   <= L_f^2 (delta theta)^2 eps / 4 -> 0 (f(x*) = 0); tensor orientation
   jump is exactly free by orbit-fixation. Hence S2 mismatch coefficient
   of (w1 - w2)^2 is EXACTLY zero; bulk energy densities are single-region
   exact (locality); no volume term.
8. Envelope identities (general symbolic proofs): fixed-Q family E~(Q) =
   E(w) + Q^2/(2I) satisfies dE~/dQ = w exactly; interior inertia
   iota_int = -2 dV_min/dw^2 exactly (chain rule + stationarity).
9. Thin-wall bag: dE_w/dR = 0 at fixed w gives the Laplace/volume-surface
   law R = 2 sigma/p exactly within the reduced family; d^2E_w/dR^2 < 0
   there (spinodal at fixed w); fixed-Q stability carried by the Q^2/2I
   term; dE/dQ = w by 8. Limits: p -> 0 as w -> w_c gives R -> infinity
   (wall-to-infinity, the P249 exterior picture).

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

## Exact frequency-isolation route history (one verdict per route)

- Route A, lex Groebner elimination to a minimal polynomial: blocked
  computationally (kernel timeout after 20 min); superseded, not refuted.
- Route B, resultant chain MW = Res_f(Res_c(pA,pB), Res_c(pA,pC)) (degree
  108 in omega^2): computed exactly in ~31 s; factorization and degree-108
  exact Sturm counts are computationally prohibitive in SymPy rational
  arithmetic. Route-scoped verdict: exact global root enumeration NOT
  EARNED this route.
- Route C, PSLQ minimal-polynomial reconstruction: degree-5 relation at
  60 dps FALSIFIED by exact 90-dps re-evaluation (F(w2) = -1.4e-18, far
  above the 1e-80 evaluation noise of those coefficients); no relation up
  to degree 12 at 160 dps. Refuted as a shortcut; no claim rests on it.
- Route D (ADOPTED, certified): Krawczyk uniqueness of the deep-branch
  Maxwell solution in a rational box (halfwidth 1e-12, strict inclusion,
  interval arithmetic 40 dps) plus interval Gershgorin positive-definite
  fixed-omega Hessian (lower bounds 8.17, 4.88, 13.13), see certify.py and
  maxwell_certificate.json. Exact local isolation and nondegeneracy.
- Exact rational witnesses (independent of all routes): V_w at
  (0, 31/100, 13/20, 41/50) equals -13739/18750000 < 0 at w^2 = 5/3 and
  -33854777/25000000 < 0 at w^2 = 45/16, so w_c^2 < 5/3 < 45/16 exactly.

Numeric landscape scans in the kernel were exploratory only (labeled at
run time); they generated the deep-branch hypothesis and appear in no claim.

## Route verdicts this attempt

- WALL_STATIC_PROFILE: refuted for the objective's purpose (no interior
  state at w = 0: exact; radial static shell: exact Derrick). Route-scoped;
  planar-loop existence recorded as frontier, not a resting verdict.
- WALL_ROTATING_MAXWELL: established; exact w_c^2 characterized by Route D
  (unique Maxwell solution in the certified rational box, nondegenerate
  interior minimum) plus exact rational upper witnesses.
- SLIP_LOCUS_DECOUPLING: established exactly (items 6-7).
- BAG_THIN_WALL: established exactly within the reduced family (Laplace
  law, envelope closure, wall-to-infinity limit).

## Licenses earned / missing

Earned: slice-reduction license, first-integral license, orbit-invariance
license, zero-mismatch license, envelope licenses, Laplace-law license,
Maxwell local-isolation license (Krawczyk), nondegeneracy license
(Gershgorin). Missing: global SOS certificate identifying the deep branch's
crossing as the global inf of V0/iota (named closing route, honestly scoped
in C-M5W-004); sigma_c value enclosures (variational functional exact;
value not needed by any claim).

## Next routes

- Global SOS certificate for V0 - (w_c^2/2) iota >= 0 (closes the global
  identification in C-M5W-004's exclusion).
- sigma_c value: geodesic-path enclosures or an h-converged profile solve,
  only if a downstream consumer needs the number.
- Promotion transaction: claims C-M5W-001..005, review, release, docs,
  memory.
