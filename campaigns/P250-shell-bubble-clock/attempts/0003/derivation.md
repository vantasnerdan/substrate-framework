# Attempt 0003 — G2: the bag (fixed-omega radial BVP), derivation record

Status: established (numeric evidence, BAG_EXISTS slice scope)
Parent: P250 gap G2 (owner-directed continuation; same branch/PR #196)
Route: receipt 0003 (frozen before production) -> exact radial operator ->
domain-reduced clamped-window collocation with amplitude-matched clamps ->
delta-continuation from the thin-wall limit -> mesh-stability acceptance.

## Exact analytic layer (verified symbolically, pre-solver)

- Radial EL on the no-winding slice S(r) = diag(m, c+b, c-b), psi = f(r):
  (2K)(u'' + (2/r)u') = grad V_omega with 2K = diag(1/2, 1, 1, 1).
  Verified: radial residual == planar residual + (2/r)(2K)u' as a zero
  polynomial on the m, c, b rows; the f row differs only by
  4 f^2 f' (sign f - 1), which vanishes exactly on the licensed slice
  f >= 0 (the Abs(f) suspension — same phenomenon as attempt 0002).
- Monotone identity (exact): d/dr (T - V_omega) = -(4/r) T on solutions.
- Virial identity (exact): int_0^L r (T + V_omega) dr = 0 for solutions
  (V_omega(A) = 0 identically in omega kills both boundary terms).
- iota = f^2 + 4b^2 exact; p(w) = (w^2 - w_*^2) iota_B(w)/2 exact
  (V_omega(A) = 0 for all omega; V_0(B) = (w_*^2/2) iota_B by Maxwell).
- Deep-branch bulk B(w) = (0, c, b, f)(w): 40-dps Newton continuation of
  the exact m = 0 system dV/dc = dV/db = dV/df = 0 from the certified
  Maxwell point; branch residual < 1e-30 per rung; det J bounded away
  from 0 past w^2 = 6.2 (no bulk fold below the exterior ceiling).

## Route history (all recorded, none hidden)

1. SINGULAR-TERM FORMULATION (S = blkdiag(0, -2I)) on the full domain
   [0, R+L]: Newton stalls — the residual metric scales as col_res/h, so
   refinement past the Newton stagnation point RAISES the metric; also
   the trivial A-everywhere solution competes in its basin. Blocked.
2. ROOT CAUSE of one stall generation: make_rhs omitted the -(2/r)u'
   term after dropping the S matrix — solved the planar ODE by mistake.
   Fixed; residual floor dropped ~100x. (Lesson recorded.)
3. KINK-SHAPE GUESS orientation bug (bag(r) = kink(r-R) instead of
   kink(R-r)) and cubic-spline extrapolation outside the kink's data
   range: produced garbage guesses (clamped to [-12, 12] and flipped).
4. TRIVIAL-BASIN collapse: tanh/kink guesses can fall into u = A
   (an exact solution of the BVP); guarded by interior-agreement checks.
5. ADOPTED: DOMAIN REDUCTION — the profile differs from the bulk states
   by e^{-m_B (R-r)} inward and e^{-m_A (r-R)} outward, so the BVP is
   solved on [R-half, R+L_ext] with clamped BCs u(R-half) = B(w),
   u(R+L_ext) = A; the flats enter ANALYTICALLY (interior T = 0,
   V(B) = -p exactly; exterior V(A) = 0 exactly). Clamp distances are
   set by FIXED MATCHING AMPLITUDE m*dist = 9.2 (amplitude e^{-9.2} ~
   1e-4; collocation conditioning e^{18.4} ~ 1e8, inside float64 —
   clamping deep in the flat zone pushes conditioning past 1e16 and
   stalls Newton, observed and recorded). No singular term needed
   (r >> 1 on the window). Exact analytic Jacobians supplied.
6. DELTA-CONTINUATION: rung 1 at delta = 1e-3 seeded by the attempt-0002
   kink shape at R = 2 sigma_0/p (O(delta)-accurate); later rungs warm
   seeded from the previous bag rigid-shifted to the new R, with fresh
   kink-guess and recenter retries. Thick rungs (R < 8) switch to the
   full-domain S-formulation (not needed below the delta = 0.008
   terminal).
7. ACCEPTANCE: mesh-stability — each rung is solved at n_wall = 1441
   and 2881 and accepted when |R(1441) - R(2881)| <= 1.5e-3 R with
   res_rms <= 2e-3 (the residual metric scales as col_res/h and
   saturates under refinement, so it is QUOTED as b1, not used as the
   acceptance criterion; receipt b1 budget 7.3e-4 covers the wall-local
   quantities through the measured dE = 1.2e-5 relative).

## Result

- R(w) family over delta = w^2 - w_*^2 in [1e-3, 7e-3] (7 rungs),
  R: 1217.38 -> 171.15; chi = R p / (2 sigma_0):
  0.99956204, 0.99926234, 0.99854712, 0.99748919, 0.99538755,
  0.99225955, 0.98760122 — the thin-wall law R = 2 sigma_0/p holds
  asymptotically (chi -> 1 as delta -> 0) and deviates linearly in
  delta beyond it (chi - 1 ~ -1.75 delta), consistent with finite-wall
  corrections.
- Envelope identity dE/dQ = omega: verified numerically along the
  family at 1.8e-4 - 2.2e-4 relative (finite-difference limited;
  the identity itself is exact — accepted S3/8).
- E at rung 1 = 4.53902671e6 vs the exact thin-wall critical energy
  16 pi sigma_0^3/(3 p^2) = 4.53197696e6: 0.16% (the expected O(delta)
  beyond-thin-wall correction).
- Virial identity: |int r (T+V)| / (int r (T+|V|)) <= 1.2e-3 across
  rungs; monotone identity at the same level (both dominated by the
  spline-differentiation error of the check itself).
- Mesh refinement (b2): dR = 4.1e-1 (1.3e-3 relative) between n = 1441
  and 2881; window-width sensitivity (b3): dR = 7.2e-1; restart (b9):
  dR = 6.8 (window-following artifact of the restart, recorded).
- Quadrature (b5): grid-halving energy difference ~1e-5 relative.
- Terminal point: delta = 0.008 (the residual floor ||J|| tol / (1.5 h)
  grows with the Hessian norm along the family; NOT a physics
  obstruction — the family continues by construction). Thick-wall
  continuation (R ~ few) needs the regular-origin formulation with a
  different guess strategy: named frontier, not pursued here.

## Does not license

Global minimality, hedgehog winding, stability spectrum (G4), any
statement beyond the slice ansatz, claim promotion (candidates
C-M5W-006/007 formalize the G1/G2 numbers for individual review).

## Reproduce

PYTHONPATH=src python proposals/P250-shell-bubble-clock/attempts/0003/bag_bvp.py
(threads pinned to 1; deterministic; ~8 min)

## Artifacts

bag_results.json (rungs + refinement + envelope), profile_bag_thin.csv,
profile_bag_thick.csv, run_stdout.txt (first-execution capture).
