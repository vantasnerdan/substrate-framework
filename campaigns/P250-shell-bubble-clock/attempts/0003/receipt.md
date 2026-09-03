# Attempt 0003 — G2 receipt: the bag (fixed-omega radial BVP), frozen before production output

Parent: P250 gap G2 (owner-directed continuation; same branch/PR #196).
Skill bindings: small-ratio-numerics (small quantities: sigma, p, the
surface-vs-volume balance; split pieces, crossed refinement, itemized
budget, extrapolate only with verified observed order); physics-erdos-loop
typed chain. Proposal S1 numerical_representation binds: exact algebra and
operator derivation FIRST; production numerics only after this receipt.

## Typed chain (frozen)

mathematical object: static spherically symmetric bag = stationary point of
  E_w = int d^3x [Tr(d_iS d_iS)/4 + |grad psi|^2/2 + V_omega(u)] on the
  no-winding slice ansatz S(r) = diag(m, c+b, c-b), psi = f(r) >= 0,
  u = (m,c,b,f)(r), at fixed omega with omega^2 in (omega_*^2, omega_edge^2).
symmetry/conservation license: spherical symmetry (no angular winding;
  fixed-basis slice ansatz, C-M5W-002 scope); rotating-frame stationarity
  at fixed omega (Legendre partner of the fixed-Q ensemble).
exact variational functional: E_w above with V_omega from
  m5_wall_clock.wall_slice_potential (Abs(f)->f licensed by C-M5W-001
  slice, verified against the certified Maxwell system in attempt 0002).
admissible function space/representation: u in C^2[0,L], regular at the
  origin (u'(0)=0), u(L)=A=(1,0,0,0); collocation (SciPy solve_bvp) with
  the Lommel singular-term matrix S = blkdiag(0_4, -2 I_4).
analytic scale and asymptotic structure:
  - radial operator (EXACT, verified symbolically this attempt):
      (2K)(u'' + (2/r) u') = grad V_omega, 2K = diag(1/2,1,1,1),
      K = diag(1/4,1/2,1/2,1/2)  [radial residual == planar residual +
      (2/r)(2K)u' exactly on rows m,c,b; f-row identical on the slice
      f>=0, difference off-slice = 4 f^2 f' (sign f - 1), zero on slice].
  - monotone identity (exact): d/dr[T - V_omega] = -4T/r on solutions.
  - virial identity (exact): int_0^L r (T + V_omega) dr = 0 for solutions
    (boundary terms vanish because V_omega(A) = 0 identically in omega).
  - interior bulk: deep-branch B(omega) = (0, c(w), b(w), f(w)) from the
    exact m=0 system dV/dc = dV/db = dV/df = 0 (Newton, 40 dps, warm
    continuation from the certified Maxwell point; det J monitored).
  - pressure (exact): p(w) = (w^2 - w_*^2) iota_B(w)/2 with iota_B =
    f^2 + 4b^2 at B(w); V_omega(A) = 0 for all omega identically.
  - exterior decay masses at A(omega): sqrt of (10, 10, 22-4w^2, 6-w^2);
    domain ceiling omega^2 < 5.5 (b-direction positivity, exact).
observable: R = bag radius read from the profile (m = 1/2 crossing;
  cross-checked by argmax T); sigma-balance quantities: E(w) =
  4 pi int r^2 (T + V_omega) dr, Q(w) = 4 pi w int r^2 (f^2 + 4b^2) dr,
  I(w) = 4 pi int r^2 (f^2+4b^2) dr.
irreducible numerical remainder: the collocation solution of the singular
  BVP and the quadratures of E, Q, virial; R-extraction ambiguity.
numerical approximation: solve_bvp, Lommel S-term, warm omega-continuation
  from thin-wall start R_guess = 2 sigma_0 / p(w) (sigma_0 = 0.72929841787
  from attempt 0002), tanh-shape guess per rung; compensated trapezoid
  (math.fsum) + adaptive Gauss-Kronrod cross-check.
permitted verdict: BAG_EXISTS (slice scope): constructed R(w) family with
  itemized budget, thin-wall law check chi = R p / (2 sigma_0) -> 1 as
  w -> w_*, envelope check dE/dQ = w along the family, virial/monotone
  identity checks. NOT licensed: global minimality, hedgehog winding,
  stability spectrum (G4), any claim promotion (candidates formalized
  after results exist).

## omega_edge (operational, frozen)

omega_edge^2 := the terminal point of the bag family = first of:
(a) analytic ceiling omega^2 = 5.5 (exterior b-direction mass^2 22-4w^2
    hits zero; BVP decay fails), or
(b) loss of the wall saddle / critical slowing of the BVP continuation
    (family fold), whichever comes first.
The deep bulk branch B(omega) itself does NOT fold below 6.2 (verified by
40-dps continuation, det J bounded away from 0), so the ceiling is set by
the wall/exterior side, not the bulk.

## Budget items (acceptance: each wall-localized item <= 1e-3 sigma_0 = 7.3e-4)

b1 collocation residual (scaled), b2 mesh refinement (n-ladder on rungs),
b3 domain tail (per-rung e^{-2 m_min L} bound + L-ladder spot check),
b4 certified-input drift (omega_*^2 box 1e-12 -> relative p drift),
b5 quadrature (GK vs compensated trapezoid, relative on E and Q),
b6 evaluator noise (fsum permutation), b7 R-extraction ambiguity
(|R_m - R_T|), b8 identity residuals (virial, monotone; normalized),
b9 branch identity (two R_guess starts -> same solution).
The thin-wall deviation chi(w) - 1 is a RESULT (convergence trend in
delta = w^2 - w_*^2), not a budget item. dE/dQ finite-difference
truncation is quoted with the check, not in the sigma-scale budget.

## Non-goals

No hedgehog/no-winding comparison, no off-slice minimality, no G4
spectrum, no claim text changes. Bag profile figure data + R(w) table are
deliverables; claim candidates C-M5W-006/007 are named only after results.
