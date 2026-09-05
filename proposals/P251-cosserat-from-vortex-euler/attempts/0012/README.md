# Attempt 0012 — N2: canonical receipt for m=2; m=1 finite-k pose diagnosis

## Receipt analysis (canonical wave energy, corrected system)

For solid rotation the Eulerian perturbation and the Lagrangian
displacement satisfy v' = D_t(xi) - (xi.grad)v0 with the basis-advection
exactly cancelling (xi.grad)v0 componentwise:

  v'_r = -I wt xi_r,  v'_th = -I wt xi_th,  v'_z = -I wt xi_z,

so xi = v'/(-I wt) (no degeneracy at wt = -Om; the earlier degeneracy
suspect was an artifact of mis-forming D_t xi).

Action and canonical energy (mode e^{i(m th + k z - w t)}):

  A = (1/2) Im int rho (xi* . v') dV,  E_canon = -wt A
  (sign convention fixed by matching the conserved mode energy).

m=2 channel, k -> 0 (fields of attempt 0011):

  xi = (r eta/a)(r_hat + I th_hat),  A = pi rho Om a^2 eta^2 / 2,
  E_canon = pi rho Om^2 a^2 eta^2 / 2 = T_max = Delta V_static.

- POSITIVE canonical energy (stable branch); the 0011 virial closure is
  CONFIRMED on receipt-grade bookkeeping. C_tw = 0 (single tube, O(k^2))
  REINSTATED for m=2.
- Why the plain functionals are exact for m=2 but not m=1: the
  second-order mean-flow sector v''_0 vanishes identically for the
  quadrupole polarization (source-free exterior: no axisymmetric
  harmonic; traceless interior: no O(eta^2) interior field), while for
  m=1 the displaced-vortex translation test (Delta E = 0 exactly)
  FAILS without the v'' sector (the O(xi^2) dipole-expansion field is
  required). The receipt-pending scope note appended to 0011 is
  therefore NARROWED: receipt-pending applies to the m=1 channel only;
  the m=2 claims stand as stated.

## m=1 finite-k route diagnosis (blocked with missing construction named)

The 3-residual truncated system (kinematic pair via amplitude
elimination, pressure continuity with advected Bernoulli, frozen-
vorticity tangential condition) closes EXACTLY for m=2 (F(root) ~ 4e-29,
attempt 0010) but does NOT close for m=1 at finite k:

  F_min(LIA scale, k=0.001) ~ 7e-3;  F_min(~0.35 k branch) ~ 1e-10 at
  k=0.001 rising to ~2.5e-7 at k=0.01 — finite minima, not zeros.

Missing construction (named): the m=1 helical sheet geometry requires
the axial sheet-strength sector of the matching conditions (the tilted
material sheet carries a v_theta AND v_z strength structure; the
truncated poloidal-only matching omits it, which first matters at
leading order precisely for the m=1 helix). The m=1 long-wave channel
(LIA branch omega ~ Gamma k^2 ln(2/(ka))/(4 pi) + c_inner and the
companion inertial branch) must be re-derived on the full helical
matching, with the v''-aware (or action-based) energy functional for
the modulus extraction.

## Status

- m=2 twist channel: COMPLETE (dispersion three-route, virial receipt-
  confirmed, single-tube C_tw = 0 exact at O(k^2)).
- m=1 bend channel: B(k) leading-log structure stands (attempt 0009
  exterior identity is receipt-independent); c_inner BLOCKED on the
  full helical matching construction named above. N2 remains active.
- Next route (0013): full helical sheet matching for m=1 (axial
  strength sector), then c_inner on the action-based functional;
  verify_cst002.py mutations m1-m5 as enumerated in 0011 plus
  (m6) truncated m=1 BC system must FAIL the closure test.
