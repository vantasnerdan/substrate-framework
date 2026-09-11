# T2 back-reaction mini-freeze — FROZEN PRE-COMPUTE (beacon, 0157/tilt)

Status: FROZEN 2026-09-12. Charter: shepherd T2 BACK-REACTION
(source demonstrated edge-T1'; back-reaction next). Freeze commit
precedes compute commit. Fences inherit (static tilt, FB-C1 never
fires, 0160 out, no measurement/carrier claim).

## Object

Tilt background -> edge defect: (i) mobility shift via
Peach-Koehler force from tilt-induced prestress; (ii) core-energy
shift vs B1 bars (P2 scoped to edge); (iii) Peierls-type barrier
from a tilt wave. Reciprocal of Vikulin via Maxwell structure on
the banked coupling — no new formalism.

## Frozen derivation structure (pre-compute)

Coupled density (banked J(eps) + tilt): w = (1/2)J(eps).th^2 +
W_FA(eps), J(eps) = J0(1+4.eta_in-3.eps_zz), J0/ell^3 = energy
density scale (per-cell volume ~ ell^3, FIDUCIAL geometry below).
Tilt-induced prestress: sig_tilt = -d w/d eps = -(J0/ell^3).th^2
x diag(2,2,0)... exact tensor form derived in compute (SymPy):
in-plane dilatational prestress propto th^2.
Peach-Koehler on edge (b=xhat, t=zhat): F = (sig.b)x t.
FROZEN PREDICTION: prestress diagonal -> F_glide = F_x IDENTICALLY
0 at O(th^2); F_climb = F_y = +2(J0/ell^3).th^2.b (sign receipted
in compute). Mobility verdict: glide UNAFFECTED at O(th^2); climb
activated but needs mass transport (NOT priced — displacement
quasi-static only, stated limit).
Core shift per length: dE_core ~ sig_tilt.b^2 (geometric O(1)
receipted) -> relative shift = beta.th^2 x O(1), beta =
(J0/ell^3)/mu.
Peierls: imposed tilt wave Th.cos(qx) -> barrier/length
dV ~ (J0/ell^3).Th^2.b/q x O(1) (estimate-level, labeled).

## Fiducial geometry (DECLARED assumptions, not banked)

Cell ~ ell^3/ring; line density Lam ~ 2pi.R/ell^3;
mu_aff = G^2.Lam.ln(ell/a)/(40pi) (banked form);
R/ell = 0.1 (multipole footnote ell>>R); ell/a = 100; th = 0.1
(small-tilt remainder O(th^3) respected).
Parametric beta = 10pi/ln(ell/a).(R/ell)^3 ~ 6.8e-3; effect
beta.th^2 ~ 7e-5. THE verdict is the parametric inequality; the
fiducial number illustrates. Reopen iff geometry differs >10x
(stated sensitivity).

## Frozen falsifiers + stop (scoped to edge)

- F1 (guard): SymPy F_glide != 0 at O(th^2) -> implementation bug,
  never physics. F2 (physics): fiducial relative core shift >
  B1 bars (7e-4) -> P2 FIRES FOR EDGE (core unbound by tilt);
  below -> BOUND, mobility verdict stands.
- F3 (physics): Peierls barrier >= core energy at fiducial
  (pinned by infinitesimal tilt -> mobility kill) or <= 0 with
  wrong sign (bug) -> disclose either; barrier << core expected.
- STOP: bank T2 verdict (bound/unbound + mobility + barrier) ->
  T3 readout-staging proposal. No formalism beyond PK + banked
  J(eps) + Volterra + fiducial cell. T2-edge needs no further
  charter once this freeze is acknowledged.
