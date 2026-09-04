# Attempt 0027 — two-region Rankine eigenproblem: naive sheet system REFUTED numerically; frame-argument target established

## Review context

PR #199 review directs the decisive construction: an Euler-derived interaction
whose EXACT second variation locks macro- and micro-rotation, replaying N2 on
the momentum-residual-correct forms of attempt 0019. This attempt opens the
correct two-region formulation (core rotation Omega_i, ambient solid rotation
Omega_o) — the minimal problem containing the frame-locking physics.

## Formulation (`parts/two_region_system.py`)

Linearized Euler (0019-correct forms) in EACH region with its own Doppler:
  interior: wt_i = w - m Om_i, D_i = wt_i^2 - 4 Om_i^2, species J_m(lam_i r),
            lam_i^2 = k^2 (4 Om_i^2 - wt_i^2)/wt_i^2
  exterior: wt_o = w - m Om_o, D_o = wt_o^2 - 4 Om_o^2, species K_m(lam_o r)
Base state carries a vortex sheet at r = a with jump (Om_o - Om_i) a.

## Negative result: the naive 3-condition sheet system is REFUTED

System {kin_i, kin_o, displaced-pressure}:
  u_r^i(a) = u_r^o(a) = -i w eta
  p'_i(a) - p'_o(a) = eta rho a (Om_o^2 - Om_i^2)   [displaced base jump]
Exact-Bessel numeric root-finding (`parts/numeric_arb.py`, mpmath 30 digits):
  - single-tube limit (Om_o = 0) does NOT reproduce the recorded neutral
    branch omega = Om_i[(m-1) - (ka)^2/6 + ...]: roots come out COMPLEX with
    Im(w) growing ~ linearly in ka (e.g. ka = 0.1: w = 0.9962 + 0.0705 i),
    i.e. a spurious growing branch;
  - the sympy 4-condition variants (adding the quasi-static frozen-strength
    row [v_t'] = 2 (Om_i - Om_o) eta) are linearly DEGENERATE (all 3x3 minors
    vanish identically), confirming that condition as stated is dependent —
    as expected, because for material vorticity the strength dynamics is the
    curl of the momentum system, not an independent scalar law.

## Mechanism named

The naive system lacks the sheet-strength DYNAMICS: for a sheet separating
two solid-rotation regions, the perturbation strength evolves by advection
on the sheet with the mean angular velocity (Om_i + Om_o)/2 and couples to
the displacement through the base jump; neither the quasi-static jump nor
its omission is correct. Equivalently: the Rankine-sheet eigenproblem needs
the strength-evolution equation as an independent relation, which reduces
the classical problem to a contour-dynamics (Kirchhoff-Love) formulation.

## Frame-argument target for the corrected system

In the frame rotating with the ambient, the exterior is quiescent and the
patch carries vorticity contrast 2 (Om_i - Om_o); the k->0 m=2 limit must be
    w(k->0) = 2 Om_o + (m-1)(Om_i - Om_o)|_{m=2} = Om_i + Om_o,
reducing to the recorded Om_i (m-1) at Om_o = 0 and to 2 Om at Om_o = Om_i
(solid body: the patch boundary is fictitious). Any corrected matching
system or contour-dynamics discretization must reproduce this limit.

## Route forward (0028+)

1. Contour-dynamics with ambient rotation (the validated 0009/0010 route
   family): boundary-integral equation for the displaced contour in the
   rotating frame, mode dispersion from the contour kernel — bypasses sheet
   matching entirely and closes on the frame-argument limit above.
2. Exact canonical (wave-action) energy of the m=2 mode as a function of
   (Om_i, Om_o) on the corrected branch: second variation in the relative
   rotation = the Euler-derived locking stiffness (candidate
   E_lock/L = (pi rho/2)(Omega_amb - Omega_self)^2 a^2 eta^2 per 0026 part 2
   is the static-frozen special case to be recovered).
3. Then the 0028 ensemble contraction and N3 rebuild.

## Status

- route: two-region eigenproblem OPENED; naive sheet conditions REFUTED
  (numeric, mechanism named); corrected target established exactly.
- evidence_scope: EXACT numerics (mpmath 30-digit) for the refutation;
  frame-argument limit exact.
