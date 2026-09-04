# Attempt 0022 — N5 symbolic two-branch dispersion: verify_cst005 10/10

## Route

Plane-wave substitution in the identified N4 micropolar system, sympy-exact
(parts/n5_part*.py), consolidated in `verify_cst005.py` (10 checks, 3
mutations, exit 0).

## Result — the dispersion of the closed linear system

Transverse sector (coupled circular polarizations, 2x2):
    (rho w^2 - (mu+alpha) k^2)(j w^2 - (c_s+c_a) k^2 - 4 alpha) = 8 alpha^2 k^2
  Branch 1 (acoustic):   w^2_-(0) = 0
  Branch 2 (optical):    w^2_+(0) = 4 alpha / j      (gapped spin wave)
  Vieta sum/product consistent with the determinant.

Longitudinal sector (decouples):
    rho w^2 = (lam + 2 mu) k^2                (P-like wave; the alpha corrections
                                               of the two gradient terms CANCEL)
    j w^2   = 2(c_a - c_tr) k^2 + 4 alpha     (longitudinal spin wave)

Fluid limit (declared premise): at L_v -> 0 every modulus vanishes and all
branches collapse to w = 0 — Euler recovery.

With the N3/N2 values substituted (alpha = L_v rho Gamma^2/(24 pi) ln(R/a),
mu = L_v rho Gamma^2/(60 pi) ln(R/a), j = L_v pi rho a^4/3 via M_eff =
pi rho a^2), the optical gap is
    w_gap^2 = 4 alpha / j = (2/3) * Gamma^2/(pi^2 a^4) * ln(R/a) / pi
             = 2 Gamma^2 ln(R/a) / (3 pi^3 a^4)
— dimensionally 1/s^2, a pure tangle-scale spin-wave gap.

## Checks and mutations

- Transverse determinant structure with the 8 alpha^2 k^2 circular coupling.
- Branch limits k -> 0 (0 and 4 alpha/j); Vieta sum/product exact.
- Longitudinal sector projection identities (both diagonals).
- Fluid-limit collapse (Euler recovery premise).
- M1 dropped spin stiffness (gap vanishes): rejected. M2 wrong coupling
  strength (8 -> 4 alpha^2): Vieta breaks: rejected. M3 wrong c_tr sign:
  longitudinal spin k^2 coefficient changes: rejected.

## In-run defect

- First "P-wave without alpha correction" mutation was vacuous: the alpha
  correction cancels in the longitudinal sector by construction. Replaced by
  the c_tr-sign mutation on the longitudinal spin branch (a real error the
  verifier must catch).

## Status

- route_verdict: established (L4: symbolic dispersion fixes the production
  IVP observable — the branch speeds and the optical gap are the declared IVP
  inputs)
- evidence_scope: EXACT (sympy)
- Next unlocked: N6 no-Cosserat contrast side (orientation-ergodic ensemble
  yields vanishing couple moduli), then N7 EPS bridge.

## Erratum (2026-09-04, post-record re-derivation)

The inline closed-form simplification of the optical gap above contains an
arithmetic slip. From this record's own stated values
(alpha = L_v rho Gamma^2 ln(R/a)/(24 pi), j = L_v pi rho a^4 / 3):

    w_gap^2 = 4 alpha / j = Gamma^2 ln(R/a) / (2 pi^2 a^4),

NOT 2 Gamma^2 ln(R/a)/(3 pi^3 a^4) as written above. The verifier is
unaffected: verify_cst005 checks w^2_+(0) = 4 alpha / j symbolically in
(alpha, j) and never uses the slipped closed form.
