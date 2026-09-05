# Attempt 0021 — N4 balance-law identification: EL equations == Comparsi micropolar form, verify_cst004 8/8

## Route

Euler-Lagrange identification (obligation N4): derive the balance laws from the
averaged action with the N3 energy, match component-by-component against the
Comparsi intake form (issue #198, 2026-09-03). All algebra sympy
(parts/n4_part*.py), consolidated in `verify_cst004.py` (8 checks, 3 mutations,
exit 0).

## Result — the derived micropolar system

Linear momentum:
    div sigma = (lam + mu - alpha) grad div u + (mu + alpha) lap u
                + 2 alpha rot Phi

Angular momentum:
    div m + eps:sigma = j Phi_tt
    eps:sigma = +2 alpha rot u - 4 alpha Phi        (spin coupling)
    dW/dPhi   = -2 alpha rot u + 4 alpha Phi        (intake pair, exact)
    div m = (c_s + c_a) lap Phi + (2 c_tr - c_a + c_s) grad div Phi

Microinertia (spin self-energy, orientation-averaged):
    j = L_v M_eff a^2 / 3,   M_eff = pi rho a^2 (0011 virial),
    I_axis = M_eff a^2/2, I_diam = M_eff a^2/4.

Standard micropolar correspondence: the couple modulus kappa of the textbook
form maps to the ensemble quantity alpha = L_v rho Gamma^2/(24 pi) ln(R/a)
(attempt 0020).

## Checks and mutations

- LINEAR balance (3 components): residual 0 against the standard form.
- SPIN pair (3 components): residual 0; action-reaction with the part-1
  dW/dPhi pair verified.
- div m wryness structure: (c_s+c_a) lap + (2c_tr-c_a+c_s) grad-div verified
  symbolically (first hand-derived coefficient 2c_tr-c_a was incomplete; the
  c_s grad-div piece is real and now part of the receipt).
- j: orientation-averaged spin energy == 1/2 j |Phi_dot|^2 exactly.
- M1 flipped coupling sign: rejected. M2 wrong diameter inertia: rejected
  (j != L_v M_eff a^2/3). M3 dropped alpha term: rejected.

## In-run defect

- .coeff() on Derivative patterns silently under-extracted (dropped the c_tr
  term in div m); replaced by direct full-expression symbolic targets.

## Status

- route_verdict: established (balance laws == Comparsi form, moduli from N3,
  microinertia from the tube spin self-energy; no fitted constant)
- evidence_scope: EXACT (sympy)
- Next unlocked: N5 symbolic two-branch dispersion (L4) — plane-wave
  substituion of the identified system; then N6 no-Cosserat side, N7 EPS bridge.
