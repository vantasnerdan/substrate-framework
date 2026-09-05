# Attempt 0020 — N3 ensemble Cauchy-Born moduli: (lambda, mu, alpha, c_tr, c_s, c_a) exact, verify_cst003 21/21

## Route

Ensemble Cauchy-Born coarse-graining with triads (obligation N3), built in
sympy parts with residual checks (parts/n3_part*.py), consolidated into
`verify_cst003.py` (CheckLedger, 21 checks, 3 mutations, exit 0).

## Declared model (per the N3 license chain)

- ensemble: isotropic tube ensemble, length density L_v, fixed Gamma, a, R
- segment energy: objective relative Green-Lagrange measure on
  eps_gamma = grad u - skew(Phi) (material frame indifference);
  tension T from Biot-Savart (straight_line_tension, P242 reuse)
- moments: exact sphere moments (P242 homogenization reuse), extended by
  joint tangent-projection contractions for the wryness sector

## Result — every modulus moment-exact, no fitted constant

Stretch sector (P242 affine matching applied to T):
    lam_eff = mu_eff = L_v T / 15,
    T = rho Gamma^2/(4 pi) ln(R/a),
    pre-stress P = L_v T / 3 I  (recorded for the N4 tangent operator)

Microrotation coupling (the Comparsi intake pair):
    W_alpha = (alpha/2)|rot u - 2 Phi|^2,  alpha_eff = L_v T / 6,
    dW/dPhi = -2 alpha rot u + 4 alpha Phi   (structure check: True)
  with T substituted: alpha_eff = L_v rho Gamma^2/(24 pi) ln(R/a)

Wryness sector (projected bend energy, joint tangent-triad moments;
C_tw^tube = 0 from 0005 gauge identity + 0011 virial closure):
    c_tr = -B L_v / 30,   c_s = B L_v / 10,   c_a = B L_v / 6
  with B the tube bend stiffness (log-running; declared N2 remainder
  B = rho Gamma^2/(4 pi)[ln(R/a) + c1], c1 = 1/2 - EulerGamma per 0015/0017).

## Checks and mutations

- Poincare-style moment reuse: second/fourth sphere moments exact;
  the axial-moment identity's general-eps correction
  (residual = -(antisym pairs^2)/15) recorded and verified.
- General asymmetric probe: ensemble energy == matched form, residual 0.
- The alpha sector carries the linear Phi coupling (dW/dPhi1 = L_v T (h12 - h21)/3).
- M1 wrong second moment (delta/4): rejected. M2 wrong fourth moment (1/21):
  rejected. M3 no MFD subtraction: rejected (probe difference -8 L_v T/75).

## Defects caught in-run

1. E_cross-free stretch assembly: the raw (non-objective) second-order stretch
   is not MFD; the relative Green-Lagrange measure is required and produces the
   alpha sector. (The vapid first M3 formulation — tr skew = 0 triviality —
   was replaced by a nonzero-Phi probe test.)
2. The two-invariant wryness matching does NOT close (three couple invariants
   needed); the three-invariant match closes exactly.

## Status

- route_verdict: established (N3 moduli receipt; per-modulus expressions above)
- evidence_scope: EXACT (sympy, moment-identities + matching + mutations)
- verify_cst003.py: 21 checks PASS, exit 0
- Next unlocked: N4 balance-law identification (Comparsi (1)/(2) form with the
  moduli above; microinertia j from the tube rotational self-energy — M_eff =
  pi rho a^2 virial, 0011 — feeds the spin-inertia normalization).
