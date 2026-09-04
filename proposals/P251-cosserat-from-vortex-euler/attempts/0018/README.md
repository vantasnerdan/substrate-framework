# Attempt 0018 — N3 canonical-energy receipt: Krein paradox resolved, k^2 coefficient reframed, C_tw located at O(k^4)

## Route

Complete the 0008 licensed route (a): canonical wave-energy receipt for
the m=2 polarization branch, replacing the fixed-eta kinetic energy
(the wrong functional, negative-definite at O(k^2)). All algebra in
sympy/mpmath per the standing directive; two independent energy forms
computed as cross-validation.

## Structural result: the canonical energy and the missing piece

Canonical modal energy (pseudo-energy of the steady base):

    E_can = E_kin + E_cross,
    E_kin  = (rho/4) Int |v'|^2 dV,
    E_cross= (rho/4) Int Re[v'* . (xi x omega_0)] dV,
    omega_0 = 2 Om Theta(a-r) zhat.

Interior fields: corrected sympy momentum solution (0017 conventions),
P = A_p J_2(lam r), lam^2 = 3 k^2 - 4k^4/3, A_p from the exact kin
condition at the sheet through O(k^2); displacement xi = v'/(-I wt)
(circular polarization: xi_th = I xi_r; xi_r(a) = eta at leading
order). Exterior: K_2 potential, kin-out normalized.

Results (rho factored, a = 1):

    E_kin_in /rho   = pi Om^2 eta^2 (209 k^4 - 1440 k^2 + 2592)/10368
    E_cross_in /rho = pi Om^2 eta^2 (521 k^4 - 1872 k^2 + 2592)/5184
    E_kin_ext /rho  = (pi Om^2 eta^2/2) (0.49983 + c2 k^2 + ...),
                      c2 = -0.061 +- 0.001 (mpmath, k = 0.03..0.3)

Anchors, all passing:
- E_kin_total(0) = pi/4 + pi/4 = pi/2 * rho Om^2 eta^2 a^2 — EXACTLY the
  0007 recorded mode energy. The 0007 functional is the KINETIC energy.
- E_can(0) = pi rho Om^2 eta^2 a^2 — the canonical energy is TWICE the
  kinetic. The missing +pi/2 piece is E_cross (the Coriolis/Arnold
  work term v'.(xi x omega_0)): precisely the piece whose absence made
  the 0008 fixed-eta sum negative. Krein paradox RESOLVED with
  mechanism: wrong functional, missing cross-term, both now quantified.
- E_cross(0) = +pi/2 rho Om^2 eta^2 a^2 > 0 for wt < 0: the step's
  Coriolis work stabilizes the sign.
- The v'-form total agrees with the Arnold/Dt-form total at k = 0
  (interior Dt xi = 0 there, sheet Arnold term = +pi rho Om^2 eta^2):
  two independent forms, one answer.

## Reframing: the k^2 coefficient is NOT C_tw

    dE_can/d(k^2) = pi rho Om^2 eta^2 a^4 * (-1/2 + c2/2) < 0.

Negative, and correctly so: the m=2 polarization wave is an ANGLE wave
(axis-angle field phi(z) = -kz/2 for the e^{ikz} mode). Its Cosserat
twist strain is d phi/dz with amplitude O(k): the couple-modulus energy
(1/2) C_tw |d phi/dz|^2 sits at O(k^4), not O(k^2). The O(k^2)
coefficient is the rotation-coupling (Krein signature of the branch,
consistent with d w/d k^2 < 0 at positive canonical energy) — NOT the
couple modulus. The 0008 framing ("C_tw receipt at k^2, flagged
negative") was a category error; superseded here.

## Consequence for the receipt

C_tw must be extracted from the O(k^4) coefficient of E_can under the
angle-wave normalization (strain amplitude |d phi/dz| = |k chi0|, with
the polarization amplitude eta expressed through chi0). This requires
the fields through O(k^4), the exterior c4, and the wt-branch
corrections — named route 0019.

Also recorded: the wave-ACTION route (integral of |xi|^2) is
critical-layer singular for this branch (coradial point
r_c = a sqrt(2 Om/w) with w ~ Om sits outside the core); the ENERGY
route is finite (E_cross has no exterior piece since omega_0 = 0
there). Route (a) proceeds on the energy form only.

## Status

- Route verdict: established as stated for the E_can structural receipt
  (Krein resolution + k^2 reframing); C_tw numeric receipt routed to
  0019 at O(k^4).
- Mutations caught in-run: spurious K2^2 factor in the exterior
  integrand (produced a fake 1/k^4 divergence, caught by the analytic
  k->0 dipole limit 0.5); series-order trap (O(t^0) cut returns 0=0).
- Owed: 0019 O(k^4) receipt and the twist-normalization; ensemble
  moduli (Cauchy-Born via triad moments); N4 balance-law identification
  against Comparsi (1)/(2).
