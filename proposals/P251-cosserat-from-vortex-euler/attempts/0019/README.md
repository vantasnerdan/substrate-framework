# Attempt 0019 — N3 canonical energy through O(k^4): e2 = 0 exactly, e4 log-exact, C_tw normalization delivered as coefficient map

## Route

O(k^4) canonical-energy receipt for the m=2 polarization branch under the
angle-wave (twist-strain) normalization, completing the 0018-reframed route.
All algebra in sympy, built and verified in parts (parts/*.py), arbitrated at
each step against exact-field mpmath numerics.

## Method (parts, each verified before use)

- `parts/mode_solution.py` — fresh normal-mode derivation about solid rotation.
  Proven compact forms (all momentum residuals zero, sympy):
      v_r = -i(wt p' + 2 Om m p/r)/(rho D),  v_t = (2 Om p' + m wt p/r)/(rho D),
      v_z = k p/(rho wt),  D = wt^2 - 4 Om^2
  and the pressure equation reduces EXACTLY to Poincare
  p'' + p'/r - m^2 p/r^2 + lam^2 p = 0, lam^2 = k^2(4 Om^2 - wt^2)/wt^2
  (residual 0). This supersedes the recorded 0006 D = (w-mOm)^2 - Om^2 form
  (an artifact of an i-factor sign slip in v_theta; sympy-arbitrated).
- `parts/fields_k4.py` — mode fields through O(k^4) on the pose-independent
  branch w = Om(1 - (ka)^2/6 + c4d (ka)^4). Checks: kin-in True, xi_r(a) = eta True.
  Defect found & fixed mid-part: branch must carry (ka) factors, not k alone
  (non-conformal a-powers in the energy coefficients made this visible).
- `parts/energy_k4.py` — interior kinetic + Arnold/Coriolis cross (functional of
  0018) + exterior potential; E_cross form ARBITRATED at k=0 against the
  +pi/2 anchor: winner Re[v*.(w0 x xi)] (README-0018's "xi x omega_0" prose had
  the cross-product order flipped).
- `parts/energy_num.py` — exact-field numeric arbitration: interior symbolic
  coefficients confirmed to 1e-6; 0018's recorded k^2 coefficients (-5/36,
  -13/36, c2 = -0.061) NOT reproduced by the exact fields — artifacts of 0018's
  truncated construction, superseded.
- `parts/energy_ext2.py` — exterior via manual Laurent inversion. sp.series on a
  reciprocal of a besselk expression silently drops log terms (reproduced twice);
  the log is REAL: 4-param fit rms 1e-12 vs 1e-9, f_log pinned.
- `parts/energy_can.py` — final assembly + validation.

## Result

    E_can/rho = pi Om^2 eta^2 a^2 [ 1  +  0         (no k^2 term)
                                    + e4 (k^4) ]
    e0 = a^2                                   (E_can(0) = pi rho Om^2 eta^2 a^2, 0018 anchor kept)
    e2 = 0 EXACTLY, all a                      (supersedes 0018's dE/dk2 = -0.53 artifact)
    e4 = -a^6 (3456 c4d + 72 log(k a) - 53 - 72 log 2 + 72 EulerGamma)/2304
       = -a^6 (3456 c4d + 72 log(k a / 2) + 72 EulerGamma - 53)/2304

Numeric validation (a = 1, c4d = -0.253, exact fields vs symbolic):

    k = 0.05: diff 2.5e-10;  k = 0.10: 1.5e-8;  k = 0.15: 1.6e-7;  k = 0.20: 8.4e-7
    (growth consistent with dropped O(k^6) terms)

The k^4 log(k a) coefficient is -a^6/32 (matches 4-param fits -0.0305..-0.0310);
the dispersion-c4d coefficient is +3456/2304 = 3/2 per |c4d|.

## Interpretation for C_tw (twist-strain normalization)

- The k^2 coefficient of the canonical energy VANISHES: the 0018 "rotation
  coupling at k^2, dE/dk^2 < 0" framing is superseded — there is no k^2 term.
  The branch's Krein signature lives in the dispersion, not in a k^2 energy shift.
- The mode's axis-angle strain is dphi/dz = -k/2 (pattern-angle ramp, exact for
  the helical m=2 mode). The couple-modulus energy is the O((ka)^4) coefficient
  with the polarization amplitude eta as the strain carrier.
- Absolute C_tw convention remains deferred to N4 (Comparsi (1)/(2) coupling
  identification fixes the phi-normalization), per frontier 0018. This attempt
  delivers the complete typed E_can(k) specification the convention maps onto.

## Anchors & regressions (all passing)

- E_kin_in(0)/rho = pi/4 Om^2 eta^2 a^2; E_cross_in(0)/rho = pi/2; E_kin_ext(0)/rho = pi/4;
  E_can(0)/rho = pi Om^2 eta^2 a^2
- Poincare reduction residual 0; momentum residuals 0; kin-in/xi sheet checks True
- e2 conformality: e2 = 0 for ALL a (an a^2(1-a^2) non-conformal artifact appeared
  mid-part with the k-vs-ka branch bug and vanished after the fix)

## Defects caught in-run

1. v_theta i-factor sign (Part 1 candidate forms) — arbitrated against solve output.
2. Prefactor pi vs pi/2 on interior integrals (k^0 anchors caught it).
3. E_cross form sign/order — k=0 anchor arbitration.
4. Winner-form built from k=0-substituted fields (general-k rebuild).
5. sp.series dropping log terms on reciprocals (twice) — manual Laurent inversion.
6. Dispersion branch missing (ka) factors — non-conformal a-powers caught it.

## Status

- route_verdict: established (E_can through O(k^4) delivered, sympy-exact,
  numerically validated; e2 = 0 supersession recorded)
- evidence_scope: EXACT (symbolic) with numeric cross-validation
- Next: N3 ensemble half (Cauchy-Born moduli via triad moments), then N4
  balance-law identification (fixes the absolute C_tw convention).

## Addendum (post-0020): forms arbitration against verify_cst002

`parts/forms_arbitration.py` — verify_cst002's check_poincare_reduction records
v_r = +i(wt p' - 2 Om m p/r)/(rho(4 Om^2 - wt^2)), v_t = (2 Om p' - m wt p/r)/(rho D).
Run against the linearized momentum residuals these FAIL (nonzero r- and theta-
residuals, sympy-exact); the part-1 forms pass all three residuals. Mechanism
named for the 0018 supersession: the 0018-lineage energy assemblies used the
wrong-signed velocity forms, which (with truncation) produced the spurious
k^2 coefficient; 0019's e2 = 0 stands on the residual-verified forms. The
verifier's own branch-structure conclusions (Poincare form, lam^2, c1 = 1/2 -
gamma, c2 = -1/6) are unaffected — they were corroborated independently by the
Kelvin-96 roots (attempt 0016/0019a).
