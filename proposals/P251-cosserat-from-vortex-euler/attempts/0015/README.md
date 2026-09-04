# Attempt 0015 — N2 frozen: verify_cst002.py 29/29 with 5-mutation battery

## Route

Frozen verifier `verify_cst002.py` for C-CST-002 (bend/twist micro-moduli),
executed per the announced route 0015. First execution: CheckFailure. The
verifier caught two further defects in the 0014-state and one in the
0006-0011 state before going green — exactly the failure modes the
mutation/audit discipline exists to catch.

## Defects caught and repaired in-run (all recorded, append-only)

1. **Interior Bessel species (caught by the frozen symbolic check).** The
   0014 record wrote the interior as P = A I_m(mu r). The symbolic
   incompressibility structure check failed on first run: the
   Poincare-consistent interior is P = A J_m(lam r) with
   lam^2 = k^2 (4 Om^2 - wt^2)/wt^2 and cylinder ODE
   P'' + P'/r - m^2 P/r^2 = -lam^2 P. The I-species field does not satisfy
   incompressibility for any sign of its parameter (the first-run
   coefficient k^2 + wt^2 mu^2/D proved it). The 0014 statement of lam^2
   was correct; the species was not.
2. **m=1 leading-log constant (caught by the sympy-exact truncation).**
   Exact rational bookkeeping (w = Om*sigma*t*L substituted, series through
   O(t)) gives r2 = -eta[2w + Om(ka)^2 (L + 1/2 - gamma)]:
   **c1 = 1/2 - gamma**, superseding the 0014 manual value -1/8. The
   proportionality r1 = rho a Om r2 is EXACT at the retained order
   (symbolic residue 0).
3. **m=2 dispersion coefficient (caught by the numeric closure check).**
   The m=2 (ka)^2-coefficient is species-dependent: the (lam a)^4 terms of
   the interior field enter the (ka)^2-order through the amplitude's
   (lam a)^-2 denominators. Measured root: c = -0.1665 +- 0.0005 across
   k = 1e-2..1e-4:
   **omega = Om(m-1) - Om(ka)^2/6 + O((ka)^4)**,
   superseding the recorded -1/12 of attempts 0006-0011 (that machinery
   ran on the I-species). The Kirchhoff k^0 limit (w/m -> Om/2) and the
   static virial are unaffected; the O(k^2) single-tube C_tw = 0 receipt
   is flagged for recomputation at N3.

## Frozen result (ALL 29 CHECKS PASS, exit 0)

- Symbolic: Poincare reduction (general m) + closure; plane-wave inertial
  dispersion wt^2 = 4 Om^2 kz^2/|k|^2; exact proportionality; root slope 2
  and constant 1/2 - gamma.
- Numeric (mpmath dps=40): m=1 LIA root tracked at k = 1e-2..1e-5 within
  25% of the asymptote, r1 = r2 at every root, residuals -> 0 across
  decades; m=2 coefficient -1/6 at three decades, |r| <= 1e-6 at roots;
  Kirchhoff pattern-speed limit; translation neutrality at k -> 0, w = 0.
- Mutations (all rejected): M1 defective 1/omega wavenumber (the
  0008-0013 defect, now a permanent named mutation), M2 single-Omega
  Coriolis, M3 flipped frozen strength, M4 wrong exterior Doppler, M5
  growing exterior harmonic.

## Status

- Route verdict: **established as stated** for the C-CST-002 branch-
  structure scope: bend channel w* = -(Gamma k^2/4pi)(L + 1/2 - gamma)
  [NUMERIC_EVIDENCE, ASYMPTOTIC], twist channel w = Om(m-1) - Om(ka)^2/6,
  translation neutrality, Kirchhoff limit.
- o(1)/O((ka)^4) remainders named (matched-asymptotics refinement, N3
  material).
- Owed: classical corroboration (c1 and the m=2 coefficient against a
  Kelvin/Saffman-level source, cited-not-recalled); C_tw receipt
  recomputation under the J-species; N2 -> N3 handoff (ensemble moduli).
