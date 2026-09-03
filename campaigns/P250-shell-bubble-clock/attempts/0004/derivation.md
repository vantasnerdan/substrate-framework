# P250 G3 attempt 0004 — exact KKT closure of ω_c² = ω*² on the slice

## Status: PARTIAL — case A and case B closed exactly; case C root counting
## pending as an independently runnable certificate run

User directive (2026-09-02): the end-to-end certificate run exceeds 30 min in
the exact number-field Sturm machinery; the PR proceeds without it. Optimizing
the script is not needed for this PR or its review. The script
`certificate.py` is committed so the remaining run (case C counting +
Newton/Krawczyk + value comparison + `OMEGAC_CLOSURE` verdict) can be
completed and rerun independently later.

## What is certified (exact, printed in run_stdout.txt)

1. `critical system == grad V verified exactly` — the KKT system in
   (m, s, d, f) = (m, 2c, 2b, f) equals ∇V of the slice potential
   V(m,s,d,f; t) symbolically (SymPy `simplify` of the difference).
2. Coercivity: leading quartic V4 = q² + 8f⁴ > 0 (q = m² + (s²+d²)/2);
   sphere-direction min 8/9 > 0 ⟹ V coercive for every t ⟹ infimum is
   attained at a critical point.
3. Case A (f = 0, d = 0): exact enumeration {(1,0): V=0},
   {(−¼,0): V=125/256}, {(0,0): V=1/2}; both non-origin values > 0 at
   t = α (α ≈ 1.66394570005915) and for all t ≥ 0 for the (0,0) value.
4. Minimal polynomial μ of ω*²: degree 32, irreducible over Q (leading
   coeff 20862942865979867136…), 4 real roots; the attempt-0001 Krawczyk
   enclosure [A0, A1] (width ≈ 1.25e-37, 44-decimal endpoints) isolates
   exactly one μ-root (rational Sturm count == 1).
5. Case B (f = 0, d ≠ 0): empty at t = α — both discriminants
   discB1 = (8t−27)(8t−3), discB2′ = 9(64t²−112t−367) certified strictly
   negative at α by exact field sign evaluation (both −1).
6. Module defect fixed: `maxwell_frequency_resultants` symbol collision
   (fresh `f` vs `Symbol('f', nonnegative=True)`) — R1/R2 were built from
   two distinct symbols printed identically. Fixed in canonical source;
   `s≠0`, `d≠0` forced in case C (else 32f²−12f+6 = t, disc 128t−624 < 0).

## What remains (future independent run)

- Case C (f ≠ 0): count real roots of squarefree R_C1 (deg 14) and R_C2
  (deg 17) at t = α in the exact field Q(α) (machinery in place:
  `_chain`, `_variations`, `count_roots_field`); solve each resultant root
  back for (s, f, m, d) (Newton 80 dps + Krawczyk 160 dps, seeds present);
  evaluate V on each critical box by interval arithmetic; compare with
  V_ω*(A) = 0. Print `OMEGAC_CLOSURE (slice scope)` when all values > 0.
- Full-space lift of the slice verdict via the accepted C-M5W-004 exclusion
  decomposition (only if exact; otherwise slice scope stands).

## Performance notes for the future run

- The Sturm chain per resultant is now cached (`_CHAIN_CACHE`); sign tests
  use an interval fast path over [A0, A1] with exact rational-Sturm
  fallback (`_SIGN_CACHE`); reduced field vectors (deg < 32) are tested for
  zero by emptiness, not by gcd.
- Remaining cost driver: isolation bisection calls `count_roots_field` per
  node; each node is 2 variation evaluations over a ~20-element field chain.
  If still slow, evaluate variations at ±∞ once and reuse across nodes
  (they do not depend on the interval), and bisect only interval endpoints.
