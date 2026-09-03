# Attempt 0004 — G3: global nonnegativity of V_omega* (omega_c^2 = omega_*^2)

## Typed chain (analytic receipt, frozen before any output)

- mathematical object: the rotating-frame slice potential
  V_w(m,c,b,f) = V0 - (w^2/2)(f^2 + 4b^2), a real quartic polynomial on
  R^4 (the aligned real-psi slice S = diag(m, c+b, c-b), psi = f >= 0
  real part; the imaginary part enters V0 only through positive terms per
  accepted C-M5W-004 exclusions).
- symmetry/conservation license: V_w is S1-invariant (orbit invariance,
  accepted); the Maxwell point w_*^2 is the unique deep-branch frequency
  with V_w(B) = V_w(A) = 0, A = (1,0,0,0), B = (0,c*,b*,f*) certified
  (attempt 0001, 55-digit coordinates, Krawczyk enclosure of w_*^2 of
  width 1.2e-43).
- exact variational functional: V_w itself (no ensemble; a pointwise
  polynomial inequality).
- admissible space/representation: the exact polynomial in four real
  variables with coefficient field Q(w_*^2); gradient relations computed
  in SymPy from the canonical module
  `substrate_framework.m5_wall_clock.wall_slice_potential` (never
  hand-coded).
- analytic scale/asymptotic structure: leading homogeneous part
  (m^2+2c^2+2b^2)^2 is positive definite, so V_w -> +inf uniformly at
  |u| -> inf; the global infimum of V_w is attained at a critical point
  of V_w (compactness of the negative sublevel sets).
- observable: min_u V_{w_*}(u); target statement: min = 0, attained
  exactly at {A, B}. Equivalently omega_c^2 = 2 inf V0/iota = w_*^2
  (using V0(B) = (w_*^2/2) iota(B) and V0 >= V_{w_*} pointwise).
- irreducible remainder: the sign of finitely many exact algebraic
  numbers (critical values of V_w). This is exact algebra (Grobner /
  resultant / Sturm / interval sign certificates), NOT production
  numerics; no discretization enters.
- numerical approximation: interval-arithmetic sign evaluation at the
  certified w_*^2 enclosure (Rigor via mpmath iv with guaranteed
  outward rounding); fallback: exact rational reconstruction +
  Sturm sequences if a value sits within the interval width of zero.
- permitted verdict: OMEGAC_CLOSURE (slice scope): V_{w_*} >= 0 on the
  slice with zeros exactly {A, B}, hence omega_c^2 = w_*^2 on the slice;
  full-space lift only via the accepted exclusion decomposition (checked
  exactly; if the lift is exact the claim upgrades to full space, else
  the lift is stated as the accepted C-M5W-004 exclusion scope).
- requires: C-M5C-001..004 (canonical action, exclusions), attempt-0001
  Maxwell certificate.
- pass_licenses: L-S1-OMEGAC (slice; full space if lift exact).
- does_not_license: dynamic stability (G4), wall tension values, bag
  existence, any statement off the accepted action.
- maximum_verdict: exact proof (symbolic) on the slice (+full space if
  the exclusion lift closes exactly).
- failure_scope: if a critical point with V < 0 is found, then
  omega_c^2 < w_*^2 strictly and C-M5W-004's exclusion upgrade is
  refuted on the slice (major canon-relevant finding; record and
  continue per continuation ladder).
- unlocks: G4 (stability about a provably nonnegative-clock-potential
  background), the interpretive claim that the shell-bubble sits at the
  true Maxwell frequency.

## Route

1. Exact gradient structure: derive and verify in SymPy that
   grad V_w = 0 in variables (m, s, d, f) = (m, 2c, 2b, f) reads
     (1) m (4q - 1 - 3m) = 0,        q = m^2 + (s^2 + d^2)/2
     (E1) d (4q + 1 - 3s - 2 w^2) = 0
     (E2) s (4q + 1) = 3 (q - m^2)
     (4) f (w_of(f,b) - w^2) = 0,    w_of = 32 f^2 - 12 f + 6 - 24 b
   This is a small case enumeration, not a generic Grobner blow-up.
2. Enumerate all cases exactly (m = 0 vs m = (4q-1)/3; d = 0 vs
   4q+1-3s = 2w^2; f = 0 vs w_of = w^2), solve each finite subsystem,
   and evaluate V_w at every real solution with w^2 = w_*^2.
3. Sign certificates: interval arithmetic (mpmath iv, dps scaled,
   guaranteed outward rounding) at the certified w_*^2 enclosure; any
   value within 1e-30 of zero is resolved exactly instead (substitution
   via the Maxwell relations; Sturm if needed).
4. Zeros: verify V_w = 0 exactly at A (trivial) and B (from the Maxwell
   system pA = pB = pC = 0, w^2 = w_of).
5. Full-space lift: compose V_full - V_slice from the canonical module
   and check exact nonnegativity of the difference (the C-M5W-004
   excluded directions); record the result as earned.

Non-goals: no SDP/numerical SOS relaxations (unverified heuristics) are
used as evidence; no claim about sectors outside the accepted exclusion
decomposition if the lift does not close.
