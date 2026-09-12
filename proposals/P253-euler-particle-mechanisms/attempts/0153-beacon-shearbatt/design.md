# 0153-beacon-shearbatt - C3 shear-scan battery DESIGN (beacon)

Charter: shepherd C3 battery (drift redirect; executable with
banked parts). Design frozen BEFORE compute (this file). M4-genus
flag travels (background shear imposed+labeled; 0043 untouched).

## Construction (all parts banked)

- Shear rhs: rhs_s(X, S) = rhs3(X) + S·X[...,0]·ê_z (planar
  shear, genuine strain S everywhere, no rigid-carry degeneracy;
  imposed+labeled background).
- Re-shoot under shear: local copy of shoot/res3/flow_frac with
  rhs_s (A3 shoot takes no kw — copy noted, not edited);
  warm-start chained upward from S = 0 banked orbit.
- m-Floquet per S: local monodromy copy (basis_field/project/
  dir_index/ndirs/axisym_of reused read-only), m = 1..2,
  growth(S) = max |multiplier|.
- Available shear S_av: trust-r3 velocity via 0077 identification
  W = −(P_z/r)ê_r + (P_r/r)ê_z with P := u (c-shift drops out of
  strain); S_av = max Frobenius sym-grad-U over r ≥ 0.15.
  Scale fallback on record: U/L ~ 1.3/1 if reconstruction
  disputed (order-unity either way; adjudication uses robust
  thresholds below).

## Falsifier (frozen)

S* = smallest S in {0, 0.05, 0.1, 0.2, 0.4, 0.8} (one refinement
round allowed around an interior transition) with growth ≤ 1.05.
- Baseline growth(S = 0) ≤ 1.05 → S* = 0 → C3 premise VACUOUS
  (nothing to stabilize; report moot, not alive).
- S* > S_av (with S_av < S*/3 robust margin) → DEAD.
- S* ≤ S_av → ALIVE-conditional (background-supplied label
  permanent; never derived stiffness).
- S* inside [S_av/3, 3·S_av] → GRAY (measure better, no verdict).

## Stop + cost

Miss → C3 dead with numbers. Gray → one refinement (S grid +
S_av mask), then verdict. Cost: ~6 shoots + ~12 monodromies
(N = 64), est. < 15 min. No A3 edits (local copies in driver).
No promotion. No import (shear is a test background, not physics).
