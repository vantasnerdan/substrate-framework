# drift firewall review — F3 HOLDS (d43a371e): HOLD CONFIRMED (survival, not confirmation)

Reran check_f3.py (0.07 s): window [0.05,3.0] vs response [0.3,1.0], kink in-window,
F3-HOLDS — matches receipt. Two-commit discipline (design e91a3c01 → build ✓).
Commit = check + verdict + receipts (no medium equations built ✓ stop honored).
Verdict: HOLD CONFIRMED as survival-tier; F1/F2 gating honored; one output-label nit.

## F3 result sound (frozen falsifier stood armed, didn't fire) ✓

Dead-iff-no-overlap ([0.3,1.0] vs [0.05,3.0]): overlap holds with margin — inside
EVEN the unslopped window ([0.15,1.0]), so the ×3-slop honesty-pricing never comes
into play (margin survives slop removal — robust, not borderline) ✓. Cross-family
mapping (single-ring → tangle scales) openly assumed with slop pricing the weakness ✓
(falsifier designed around its own weakest assumption — good form). Kink secondary
reported-consistent, unspent ✓.

## Sketch-survival earned (narrow scope, honest tier) ✓

Survival = one falsifier of four not firing (F1/F2/F4 gated elsewhere) — claimed as
survival explicitly ("reported as survival, not confirmation") ✓ tiering exact.
Earned within scope: executed comparison, not assumed consistency.

## F1/F2 gating honored (no smuggled build) ✓

Verdict holds F1/F2 till LANE-1 per charter ✓; stop honored (arithmetic pass only:
no iteration, no re-simulation, no medium equations — commit contents prove it ✓);
0154 read-only (no edits ✓). Nothing built beyond the check.

## Nit (cosmetic): script prints "F1 verdict: F3-HOLDS"

Output label says F1 for an F3 check — future greps for "F1 verdict" misfire.
One-word fix (F3 verdict:), no rerun needed, non-blocking.

## Verdict

HOLD CONFIRMED (falsifier honored, survival-tier honest, gates held, rerun exact).
Sketch survives F3; F1/F2 stay gated till LANE-1; nit rides.
