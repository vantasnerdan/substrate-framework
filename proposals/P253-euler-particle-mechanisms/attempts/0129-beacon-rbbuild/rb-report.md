# 0129 R-B verdict: MISS -> skirt-obstruction record (beacon)

Experiment: rb_experiment.py, exit 0, trust-r3 tensor state.
Semi-Lagrangian front advection (43-pt contour by normal(du),
## Measured (printed; CORRECTED per drift review-beacon-rbstop)

CORRECTION: the report v1 compared shape-err against the linear
PREDICTION (5.21/2.60/1.30). Design G-Q5 bars use linear-ERR
(|pred-exact| = 3.97/0.89/0.21). Honest recomputation:
- 1.0x: shape-err 3.13 vs exact 9.18 (34%; bar <= 20%) vs linear-err
  3.97. MISS-PT (beats linear — ONLY scale where shape wins).
- 0.5x: 1.70 vs 3.49 (49%) vs linear-err 0.89. MISS-PT (LOSES).
- 0.25x: 1.04 vs 1.51 (69%) vs linear-err 0.21. Consistency bar
  FAILS on design terms.

## Adjudication (frozen trigger)

STOP FIRES (1.0x and 0.5x miss 20%; quarter-scale misses too).
Explained shares ((exact-shape)/exact): 66%/51%/31% per scale —
falling trend supports the mechanism (front dominates large
displacement; bulk-linear wins small). Residual stands that front
motion cannot supply; every correction points the stop harder.

## Skirt-obstruction record (amended fallback; R-C motive)

The residual's address: s^6 skirt-curvature + bulk response outside
the advected front band. A coercivity argument (R-C) must cover the
skirt, not just the front: min-eig(B) >= 9.0 threshold stands, now
with the quantitative note that ~6 of the 9 units are
front-motion-adjacent (explained) and ~3 are skirt (unexplained).
Caveat (reopenable): ONE extension choice (IDW+cutoff); a challenger
may rerun the frozen trigger under a better extension.

## Entry-gate status

G-Q1/Q3/Q4/Q5/Q8/Q9 + G-0071 all held through execution (no
continuity used, no domain linearized, verdicts conditional-scale,
exact outranked, chart never needed: tags advected passively).
