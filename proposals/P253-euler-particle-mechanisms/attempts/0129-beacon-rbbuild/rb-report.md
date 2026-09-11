# 0129 R-B verdict: MISS -> skirt-obstruction record (beacon)

Experiment: rb_experiment.py, exit 0, trust-r3 tensor state.
Semi-Lagrangian front advection (43-pt contour by normal(du),
IDW extension 0065-fenced smooth, detJ-relative gate: zero drops).

## Measured (printed)

- 1.0x: shape-err 3.13 vs exact 9.18 (34%; bar <= 20%) vs linear
  3.97. MISS-PT (beats linear, misses bar).
- 0.5x: shape-err 1.70 vs exact 3.49 (49%) vs linear 2.60. MISS-PT.
- 0.25x: shape-err 1.04 vs exact 1.51 (69%) vs linear 1.30.
  Consistency bar (<= linear) PASSES.

## Adjudication (frozen trigger)

STOP FIRES (1.0x and 0.5x miss 20%). Front advection explains a
MAJORITY (~60-70%) of the exact response at every scale and beats
linearized transport everywhere — but a ~1/3 residual stands that
front motion cannot supply.

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
