# Attempt 0013 — N2/m=1: true slow branch found; LIA-branch pose gap named

## Correction of the 0012 diagnosis

The m=1 truncated system DOES carry a true slow branch — the 0012
"non-closure" read F minima at the LIA scale because the branch is NOT of
the LIA form. Full-range residual scans (k = 1e-3):

  |r1| and |r2| share ONE global minimum at

    omega = 0.0003522958...  (|r1| = |r2| = 2.2e-6),

  linear in k. Decade extrapolation of omega/k:
  0.3325 (k=0.02), 0.3413 (0.01), 0.3523 (0.001) ->

    ** omega = (sqrt(2)/4) Om (ka) + O((ka)^2 ...) **  [NUMERIC_EVIDENCE]

  Residuals at the branch point -> 0 with k (|r| ~ 5e-4 at k=0.01,
  2.2e-6 at k=0.001). This is a genuine Coriolis-inertial branch of the
  corrected (2 Om) system, in the wt -> -Om family (same family as the
  m=2 polarization wave), with a static-displacement continuation
  (pattern speed -> 0, the translation equilibrium).

## Named missing construction (why the textbook LIA k^2 ln branch is absent)

The classical bending branch omega ~ Gamma k^2 ln(2/(ka))/(4 pi) is NOT a
zero of the current 4-condition set (F at the LIA scale ~ 7e-3, not a local
minimum). The frozen-vorticity sheet-strength transport was imposed in its
untilted form [v_th] = +2 Om eta; for the m=1 HELICAL sheet the strength
vector tilts (O(k eta) transverse components), and the tilt correction of
the transport enters the dispersion at exactly leading order O(k^2 ln) for
m=1 (while for m=2 it is higher order, which is why the m=2 channel closed
exactly). Next construction: tilt-corrected sheet-strength transport ->
the LIA branch and its constant (the c_inner analog) on that branch.

Also owed: classical corroboration of the linear-inertial branch and of
the tilt-corrected LIA constant for the Rankine core (Kelvin 1880 /
Saffman-level source), per the cited-not-recalled discipline.

## Status

- m=2 twist channel: COMPLETE (unchanged; receipt-grade).
- m=1 bend channel: slow-branch structure established numerically
  (linear-inertial; coefficient sqrt(2)/4); LIA branch blocked on the
  named tilt-corrected transport construction.
- N2 remains active; verify_cst002 deferred until the m=1 branch set is
  complete, so the frozen verifier covers BOTH branches and the m=2
  channel with mutations m1-m6 plus (m7) the untilted transport must
  FAIL to produce the k^2 ln branch.

## Sharpening (same attempt, post-record check)

Re-examination of the omega = 0.35k structure: its residuals at the F
minimum are 2.2e-6 at k=1e-3 — SMALL BUT NONZERO, decaying as a power
(~k^2.4), and no exact common zero of r1 and r2 exists anywhere on the
scanned domain (1e-9 < w < 2). Both the linear-inertial structure and the
LIA-scale structure are therefore ASYMPTOTIC-ONLY in the current pose:
the linear problem as matched (kin pair + pressure continuity with
advected Bernoulli + frozen tangential condition) is STILL MISSING a
sector — the m=1 eigenbranch has no exact realization in it, while the
m=2 channel closes exactly (F ~ 4e-29).

Candidate missing sectors, to be adjudicated by 0014 in order:
1. The O(k eta) tilt correction of the sheet-strength transport
   (carried over from the 0013 diagnosis).
2. The exterior homogeneous sector at k != 0 beyond the decaying
   K_m harmonic (the uniform-translation sector exists only at k = 0;
   something must carry the impulse structure of the displaced
   filament at k != 0).
3. The interior regularity structure at the omega -> 0 degeneracy
   (mu -> infinity-limit of the Poincare coefficient: the r-linear
   pressure is approached only logarithmically; the matched interior
   may require the finite-mu treatment rather than either limit).

Route verdict for the m=1 finite-k pose: blocked, missing construction
named at the sector level (1)-(3). The m=2 channel and the receipt are
unaffected. N2 remains active.
