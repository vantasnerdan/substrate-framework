# 0156 F3-only build DESIGN (beacon) — frozen pre-compute

Charter: shepherd F3-only (790812aa). F1/F2 HELD till LANE-1
fires (hold condition, not drop). Sketch 0155 §3-F3 under test.

## F3 question

Is the banked backbone response localization consistent with the
F-A μ-medium's spatial window [a, ℓ]? Clash → sketch DEAD; overlap
→ F3 HOLDS (survives, not alive — F1/F2 still gated).

## Banked inputs (read-only, verdicts untouched)

- 0151: response lives in skirt annulus 0.3 < |s| < 1.0 (ring
  units R = 1); kink band ~0.2x. Challenger-report numbers.
- 0153: S_av = 19.67 strain-rate field (amplitude context only;
  rates don't enter the spatial check).
- 0117 member: R_ring = 1, ε_core = 0.15 (frozen design params).

## Stated scale mapping (assumption, honesty-barred)

F-A's tangle scales mapped onto ring units: a ↔ ε_core = 0.15
(core radius), ℓ ↔ R = 1 (single-ring outer scale). Cross-family
(single-ring maps vs tangle medium) — stated, not derived; the
×3-slop rule below prices exactly this weakness.

## Frozen falsifier F3

F3 FIRES (sketch dead) iff NO overlap between response interval
[0.3, 1.0] and window [a/3, 3ℓ] (×3 slop each end absorbs the
mapping weakness): i.e. dead iff 1.0 < a/3 OR 0.3 > 3ℓ. Else F3
HOLDS. Kink secondary: kink 0.2x vs a/3 — reported, not decisive.

## Frozen stop

One arithmetic pass over banked numbers (assert script, exit 0).
No iteration, no re-simulation, no medium equations (gated).
Verdict: F3-FIRES or F3-HOLDS. Then report.
