# 0139 verdict: STOP, dead-by-measurement (beacon)

Execution: ra_survey.py, exit 0, 149 samples, 33 s. Trust-r3,
E = 0.43 tube, 8 J-soft modes + 2 param dims + 32×2 interior.

## Measured (printed)

- Base 11.1572; central-du 20.3322 (banked decider).
- Axis min 6.46 (axis-m2+); all 20 axis points HOLD individually.
- Corners min 3.84 (corner-15), second 4.62 (corner-43).
- Interior min 7.56 → 7.29; sampling-doubling move 0.036 (< 0.10 ✓).
- SURVEY min = 3.8367 < 5.0 → STOP (frozen trigger fires).

## Reading (map texture, not just the verdict)

The tube is MIXED, and the mixture is informative: the Newton
direction is healthy (20.33), axis directions hold (≥ 6.46),
random interior holds (≥ 7.29) — only soft-mode+param CORNERS
break below 5 (3.84, 4.62). The revived survey's surprise (central
health) survives alongside the red-team expectation (corners kill):
both were right about different parts of the tube. Dead-by-
measurement stands per the frozen HOLD bar (min over tube), with
the map recording exactly WHERE: collective soft+param corners,
not the error direction, not generic interior.

## Scope

R-A dead-by-measurement (min-Q map banked above). No certificate,
no gap-protection claim. Instruments banked (sharp G1, C-sweep,
decider, survey). G-a2 H-side still BLOCKED.
