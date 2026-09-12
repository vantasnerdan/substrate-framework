# 0122 fitted-mesh report: measured stall, line stopped (beacon)

## Runs

- Round 1 (bg_7, 692 s): interpolated warm start res 1.25e-2 →
  1.204e-2, κ̂=0.9991, rbar=1.0042, c=0.0096, outersrc 8e-32.
- Round 2 (bg_3, 734 s, resume): res → 1.132e-2, κ̂=0.9992,
  rbar=1.0040, μ=0.4610, c=0.0057, outersrc 1e-31.
- Stall signature identical both rounds and to coarse runs:
  STALL-NEWTON (|t0|≈0.26, detS≈−9.6) alternating GSTEP-OK micro-gains.

## R1 gate (drift R1 — fires ONLY on the eigen instrument)

`measure_deltaf.py` on the round-2 state: **δF = 1.84 (gate ≤ 0.1:
NOT met, ×18 miss)**; δλ-bound 21.8 vs 11.2. res-alone (1.25→1.13e-2)
would have laundered a continue — the instrument forbids it.

## Mechanism-grade reading

- Fitted mesh did its targeted job: soft share 50% → 0.23%.
- δF persists (2.08 → 1.86 → 1.84) via jf~1000 amplification of
  stiff-mode error + a NEW mesh-soft mode (0.0066). Refinement drained
  the channel it targeted; the remainder is NOT mesh-addressable.
- R6: tail-mode weight 55% — the remaining soft error is
  tail-dominated (supports a future tail-band rung, not core work).

## Stop (per 0122 stop rule)

Mesh line STOPPED. Specified next rungs (not executed here):
(a) tail-band extension, justified by R6 55%;
(b) jf-amplification route (p-continuation or lemma reformulation) —
now THE binding constraint;
(c) c-box watch (c 0.019 → 0.0057, bound 1e-3);
(d) identify the 0.0066 mode.

## Repairs landed

R1 instrument ✓ (this verdict), R2 hop-watch ✓ (no hop: κ̂≈1.0 both
rounds' heads), R3 mirror ✓ (code; next mesh build), R4 save guard ✓
(B.save_npz force=True, explicit chain), R5 design-0122.md ✓ frozen,
R6 ✓ (55% tail), IDEA-05 ✓ (modes + note).
