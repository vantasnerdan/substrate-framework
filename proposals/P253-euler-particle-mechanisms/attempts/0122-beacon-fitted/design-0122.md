# 0122 frozen design: contour-fitted rebuild (beacon)

Status: FROZEN 2026-09-10 (drift R5 — written during bg_3 run; governs
all 0122 production solves from here; earlier exploratory runs labeled).

## Object

Bordered (κ, rbar) Newton (p=6, reg=1e-3) on a contour-fitted mesh,
warm-started from the trust state, with c≥0 box + gradient fallback
+ outer-source monitor (all in 0117 build_member.py, unchanged).

## Frozen freedoms

- Contour: src=0 level set of trust-r3 state, 72-ray polar bisection
  around (1,0), star-shaped assumption (measured: 43 pts, core
  r∈[0.62,1.92], |z|≤1.27).
- Band: widths 0.30 (level 1) + 0.15 (level 2), red-refinement ×2
  (861 → 1069 → 1451 nodes). MIRRORED z→−z (drift R3 — even problem).
- Box/BCs/measure: unchanged from 0117 design (6×3, Dirichlet outer,
  r-measure forms, OMP pinned, versions recorded).
- Gate: δF ≤ 0.1 re-measured by the eigen-pipeline on the fitted state
  (drift R1 — res NEVER gates; 0121 showed ×15–60 amplification).
- Hop-watch (drift R2): first fitted iters must show κ̂≈1/umax≈1.3;
  κ̂→3–10 or umax jump = 0119 hop signature → ABORT, do not burn.
- Saves: through save_npz guard (drift R4); per-rung versioned names.
- Stop rule: gate met → continue to tol; mechanism-grade stall →
  stop with specified rung (no blind variants).

## Out of scope (this attempt)

Tail-band extension (trigger: soft share persisting with tail-mode
dominance per R6 analysis); Maxwell stage (after member tol-level);
production feed (after gate).
