# Tool receipts — 0120 (beacon)

## T1 — trust rounds (exits 0; STALL/GSTEP verdicts inline)

- R9 recovery (--single-p 6, coarse): 0.735 → 2.02e-2, 4 iters;
  kap=1.018, rbar=0.996, c=+0.066, outersrc ~1e-35. Deterministic.
- Trust round 1 (from R9): → 1.17e-2.
- Trust round 2: → 8.48e-3 (kap=1.0002, rbar=1.0055).
- Trust round 3 (bg_1 at writing): target ~6e-3.
- Every Newton stall broken by GSTEP-OK at least once per round;
  c stays +0.017..0.021 (box never binds on-branch); outersrc
  1e-33..1e-35 throughout (compactness monitored, clean).

## T2 — refinement diagnostic (exit 0, 282 s)

Trust state interpolated 40×20 → 80×40 → 120×60, residual
re-evaluated (no solve): 5.85e-3 → 6.13e-3 → 3.79e-3. Observed
order ≈ 0.4 (unfitted free boundary); brute force uneconomic.

## T3 — trust-state feed (exit 0, 0.65 s)

λ_ω = 11.12635817 (exact doublet) + axial 35.01; H-rows only,
EXPLORATORY with 6.2e-3 floor noted at print time.

## T4 — basin probe (exit 0)

18%-perturbed restart escapes (lands kap=0.89/res 0.74, c pinned
at bound — box observed working). Basin narrow; R9 not global.
Overwrote member-trust via --force (own bypass, documented);
regenerated deterministically (this receipt chain).

## T5 — graph re-index + query (exits 0)

`gitnexus analyze`: 57,352 nodes / 86,178 edges. Query surfaced
`radial_harmonic_balance.py` (collocation discipline, declared
branch coordinates, channel-aware outer policy): patterns noted,
no code reused (different domain) — recorded to close the user-
asked tool loop honestly.
