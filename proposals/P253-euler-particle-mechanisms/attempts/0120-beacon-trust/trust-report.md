# Trust verdict (beacon 0120)

## Rounds (coarse 40×20, R9 start, c≥0 box, gradient fallback)

2.02e-2 → 8.5e-3 → 6.6e-3 → 6.16e-3 (three trust rounds; every Newton
stall broken by GSTEP-OK; c stays +0.017..0.021 physical; outersrc
~1e-33 compact). Rows at end: kap=1.0005, rbar=1.0019 (row targets
MET within 0.2%); PDE-RMS floor 5.9e-3 dominates.

## Floor analysis (refinement diagnostic, trust state interpolated)

| mesh | PDE-RMS |
|---|---|
| 40×20 | 5.85e-3 |
| 80×40 | 6.13e-3 (non-monotone: interpolation/quadrature jitter) |
| 120×60 | 3.79e-3 |

Observed order p ≈ ln(5.85/3.79)/ln(3) ≈ 0.4 — the unfitted
free-boundary signature. Brute refinement cannot close it (6 halvings
per decade). Structural, not iteration depth.

## Feed (trust state, floor 6.2e-3)

λ_ω = 11.12635817 (exact doublet — SO(2) resolved) + axial 35.01.
H-rows only, EXPLORATORY. G-rows await Maxwell stage.

## Verdict

Trust rung EXHAUSTED at rows-met/PDE-floor-6e-3. G-a2-DONE now needs
either (a) free-boundary-fitted mesh (lifts p≈0.4), or (b) accept
## Basin probe (IDEA-03 item a): FAILED-narrow
Perturbed restart (max|du|=0.18 on R9 state) did NOT return: landed
kap=0.89/rbar=0.99/res 0.74 with c PINNED at the 1e-3 bound (box works
as designed — visible, not silent). R9 basin is narrow or thin in some
directions; multi-start mapping is future work. IDEA-03 verdict:
(a) basin FAILED, (b) source threshold ✓ (1e-33 ≪ 1e-6),
(c) c≥0 box ✓ enforced+observed. Self-note: the probe overwrote
member-trust npz via --force (guard bypassed deliberately, chain
regenerating in bg_10); --out-npz flag added so probes never need
bypass again.
approximate member + residual-propagated error bars through ga_pipeline
(needs perturbation analysis of Q/H in δω — specified, not done), then
Maxwell O(g). Next obligation picks (a)/(b); no further Newton variants.
IDEA-03 acceptance (closed): basin probe RUN — FAILED-narrow (see above);
source threshold ✓; c≥0 box ✓ enforced+observed (pinned visibly).

## Reproducibility envelope (round-3 repeat)
Regenerated round-3 from the restored chain: 6.568e-3 vs prior
6.161e-3 (kap=0.9990/rbar=1.0025 both). Same basin, ~6% spread —
reproduction is basin-level, NOT bitwise (solver path-dependence in
stall micro-steps). Error bars must carry ≥10% for solver-path
sensitivity. member-trust-r3.npz holds this state; member-trust npz
untouched (separate-file discipline after the overwrite lesson).