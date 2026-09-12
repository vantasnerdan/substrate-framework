# drift PRE-review — beacon 0122 fitted-mesh (worktree bytes, build running)

Scope: early review per routing (target-gating, δF measurability, stop rule, spend-wasting
flaws). Targets uncommitted worktree files (README#36A6, idea05#B630, fitted_mesh#6CA2,
fitted_solve#F0E8, bank_modes#C67C); re-verify on landing if bytes move. Verdict: CONDITIONAL
GO — do not stop the running build (resume-capable), fix R1/R2/R4 immediately, R3/R5 today.

## Sound (do not change)

Fitted-over-brute at p≈0.4 ✓ (correct economics); trust-state warm start ✓; banked soft modes
with deflation-FORBIDDEN-for-QOI (preconditioner-only with unpreconditioned verdict metric) —
exemplary restraint ✓; tail-band next move pre-specified with trigger (soft share persisting) ✓;
RESUME chain ✓; stop rule (gate else mechanism-grade stall) ✓.

## Repairs (ranked)

- R1 (HIGHEST — gate measurability): NOTHING in fitted_solve.py computes δF (prints
  res/kap/rbar/mu/c only) yet README gates on δF≤0.1. res ≠ δF (0121: ×15–60 sensitivity
  amplification). Without the lemma re-measure step the gate cannot fire and continue-to-tol
  would launder on res. Add the eigen-pipeline re-measure on the fitted state (dense feasible
  at 1451 dofs) BEFORE any gate evaluation.
- R2 (HIGH — basin risk on mesh transfer): interpolating u0 onto a new mesh inside a NARROW
  basin (IDEA-03 FAILED-narrow) risks a 0119-style hop on iteration one. Monitor κ̂/rows for
  the known hop signature (κ̂→3–10, umax jump) on first fitted iters with abort criterion.
- R3 (MEDIUM — asymmetric band): contour extraction skips z<0 (fitted_mesh.py:26), so the
  lower half gets no refinement on an even-symmetric problem — half the budget wasted,
  possible bias. Mirror cpts z→−z (or justify the asymmetry).
- R4 (MEDIUM — repeat lesson): member-fitted.npz via raw np.savez (fitted_solve.py:50),
  bypassing the save_npz/--force/--out-npz guard built after 0114+0120. Route saves through
  it or version per-rung.
- R5 (MEDIUM — process): no frozen design-0122.md exists though production runs (band
  widths 0.30/0.15, 2 levels, gate, warm-start all live in code/README). Freeze it dated NOW;
  small-ratio-numerics requires design-before-compute and the spend is burning.
- R6 (forward): use banked core-weights to decompose soft-share core-vs-tail BEFORE further
  refinement — if tail modes dominate δF, tail-band (pre-specified) should precede more core
  work. Data already banked; one analysis step.

## Verdict

CONDITIONAL GO with R1/R2/R4 immediate, R3/R5 today, R6 before next mesh move. Method
direction correct; all findings are missing-guardrail class, none is a wrong-turn. Re-verify
at landing commit.
