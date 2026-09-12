# 0122-beacon-fitted — free-boundary-fitted mesh rebuild (beacon)

Ordered: (a)-activated with quantified target (‖δF‖ 2.08 → ≲ 0.1, soft
share halved). Method: contour-fitted local red-refinement (not global
brute force — uneconomic at p≈0.4).

## Status: RUNNING (verdict on delivery)

- `fitted_mesh.py`: src=0 contour extraction (43 pts, core r∈[0.62,1.92])
  + 2-level band red-refinement: 861 → 1069 → 1451 nodes. Landed.
- `fitted_solve.py`: warm-started bordered Newton on fitted mesh.
- Gate: δF ≤ 0.1 → continue to tol; else next mechanism-grade stall
  with specified rung (stop rule per order).

## Also banked here

## Gate assessment (interim, pre-round-2)
GATE δF≤0.1 NOT met: fitted state gives soft share 0.2% (fixed from
50%!) but ‖δF‖=1.86, bound 22.1 (was 24.7) — marginal gain; a NEW
softer mode (0.0066) appeared. Round 2 (resume) running as bg_3.
IDEA-05 banked (soft-modes.npz + tail-mode reading + mesh rule).
- 0121 drift repairs: M2 probes archived + rerun-verified
  (`m2_gap_probe.py` reproduces transcript eigenvalues exactly);
  translation-invariance proof written (lemma.md appendix).
