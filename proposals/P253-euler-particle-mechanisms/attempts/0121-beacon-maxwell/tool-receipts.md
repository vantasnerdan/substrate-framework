# Tool receipts — 0121 (beacon)

## M1 — gap probe (exit 0, decisive)

ARPACK `eigsh(SM)` k=6: No convergence, 8001 iters, 0/6 (5.96 s).
Dense `eigvalsh` (800 dofs): 1 negative (−0.506), soft cluster
0.016–0.075 (8 modes), bulk to 43.2 (1.3 s). Method note: SM-Krylov
is the wrong tool for indefinite clustered spectra; dense is exact
at this size — recorded for the (a) build.

## M2 — lemma computation (exit 0, with one scrapped attempt)

First attempt tangled δF-mapping scaffolding (broadcast error, killed
before any claim). Clean rerun: ‖δu‖_rms = 0.0472, ‖δF‖_2 = 2.0831,
δλ ≤ 24.66 vs λ = 11.16, margin 0.45, soft share 50%. All quantities
above are printed outputs, rounded here; the 4.72-vs-11.2 state-
sensitivity note stands separately.

## M3 — graph query on fresh index (exit 0)

`gitnexus query "Newton solver..." -r substrate-framework`:
surfaced `radial_harmonic_balance.py` (collocation residual,
declared branch coordinates, channel-aware outer policy). Reviewed:
patterns noted, no code reused (different domain). ripwire map on
0117 dir returns legend-only output (tool immaturity on my side,
recorded; grep covered the no-prior-solver finding across 199 files).
