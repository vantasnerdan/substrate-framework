# drift firewall review — cipher 0120 M2-B1 exploratory (9da41e54, A1+A2; A3 queued)

Rule firewalled: attempt=result with banked mechanism; weakest-first order; no proof demanded.
Both runs reproduced digit-for-digit by drift (A1 8.23/9.20/0.000; A2 ell 1.0166/1.011 +
0.986a).

## A1 (H5 exclusion) — PASS-in-model

Construction minimal and exact (orbit integration + d(t)/a series on reproduced rk4).
min 8.23, never <4 ✓; "pending Euler orbit" correctly bounds the transfer claim (distance
is orbit property; live version needs Euler orbit). Precision note (non-fatal): "margin
widens ∝a⁻¹ (2× at a=0.025)" is ANALYTIC scaling at fixed orbit, not a run — label it so;
no receipt claims otherwise, but a reader could misread it as measured.

## A2 (H1 rigidity) — repaired PASS-in-model; repair LEGITIMATE, failure banking COMPLETE

Attempt-1 failure correctly diagnosed as METHOD (regularized-BS interior void for r<a —
the model has no core solution there), not physics; ell 4.2/16.3 discarded with reason. ✓
Rankine graft principled, not ad-hoc: solid rotation Ω=Γ/2πa² via EXACT matrix split
(explicit-Euler 8e13 blowup banked as integrator trap), ring translation from model,
other-ring strain differential (m−base) with self-singular excluded — the correct
rigid-interior + external-strain decomposition. Result 1.0166/1.011 vs O(1.05) expectation:
same ballpark, credible. Observation (not repair): mean radius 0.986a = 2.8% area drift,
same order as the signal — ellipticity ratio is scale-invariant so the reading stands, but
the drift is consistent with discrete differential advection, worth one line in the receipt.
Banking complete: void-interior + broadcast + blowup + paren-typo all recorded with lessons;
dye cousin-caveat attached (interior untested). Remaining gaps stated (graft≠Euler core;
multi-period secular open per 0113 item 2). ✓

## Scoping honesty + order — HELD

Weakest-first honored; "PASS-in-model at Rankine scope" never exceeds evidence; H2/H3
blocked-scoped with revisit condition (G-a2/solver); A3 queued, combo still gated, R-EM2
owner-side restated. No proof demanded or smuggled.

## Verdict

PASS A1 + repaired-PASS A2 as exploratory construction (one labeling note on a⁻¹ scaling,
one area-drift observation). M2-B1 advances two rung-ladder steps in-model; live-field H1/H2
and 3D A3 remain the open constructions.
