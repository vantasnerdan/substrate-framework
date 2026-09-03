# Attempt 0002 — G1: sigma_0 value and the wall profile

Status: established (numeric evidence, SIGMA0_SLICE_VALUE scope)
Parent: P250 gap G1 (owner-directed continuation of the accepted
C-M5W-001..005 campaign; same branch/PR #196)
Route: analytic receipt (receipt.md, with the production-route amendment)
-> landscape check -> L-continuation collocation -> crossed refinement ->
independent BC and quadrature cross-checks -> itemized budget.

## Input data (certified)

- omega_*^2 = 1.6639457000591502988561930002961614442197723 (55-digit
  certified deep-branch Maxwell crossing; interval enclosure width 1.2e-43).
- Bulk states: A = (1,0,0,0); B = (0, c*, b*, f*) with the 55-digit
  certified coordinates; Krawczyk uniqueness box 1e-12.
- Exact potential/EL primitives imported from
  src/substrate_framework/m5_wall_clock.py (no restatement); first-integral
  identity re-verified exactly at run start.

## Routes tried

1. Direct long-box collocation (Dirichlet or exponential-fit BC, tanh /
   linear / step / string initial guesses): FAILS into multi-kink trains —
   Newton basin structure, documented as basin probes (status=1,
   non-monotone, r1 in {0.93, 2.19} garbage).
2. Hand-rolled field-space string method (normal descent): stalls on a
   discretization shelf; sigma drifts in a ~0.5%-wide flat valley
   (0.743..0.749). NOT part of the budget; used only as the L=2 initial
   guess and an order-of-magnitude cross-check (it brackets 0.7293 from
   above, as it must: the string sigma overestimates the kink tension).
3. NEB with spring+tangent-switching forces: DIVERGES here (overflow);
   abandoned — recorded as a route verdict, not hidden.
4. L-CONTINUATION (production, adopted): short boxes admit only the
   1-kink; warm-started extension L = 2, 3, 4.5, 6, 8, 10, 12 with tol
   1e-6..1e-10 tracks the monotone branch to L_FIX. Converged cleanly.

## Result

sigma_0 (slice-sector wall tension at omega*^2) = 0.72929841787 with
itemized error budget below; the wall profile u(x) = (m,c,b,f)(x) is
constructed (profile_L12_n3201.csv, monotone, first integral conserved).

Budget (finest L=12 rung; acceptance threshold 7.3e-4):

- b1 collocation residual (scaled 2L): 2.0e-7  [dominant, conservative]
- b2 h-truncation half-spread:           5.3e-12
- b3 domain (|sigma(L=12)-sigma(L=10)|): 1.5e-11
- b4 certified-box drift (+/-1e-12):     2.5e-11
- b5 quadrature (GK vs compensated trap):9.1e-12
- b6 evaluator noise:                    4.4e-16 (fsum permutation exact 0)
- b7 boundary-treatment agreement:       1.1e-13 (Dirichlet vs expfit)
- b8 first-integral drift:               2.0e-11 (4.0e-15 at the clean rung)
- total: 2.0e-7 -> PASS (0.00003% relative)

Route agreement: four tension routes (T+V, 2V, 2T, Gauss-Kronrod through
the collocation spline) agree to 1.6e-11; the clean rung (L=12, n=1601,
status=0) has residual RMS 2.4e-11 and fi-drift 4.1e-14.

Landscape facts established on the way (recorded for G3):
V_w >= 0 on the sampled A->B corridor at omega*^2 (interior min 0.015 at
m=0.1 over a 60^3 grid in (c,b,f) x 21 m-slices); mountain-pass height
max V_w on the connecting path ~ 0.309; tail masses at A: 2.082 (softest,
f-direction), 3.162, 3.162, 3.917; at B: 2.528, 2.899, 3.277, 7.012.

## Reproduce

PYTHONPATH=src python proposals/P250-shell-bubble-clock/attempts/0002/sigma0_wall.py
(threads pinned to 1 at import; deterministic; ~13 min)

## Does not license

Full-space (off-slice) minimality; bag physics (G2); stability (G4); no
claim promotion by itself. Claim candidate C-M5W-006 (numeric tension +
constructed profile) formalizes this result for individual review.
