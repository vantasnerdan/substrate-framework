# D3b F-bar freeze (pre-compute; D-08 form) — adiabatic-Berry closed loop

Construction: Γ schedule closed loop 1.0 → 1.01 → 1.0 over 10 periods
(there-and-back; max excursion INSIDE N4-proven-smooth zone ≤1%; N=64).
Per period p: re-shoot (damped, banked/continuation hot start) → section-flow
3x3 FD → transverse eigvec; overlap-chain across periods (Berry transport,
genuinely nontrivial: base orbits differ); per-period Floquet phase φ_p.
Observable: geometric estimator γ_geom = Γ_tot − Σφ_p (total chain phase
minus summed per-period dynamical); floor ~1e-3 (overlap resolution × links).
- ALIVE-leaning: |γ_geom| > 0.05 AND eps-halving leg agrees within 2× →
  geometric accumulation EXISTS (P2-adjacent leg opens real).
- DEAD: |γ_geom| ≤ 0.01 → no detectable geometric accumulation.
- Gray otherwise → UNRESOLVED + named leg (slower ramp, N=128).
- STOP (frozen, N4 fold warning carried — ENTRY CONDITION): per-period
  branch guard (re-shoot res < 1e-8 AND radii within 2× banked AND |ΔT|/T <
  1% adiabaticity) — violation → STOP + record fold/break location (maps
  smooth-patch boundary; informative, not failure).
- INVALID (no verdict): overlap link < 0.99 (level crossing — re-examine).
Run reports (γ_tot, Σφ_p, γ_geom, guard log, verdict), nothing else.
P2-adjacency scope (pre-computed): existence only; dynamical charge still
missing-5; no promotion.
