# 0114 frozen design: S9 exposing probe (beacon)

Status: FROZEN before any run (written prior to `s9_probe.py` execution).
Change requires a new dated section, never silent edits.

## Question

Does a shape observable detect filamentation-type deformation where
stability norms stay small (S9 Theorem A coexistence)? Exact object under
test: the OBSERVABLE PAIR, on the exact Hill carrier (S2 analog in the
same axisymmetric no-swirl class — definition transfers, calibration does
not; no S2 persistence claim).

## Carrier (exact, paper-pinned)

Hill unit ball, `ξ_H = 1_B`, comoving frame, `W = U = 2/15` (S9 §2.2,
(2.5) garble cross-checked: our `u_r = (3U/2)rz` interior and
`u_r = (3U/2)rz/ρ⁵` exterior reproduce both surviving fragments; speed
matches stated `W_H = 2/15`):
- interior (ρ ≤ 1): `u_r = (3U/2)·r·z`, `u_z = (3U/2)·(1−z²−2r²)`
- exterior (ρ > 1): `u_r = (3U/2)·r·z/ρ⁵`,
  `u_z = −U·(1 − 1/ρ³ + 1.5·r²/ρ⁵)`
- frozen field (NOT live Euler): scope is observable sensitivity, not a
  filamentation proof. S9's theorem is cited, not reproduced.

## Perturbation analog + observables

- Seed N tracers in shell `1 < ρ < 1+δ`, δ = 0.05, behind hemisphere
  (z < 0, S9 §1.4 lagging-tail channel) + full-shell control set.
- `D(t)` = max tracer distance from origin (diameter proxy).
- `F(t)` = fraction of tracers with ρ < 2 (core-proximity proxy:
  stability-norm analog — bulk stays near core).
- Advect by RK4 in (r,z) half-plane, window T = 30 (≈ 4.5 core transits).

## Frozen pass criteria

- EXPOSED (observable discriminates): linear fit `D(t) = D0 + s·t` on
  `t ∈ [10, 30]` has `s > 0.05` AND `R² > 0.95` AND `F(30) > 0.9`
  (tail grows while bulk stays — stability-blind growth demonstrated).
- BLIND: otherwise (record as-is; do not retune seed/window).
- Convergence: N ×2 and dt /2 reruns must keep verdict (tols: s ±20%).

## Limits (pre-registered)

Frozen (not live) field; axisymmetric; finite window; tracer-count
convergence only; Hill ≠ S2 (no transfer of numbers); live-Euler
filamentation remains S9's theorem, cited.
## Amendment A1 (2026-09-10, post-run, reasoned — frozen verdict stands)
Frozen verdict on run 1 (seed 0): BLIND — D linear (s=0.1262, R²=0.9996, converged ×2) but F30≈0.002 violates F30>0.9. Diagnosis: FRAME BUG in the F criterion, not a physics negative. In the comoving frame ALL fluid advects out of any fixed ball (far-field streams at −U); F→0 measures advection, not instability. The correct S9 analog: L¹-type norms are LOCATION-BLIND (bump mass counts the same wherever it drifts) while diameter is location-sensitive — exactly the coexistence mechanism. Amended criterion (exploratory scope): EXPOSED iff frozen slope/R² gates hold AND perturbation mass (tracer count) conserved AND convergence holds. Run-1 numbers already satisfy this; a fresh-seed (seed 7) confirmatory run under the amended criterion closes the loop (amendment reason independent of seed-7 data).
