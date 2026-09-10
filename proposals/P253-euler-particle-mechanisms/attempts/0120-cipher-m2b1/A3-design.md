# A3 frozen design — resonance scan under drift C1–C8 (frozen BEFORE compute)

Precedent: A1/A2 PASS-in-model disciplines (analytic-vs-measured; core honesty).
Base orbit: PoC-2 Newton T=4.088, a=0.05, two-ring leapfrog (receipts/poc2-filament/).
Every number below is a FROZEN INPUT except those labeled MEASURED.

## D1 — construction (3D filament extension)
Two rings × N nodes each (N=64 base, N=128 convergence), Rosenhead–Moore core δ=a,
Saffman self-law retained per ring. State = 6N node coords; periodic orbit =
axisymmetric PoC-2 solution revolved (m=0 seed). Variational system integrated
alongside (tangent-linear), monodromy M over one period T; per-m blocks via
discrete-Fourier projection of perturbations (m=0..6 base, m≤8 convergence leg).
Floquet ρ per m from eig(M_m); growth verdict per m, never aggregated (J1–J4).

## D2 — constraint resolutions (C1–C8)
- C1: W is MEASURED — mutual-strain rate at closest approach from the base orbit
  diagnostics (σ from passage run §3: strain 5.905/T-units scale; exact number
  printed by scan script, frozen as measured). Ruban bands = TESTS; falsifiers
  pre-attributed: strong m=1,2 growth ⇒ break #1 (core-stiffness) confirmed,
  band-edge shift quantified by core-variant leg.
- C2: BOTH conventions carried through the whole scan; the a-ladder
  a ∈ {0.10, 0.05, 0.025} gives Λ_Saff ∈ {4.38, 5.08, 5.77} and Λ_ln ∈ {2.30, 3.00, 3.69}
  — the ladder CROSSES overlap→band in both conventions, so the straddle is
  resolved empirically, not by convention choice.
- C3 (D-08): primary norm = filament-displacement L2 (node-averaged, both rings);
  control norm = kinetic-energy perturbation. Growth ⇔ |ρ| > 1+tol, tol=1e-6
  (floor: Newton residual 6.42e-11; variational integrator tolerance 1e-9;
  tol = 1e3× integrator floor, stated). Inequality-form verdict per D-08.
- C4: core model EXPLICIT — primary: Rosenhead–Moore δ=a (orbit-consistent);
  variant leg: Rankine-graft core (A2 construction) at a=0.05 only. Band-edge
  shift between legs = break-#1 quantification. No hollow-core λ=0.5 leg
  (superfluid-only; named exclusion).
- C5: m-range 0–6 + convergence leg m≤8 at N=128; per-m ρ=exp(±√μT) reported
  (μ printed, ρ derived — analytic-vs-measured labeled). Orbits: base PoC-2
  (R=1,Γ=1,a=0.05) + Buttà control (R=3,Γ=3,a=0.05: |log a|≈3 regime direction,
  radius/circulation scaled toward O(|log ε|); prediction tested THERE, not at base).
- C6 (IDEA-07): DENSE monodromy (full 6N×6N at N=64 = 384×384, direct eig) +
  NAMED soft subspace SOFT3 = {time-shift tangent Ẋ, x-translate, y-translate}
  projected out before growth verdicts; near-unit multipliers (|ρ−1|<1e-4)
  attributed to SOFT3, never to physics.
- C7: all verdicts suffixed -in-model; uncrossed gap NAMED in every verdict line:
  "live-field Euler check (filament law ⇒ Biot–Savart continuum) not crossed".
- C8: receipts/poc2-filament/ style — run_a3.py + run.log + this design frozen;
  seed/param sensitivity: N-doubling + a-ladder ARE the sensitivity legs;
  rerun-verified per 0120 discipline.

## D3 — ideas this design generates (update ⇒ ideas, primary job)
1. a-ladder-as-convention-resolver: any future Λ-straddle is settled by a ladder
   crossing the disputed boundary in ALL conventions simultaneously.
2. Buttà-control pattern: test imported-regime predictions at scaled parameters
   honoring the source regime, keep base parameters for the native claim — one
   orbit per regime, never one orbit serving two regimes.
3. SOFT3 naming pattern: every future monodromy names its soft subspace BEFORE
   compute; unnamed soft modes are where false growths hide.
4. Core-variant leg as break-quantifier: when a transfer names a stiffness break,
   the scan carries the two cores side by side and the DIFFERENCE is the result.
