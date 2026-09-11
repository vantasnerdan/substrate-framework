# 0128 pre-computation freezes (drift F1/F2 — frozen BEFORE any run)

## F1 — κ misfit bar (S1 falsifier, D-08 form)
Setup: A3 filament code + Magnus term ε·ρΓ dl×(U−v) as O(ε) perturbation over
two distinct leapfrog passages (A: first overtake, B: second overtake at
different separation). Observable per passage: curvature-response field
R(s,t) = (perturbed − bare)/ε at fixed (s,t) grid. Model: R = κ·M(L) with
M from the single parameter κ (q/m slaved to κ via frozen R-EM3 relation).
- Metric: r(κ) = max over passages P∈{A,B} of ‖R_P − κM_P‖₂/‖R_P‖₂.
- FIT-PASS (charge-like): min_κ r(κ) ≤ 0.10 AND independently-fitted κ_A, κ_B
  agree within 2× (same charge across passages, not per-passage tuning).
- KILL (not-a-charge): min_κ r(κ) > 0.25 on either passage, OR κ_A/κ_B differ
  by >3×, OR best-fit κ < 0 for L > 0 (wrong-way force, mechanism named).
- Gray (0.10, 0.25]: UNRESOLVED — finer run (smaller ε, more passages), NO
  verdict either way. Tolerance 0.10 = 10× the A3 FD/entry floor (~1e-2 at
  perturbation scale, generous by design); kill margin 2.5× above pass bar.
Frozen pre-compute; the run reports (κ*, r_A, r_B) or the kill line, nothing else.

## F2 — core-singularity treatment (S3 bridge spec amendment)
Clebsch (α,β) are smooth only AWAY from vortex lines — the very carriers of
the integer charge. Bridge-spec domain rule (frozen):
- Domain: Ω_δ = R³ minus closed tubular neighborhoods (radius δ) of all
  filament cores, δ ~ a (core scale, stated per run).
- (α,β) smooth on Ω_δ; vorticity distributional (δ-supported on filaments);
  F²-bridge computed on Ω_δ PLUS separately-stated distributional core terms
  (no silent dropping: core terms printed even if zero).
- Lin constraint enforced weakly: test variations compactly supported in Ω_δ
  (vanishing near cores); core-adjacent variations EXCLUDED from the bridge
  (named, not smuggled).
- δ-dependence gate: repeat at δ/2; if the F² verdict (form/non-form) CHANGES,
  bridge = INCONCLUSIVE by construction (singular δ→0 limit, stated in
  advance — not a failure of effort, a boundary of the representation).
- If distributional core terms dominate the A-kinetic term (core carries the
  dynamics, bulk is spectator), S3's "E,B" localize to filaments — record as
  S3→S1 REDUCTION (convergence of sketches, not a kill).
