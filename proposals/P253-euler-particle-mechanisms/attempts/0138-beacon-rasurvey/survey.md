# 0138-beacon-rasurvey - R-A execution paper survey (beacon)

Charter: shepherd (assumption-hunt top breaker). Paper ONLY: no
build, no assemblies — banked numbers only. Question: would R-A die
as red-team predicts, or does the prediction have a gap?

## Red-team prediction under audit

"Banked G2' measures variation 9.1/3.5/37.6 over the E-scale tube,
so min-Q ≈ 11 − 9 ≈ 2 < 5.0 HOLD: certain kill, zero cost, before
one assembly is spent." (review-beacon-0127scope, R-A break.)

## Assumption hunt (three load-bearing points)

### A1. "min-Q ≈ λ − 9.1" assumes worst-alignment — GAP (real)

Banked: max|eig(ΔQ_exact)| = 9.18 at central (trust-r3). NOT
banked: the EIGENVALUES/SIGN STRUCTURE of ΔQ, hence NOT banked:
λ_min(Q + ΔQ_central) itself. Weyl gives a LOWER bound
(λ ≥ 11.16 − 9.18 = 1.98); the actual central value could be
higher if the perturbation concentrates in the axial (35) direction
or splits the doublet favorably. "≈ 2" is the worst case
presented as the value. No banked number upper-bounds min-Q.

### A2. Revive needs the WHOLE tube, not the center — sustains expectation

Min over B(u_h, 0.43) ≤ value at EVERY tube point, including 2×
corners with |ΔQ| ~ 38. Revive (min ≥ 5 everywhere) needs
favorable alignment at MULTIPLE independent large perturbations
simultaneously — alignment miracle, probability negligible.
This is why the prediction stands as EXPECTATION despite A1.

### A3. Tube must contain u* (Q1) — cuts against revive, recorded

E = 0.43 inherits 2‖du‖; if true error exceeds it, the survey is
meaningless (protection-cost note already in 0127 amendments).
Smaller E weakens the claim it could possibly support. No
direction here favors revive either.

## Confirm-or-revive verdict: EXPECT-CONFIRM, formally UNCONFIRMED

- Expectation: R-A dies on execution (A2: 2× corners alone make
  revive miraculous; A3 removes the shrink-the-tube escape).
- Gap (A1): the central step is worst-case arithmetic, not a
  measurement — a parameter-free confirmation needs ≥ 1 executed
  survey point.
- Decisive next (ONE assembly, specified not executed): reassemble
  Q at u_h + du_central (gauged step, tested pipeline) and print
  λ_min. If ≈ 2 → prediction converts to measurement, full survey
  unnecessary (kill certain). If ≥ 5 → surprise: run the full
  amended survey (E + named coordinates) because the alignment
  structure is favorable and A2's miracle needs re-pricing.
  Cost: seconds. This is the cheapest verdict in the program and
  it is still unspent.

## Salvage (all outcomes)

- Central ≈ 2: R-A dead-by-measurement; min-Q map unnecessary;
  program record complete (prediction + confirmation, no compute
  beyond one assembly).
- Central ≥ 5: R-A revived as live; execute amended survey.
- Either way the assumption-hunt method banks: worst-case bounds
  distinguished from measurements throughout.
