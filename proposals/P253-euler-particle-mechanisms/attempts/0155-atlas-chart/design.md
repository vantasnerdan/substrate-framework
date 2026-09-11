# 0155-atlas-chart — joint orbit chart + H⁻¹ coercivity BUILD (atlas)

Charter: shepherd 0071-chart build (shared X1+R-B object; joint cost serves
both consumers). Falsifier + stop FROZEN here pre-compute. Design only
above the freeze line; build record below it.

Missing spec (0071 review lines 199–215): joint orbit chart and relative
energy–Casimir/electric H⁻¹ coercivity for persistent charge localization.
Stationarity needs W₀·∇χ₀=0 for locked profile χ₀=F(I); persistence on the
joint orbit needs [ξ,ω₀]=0 ⟹ ξ·∇χ₀=0 (every allowed vorticity stabilizer
preserves the tag). Construction must quotient translations, fix
center/impulse and neutral rows, identify the stationary nonlinear map and
its spaces, invert the KKS/Hessian complement, and prove the translation
adjoint pairing. Known-hard wall (0071 hit it; C2 priced identically).

## Consumers (one construction, three verdicts)

- X1 D-0071: slow-collective chart at working order for the Hannay circuit
  (0140 design). Chart-without-coercivity SUFFICES for X1's sign test.
- R-B revival: chart + coercivity for persistent localization.
- R-C C2: Casimir Hessian definiteness on the member orbit (0150 design:
  indefinite → C2 dead, direct computation, decisive either way).

## Inputs (all banked, no new solves)

- A3 leapfrog: PoC-2 Newton orbit T=4.088, a=0.05; N=64/128 banked modes
  m1–m6 + section-m0 arc (0120-cipher-m2b1 receipts, unmodified).
- Banked monodromy blocks per m (0120 A3 variational code).
- 0071 review (missing-spec authority) + 0150 C2 falsifier spec.

## Construction (frozen)

1. COMOVE: quotient translations by fixing center/impulse on the banked
   PoC-2 orbit; neutral rows = translation modes (adjoint pairing checked
   against banked monodromy neutral directions, not assumed).
2. CHART: slow-collective coordinates (filament position/shape, same
   variables as the κ-fits) in the comoving frame over one period T —
   the working-order chart X1 circuits on.
3. HESSIAN: assemble Casimir second variation on the member orbit from
   banked linearization; test definiteness (C2 verdict: definite →
   C2 lives; indefinite → C2 DEAD, chart survives as kinematics only).
4. STABILIZER: test [ξ,ω₀]=0 ⟹ ξ·∇χ₀=0 on banked vorticity/tag fields
   where available; where banked data lacks the tag field, record
   TAG-ABSENT (conditional, not smuggled).
5. COERCIVITY: relative energy–Casimir H⁻¹ lower bound on the chart —
   full rigor NOT promised; report bound-or-wall honestly.

## Falsifier (frozen pre-compute)

- Hessian indefinite on the member orbit → coercivity DEAD + C2 DEAD;
  kinematic chart (steps 1–2) may still stand as conditional object.
- Adjoint pairing fails (neutral rows ≠ translation modes) → comoving
  chart ill-founded → whole build DEAD, record 0071-repeat.
- Both decided on banked numbers, zero new solves.

## Stop (frozen pre-compute)

- Banked resolution insufficient for Hessian assembly at working order →
  STOP, record 0071-repeat, do NOT smooth/approximate around it.
- Adiabatic <10× at 0124 floor → X1 circuit ill-defined (0140 stop
  stands); chart steps 1–2 may still land for R-B/R-C.
- No second rounds: one falsifier test per consumer verdict.

## Cost (frozen)

Paper + banked-number compute only. No solves, no new members, no cluster.
Drift firewall review invited on this design pre-compute; build runs on
shepherd charter regardless (ruling already given), review folds on arrival.

## Build record (below freeze line)

- 2026-09-11: design frozen, drift review invited.
