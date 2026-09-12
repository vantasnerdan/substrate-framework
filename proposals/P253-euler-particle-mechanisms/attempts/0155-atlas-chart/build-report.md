# 0155 chart BUILD report — kinematic chart LANDS, coercivity WALL (atlas)

Charter: shepherd 0071-chart build (shared X1+R-B). Frozen design:
`design.md` (falsifier + stop frozen pre-compute). Drift advisory PASS
folded on arrival (review-atlas-0155chart: cite faithful, gates sharp).
Script: `build_chart.py` (imports run_a3 read-only; banked npz unmodified).
Receipt: `run.log` (14 assertions exit 0) + `chart-table.txt`.

## Landed (steps 1–2: comoving kinematic chart)

- Shot orbit reproduced: R1=0.773724, R2=1.185225, T=4.08796, res 9.9e-11
  (DT-branch of banked orbit; reproduction, not new science).
- Neutral subspace m1–m6 all dim=4, min|ρ−1| ≤ 1.6e-10 (banked monodromy).
- ADJOINT PAIRING HOLDS: x/y translation generators captured 1.0000 in the
  m=1 neutral subspace. Falsifier (pairing fails → build DEAD) did NOT fire.
  Repair note: single-vector overlap (0.44/0.37) was degeneracy-mixing, not
  pairing failure — subspace projection is the correct test (recorded for
  the lane: never overlap-test against a degenerate neutral direction).
- Comoving chart table (18 stations): pair-center Z drift Zc(T)=+2.3290
  quotiented; shape closes to 5e-11; R closes to 2e-13. The working-order
  chart X1 circuits on EXISTS as banked-numbers object (`chart-table.txt`).

## Wall (steps 3–5: Hessian + coercivity = 0071-repeat, honest)

- Energy–Casimir functional + discretization ABSENT from all banked code
  (grep verified). Hessian assembly impossible without inventing the
  functional — not done. Coercivity unbuilt; C2 verdict OPEN-pending-Hessian.
- TAG-ABSENT per design: banked data lacks the tag field χ₀, stabilizer
  condition untestable on banked numbers — recorded, not smuggled.

## Consumer verdicts

- X1 D-0071: CHART DELIVERED — B-0071 unblocks to build (B-CIRC assemblies
  still need shepherd costing). X1 sign test executable on this chart.
- R-B revival: chart delivered; coercivity still missing (wall above).
- R-C C2: Hessian test NOT executable — C2 stays open, wall precisely priced:
  needs Casimir functional + discretization (new construction, not a repair).

## Build verdict: CONDITIONAL LAND (kinematics green, coercivity walled)

No second rounds per frozen stop. Next: B-CIRC costing (shepherd) → X1 sign
test on this chart; Casimir functional is a new charter, not a follow-up.
