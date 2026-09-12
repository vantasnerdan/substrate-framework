# S4 García–Hassainia–Hmidi verified (beacon 0110)

Source: *Time-periodic leapfrogging vortex rings in the 3D Euler equations*,
arXiv:2603.21644v1 [math.AP] 23 Mar 2026. Fetched abs + full intro/Theorem
1.1 + §1.2–1.3 this session (V3; 6181-line extraction, artifact://40).

## Exact consumed scope

- Theorem 1.1: for fixed `κ > 0`, `ε ∈ (0,ε_0)`, and `λ` in a Borel set
  `C_ε ⊂ (a,b)` of asymptotically FULL measure, a GLOBAL axisymmetric
  no-swirl Euler solution of patch type `q = ε^−2 1_{D_1(t)} + ε^−2 1_{D_2(t)}`
  exists; domains `D_j(t)` are T-periodic modulations around the filament pair
  `(P_1,P_2)` in an almost-uniformly translating frame (drift
  `U_ε ∼ √(ln|ln ε|)/(4√(2κ))`, singular in `|log ε|`).
- Method: contour-dynamics desingularization of two filaments → Hamiltonian
  patch system; degenerate KAM + pseudo-differential reducibility (4
  conjugations) + Nash–Moser; singular small divisors degenerate with
  thickness ε; speed modulation kills the first sine mode only.
- Hypotheses: axisymmetric NO-SWIRL; nearly concentric pair, core area O(ε²),
  separation O(|log ε|^−1/2), ring radius O(1) — the S3 regime, explicitly
  adopted (§1.2: "will be considered in the present work").
- Shape: elliptical (not circular) cores from the `(ϱ = r²/2, z)` Hamiltonian
  variables; non-rigid periodicity `f(τ,θ) = (1/8)(2p_11)^−3/4 cos3θ + …`.

## Perturbation class / norm

KAM persistence (Cantor-set λ, Sobolev index s large for boundary profile
`f ∈ H^s(T)`), not a stability neighborhood of a named steady carrier. No
claim about general nonlinear stability — periodic existence only.

## S2-carrier transfer verdict

- Same compatible regime as S3 (O(1) radius, thin core, no-swirl); frozen-issue
  row satisfied: parameters/regularity reviewed, periodic existence does NOT
  supply general nonlinear stability (source claims none).
- What transfers NOW: nothing as a theorem (S2 member ≠ constructed pair;
  λ-Cantor qualification has no S2 analogue yet). Supplies the P3 target
  shape: an all-time two-ring motion exists in the same regime the S2 carrier
  lives in — motivation, not a lemma.
- Relation to S3: S4 is the long-time counterpart (S3 = finite window,
  S3's periodic case "very difficult" → S4 solves it via KAM in translating
  frame). Both live in the R1 interaction column, both pre-stability.
