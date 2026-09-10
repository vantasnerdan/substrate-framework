# S7 Faddeev–Niemi verified (beacon 0110)

Sources: *Toroidal configurations as stable solitons*, hep-th/9705176v1,
22 May 1997 (detailed companion to the Nature 387 (1997) 58–61 letter cited
in the frozen issue, which is its ref [7]). Full PDF fetched this session
(V5; 392-line extraction, artifact://38). The Nature letter itself was not
separately fetched — boundary recorded below.

## Exact consumed scope

- Action (1): `S = ∫[(∂n)²/2e² + (n·∂n×∂n)²/4g²]` — O(3) sigma model + FOURTH-
  order Skyrme term. Stabilization is Derrick-virial: `H = E₂+E₄ →
  ρE₂ + E₄/ρ`, virial `E₂ = E₄` (5). Without E₄ no 3D finite-energy solitons.
- Topology: `n: R³∪{∞} ≅ S³ → S²`, `π₃(S²) ≅ Z`, Hopf invariant (8)/(13);
  energy bound `H ≥ c·|Q_H|^3/4` (9) — existence SUGGESTED, not proved.
- Existence evidence: NUMERICAL gradient-flow (29)–(34) with virial
  renormalization `G(τ) → G∗`, PDE2D FEM, ~200 CPU-hours/phase, Q_H =
  0.99997, tube-shaped energy density; small-r convergence difficulties
  honestly reported (§5). Analytic solution "appears impossible" (§4).
- Scale: sole dimensionful coupling G sets size/shape; solitons at other G
  by scale transform. Length of each knot is NOT predicted from topology —
  the paper's own critique of phenomenological string models (§1) applies.

## S2-carrier transfer verdict

- Topology ≠ particle (frozen-issue warning confirmed at source): Hopf charge
  is an integer of the n-field; no Euler derivation, no quantization, no
  species map. Which terms/constraints/dynamics follow from Euler: NONE
  shown — E₄ has no Euler origin; Euler has no fixed scale G.
- Boundary: Nature-letter-specific claims beyond this companion (trefoil
  evidence) not re-verified; irrelevant to transfer (would not supply Euler
  origin either).
- R2 design constraint preserved: any Euler–FS bridge must derive a quartic
  stabilizer + scale from the substrate, or import them as NAMED hypotheses
  per issue203-frozen §6 (with foundational-proposal cost).
