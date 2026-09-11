# 02 — FBDYN D2: strain–director coupling DERIVED (bankable round 2)

Scope: 00-fbdyn-scope.md (frozen pre-compute, 92825c5e), round D2:
enumerate the symmetry-allowed (ε, n̂, ∇n̂) scalars and derive the
leading coupling. Discipline note (D1 repairs applied): this round
claims the COUPLING structure and its coefficients — propagation
consequences (ω² shifts, direction-dependence) are D3's object, and
no spectrum/PSD statement is made here. Receipts:
receipts/run_fbd2.py + run_fbd2.log (8 assertions — 5 identity +
3 mutations — exit 0).

## Derived structure (each step receipt-backed)

- **RD2-1 (COUNTING, structural):** NO bulk coupling exists at
  O(ε)(∇n̂): one ε block + one gradient block + one n̂ = 5 indices —
  ODD, zero perfect delta-pairings (combinatorial receipt); the
  6-index even candidates contract to εᵢⱼn̂ᵢ∂ⱼ(n̂·n̂)/2, IDENTICALLY
  zero on the unit sphere (n·n = 1 pointwise). The coupling starts
  at O(ε)(∇n̂)². The COUNT is a falsifiable model signature.
- **RD2-2 (angular 4-tensor):** ∫dΩ ŝᵢŝⱼŝₖŝₗ =
  (4π/15)(δᵢⱼδₖₗ + δᵢₖδⱼₗ + δᵢₗδⱼₖ) — 7 representative components by
  direct integration.
- **RD2-3 (LEADING COUPLING):** via the first-order affine
  separation map r → (I+ε)r with measure dilation (1+tr ε),
  assembled from the exact angular contractions and checked against
  the closed formula on all 9 index pairs:
  **W_coup = (K p²ξ²M₄ 4π/15)·[4εᵢⱼ + 7(tr ε)δᵢⱼ]·∂ᵢn̂·∂ⱼn̂.**
  The strain modulates the director stiffness — the Vikulin J(ε)
  analog on the medium side (F-C's derived J(ε structure).
  Allowed-structure COUNT = 2 (deviatoric + trace), receipted.
  Coefficients derived from the declared kernel; number model-level
  (same tier as K_n, stated).
- **RD2-4 (objectivity):** W_coup invariant under the JOINT rotation
  (ε, ∇n̂) → (RεRᵀ, R∇n̂Rᵀ) — exact on symbolic data for a concrete
  rotation (RB9 analog extended to ∇n̂).
- **RD2-5 (reduction/continuity):** W_coup = 0 identically at zero
  gradient — the static F-B-lite tier is recovered untouched; the
  localized pre-stress K·p·(n̂·εn̂ − tr ε/3) at uniform ẑ reproduces
  RB6's built linear term exactly.

## Mutations (negative controls — the kills are live)

- **MB-D2-1:** rotating ε WITHOUT ∇n̂ changes W_coup (nonzero on
  symbolic data) — the joint rule is load-bearing; a build asserting
  rotated-eps-only invariance is caught.
- **MB-D2-2:** an ε-coupled twist insertion flips sign exactly under
  inversion (carrying D1's parity-odd receipt) — the frozen
  achirality declaration excludes it, detected.
- **MB-D2-3:** the forbidden single-gradient candidate reduces to
  the constraint zero EXACTLY — a bulk O(ε)(∇n̂) coupling claim is
  caught as identically zero, not a fitted small value.

## D2 verdicts (PROPOSED — drift review requested)

- **FB-D-coupling:** the leading strain–director coupling is DERIVED:
  count = 2 allowed structures at O(ε)(∇n̂)², both coefficients
  given by the declared-kernel formula above; objectivity is joint
  (RD2-4); continuity with F-B-lite is exact (RD2-5). NO kill fired;
  no STOP rule triggered (the enumeration used only the declared
  direction-statistics model — no positional co-variation needed).
- Mechanism hygiene (MA-4 pattern): the coupling enters as stiffness
  MODULATION on the medium side; no force template is used or
  needed — the dead force-template family plays no role here.
- Fences unchanged: model-level family build; no carrier/electron
  claim; p fixed; achiral declared; no measurement; no LANE-1
  membership. Dated deltas: none this round (the derivation ran
  inside the frozen scope with the D1 deltas (i)–(iii carrying over).

## Next

D3 dispersion round: quadratic expansion about uniform (n̂₀, rest) —
director-wave dispersion ω²(k) = (K_n k² + pre-stress contribution +
coupling k-terms)/I_n, REAL second-order required on p ∈ [0, 1);
direction-dependence formula; the coupling's k-dependent readout
δ(ω²)(k) (medium-side readout, F-C RC5 analog). Kinetic structure:
I_n = ρξ³κ_n scaling-level (κ_n priced-not-derived — footnote
travels).
