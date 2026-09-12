# 02 — 0062-R2: resonance localization DONE (bankable round 2)

Scope: 00-scope.md (frozen, this attempt), round R2 per the frozen
plan, authorized (#110 + GOAL RESTATE routing): first source-bearing
resonance of the exact block (18) at fixed toroidal n — exact
crossing at the BLOCK level, not the raw factor (17). Consumed at
stated scopes: 0159 per-sector contour/Grushin pattern (m = 0 chain
near-diagonal pair inside the contour; |m| >= 1 uniformly exterior
via the F-C3 stiffening displacement); 0052 scaling as localization
target only; 0058 nondivisibility counterexample as the standing
raw-factor mutation anchor; the R1 seed (0161/01) as the coupling's
origin. Receipts: receipts/run_0062r2.py + run_0062r2.log (6
assertions — 4 identity + 2 mutations — exit 0). Discipline: the
receipts validate the encoded model-level (two-mode) identities
ONLY; the continuum sandwiched trace + distorted adjoint at
lambda_* is R3's object and is NOT claimed here.

## Localized structure (each step receipt-backed)

- **RR2-1 (sector displacement budget, EXACT):** the block diagonal
  per poloidal sector m carries the F-C3 stiffening displacement
  +m²δ² (sign receipted at 0159). The |m| ≥ 1 crossing would need
  the displacement inside the controlled band (contour radius
  γδ³L²): displacement/radius = m²/(γδL²) → ∞, and the crossing
  scale m* = √(γδ)·L < 1 for ALL δ ∈ (0, 1) at γ = L = 1 — NO
  integer m ≥ 1 sits inside the controlled window: **the |m| ≥ 1
  poloidal sectors carry NO source-bearing resonance there** (the
  0159 exterior verdict, consumed at the block; the localization
  reduces to the m = 0 chain).
- **RR2-2 (m = 0 two-mode Grushin crossing, EXACT):** the
  near-diagonal pair M(λ) = [[D₁−λ, g], [g, D₂−λ]] has eigenvalues
  **λ_± = (D₁+D₂)/2 ± √(((D₁−D₂)/2)² + g²)** — annihilation of
  det M IDENTICAL, splitting (λ₊−λ₋)² = (D₁−D₂)² + 4g² a SUM OF
  SQUARES: the block-level crossing is REAL for real data (two real
  resonances, exact formulas).
- **RR2-3 (source-bearing coupling, EXACT):** g is the Hodge
  element from the R1 seed's (16) structure:
  g = (ω₀·k)/|k|² · ⟨e₁, k×e₂⟩ on transverse polarizations;
  witness ⟨e₁, k×e₂⟩ = −1 ≠ 0 (e₁ = x̂, e₂ = ŷ, k = ẑ): **g ≠ 0
  exactly** — the resonance is SOURCE-BEARING (a bare transport
  divisor gives g = 0 exactly).
- **RR2-4 (transparency limit):** g → 0 recovers the bare transport
  crossings {D₁, D₂} exactly (ordered witness) — the block formula
  is the continuous extension of the raw picture; the README's
  transparency fork degenerates to the raw crossings at this trace
  when the coupling vanishes by exact symmetry.

## Mutations (negative controls — the kills are live)

- **MB2-1 (raw-factor substitution; 0058 anchor):** the raw
  crossing λ_raw = D₁ misses the coupling shift:
  λ₊ − D₁ = √(((D₁−D₂)/2)² + g²) − (D₁−D₂)/2 ≠ 0 for g ≠ 0
  (concrete witness D₁ = D₂: shift = |g|) — replacing the block
  crossing by the raw factor D(I) = 0 is caught (the
  nondivisibility class; the README's "spectrum of (18), not merely
  the raw interval of (17)" requirement is enforced by the receipt).
- **MB2-2 (stiffening-sign flip):** flipping the displacement sign
  moves every |m| ≥ 1 localized level by −2m²δ² ≠ 0 — the
  |m| ≥ 1 exclusion is SIGN-DEPENDENT: the F-C3 sign budget is
  load-bearing for the localization claim; a softening sign would
  move sector candidates back toward the controlled window (a
  different, unreceipted build).

## R2 verdicts (PROPOSED — drift review requested)

- **R2-resonance: BLOCK-LEVEL CROSSING EXISTS.** At fixed toroidal
  n, the localization is complete at two-mode model scope: the m = 0
  chain carries a REAL source-bearing resonance pair with exact
  formulas; the |m| ≥ 1 sectors are contour-killed by the receipted
  displacement budget; the coupling is nonzero by the R1 seed's
  Hodge structure. **Kill-(i) NOT fired** (a resonance exists at
  declared order at block scope); kills (ii)/(iii) not yet reached
  (V_* and the branch are R4/R5 objects).
- **What this round does NOT claim:** the continuum
  limiting-absorption sandwiched trace and distorted adjoint at
  λ_* (R3's object — the exact Green/Grushin trace with finite
  circulation rows); the continuum eigenvalue theorem (the receipts
  validate the two-mode encoded algebra; the function-space
  statements are governed by the 0062 README and the consumed
  0159/0052 scopes). The resonance is LOCALIZED AT MODEL SCOPE; its
  physical (continuum) existence rides R3's trace construction.
- Scope deltas: none (round ran inside the frozen scope; single
  fixed-n block; no all-sector claim — RR2-1 is per-sector at the
  block, consistent with STOP-1).
- Fences restated: P2/LP2 classical only; no carrier/electron/
  quantum claim; no measurement; append-only; drift firewall
  reviews every round.

## Next

R3 TRACE per the frozen plan: the source-specific sandwiched
core/interface/exterior Hodge Green/Grushin trace at the localized
λ_* — distorted adjoint ψ_*^* via (20), residue normalization (10),
finite circulation/interface/exterior rows, range condition (21) in
distorted form; raw-N(I_*) mutation live. Bankable: the trace
functional exists at model scope. Then R4 (V_* evaluation) and R5
(branch decision).
