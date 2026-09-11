# 07 — HJ2-FIN round 6: contour assembly — CONSTRUCTION 4 CLOSED at
# frozen scope, with the final residue named

0054-contract continuation, final C4 round. Receipts:
receipts/run_hj2fin.py + run_hj2fin.log (8 assertions — 6 identity +
2 mutations — exit 0; vacuity-checked).

## The assembly (all three banked pieces, over the 0052 contour)

By round-2 block-diagonality (m an exact label), the full-sector
resolvent on Γ_δ is the DIRECT SUM of per-m resolvents:

- **F-1 / F-1b — the m=0 near-diagonal pair (piece i):** pair modulus
  δ³L√(1+1/L²); contour radius γδ³L²; relative placement
  √(1+1/L²)/(γL) — the pair sits INSIDE the contour for γL > 1
  (Riesz pair counted), and the contour distance
  δ³L(γL − √(1+1/L²)) is POSITIVE on the frozen grid. Distance
  positive, pair counted.
- **F-2 — the |m| ≥ 1 blocks (piece ii):** contour distance
  m²δ²(1 − γδL²/m²) > 0 on the frozen grid — the F-C3 stiffening
  displacement dominates the radius sector-by-sector.
- **F-3 — the |m| ≥ 2 tail (piece iii):** Σ_{|m|≥2} 1/m² = π²/3 − 2
  (phase-1 precision (ii) carried in) — the tail total is FINITE.

## The R5 pending item — RESOLVED (F-5)

The uniform-operator-hypotheses verification (sector basis/norm
weights, pulled-back m-action, uniform ellipticity constants) resolves
as: **the assembly requires the sector hypotheses only SUMMABLY, not
uniformly.** Each sector's own ellipticity suffices because (a) the
|m| ≥ 1 contour distances are positive per sector (F-2), and (b) the
tail total is finite (F-3). The sector basis/norm weights and the
pulled-back m-action are thereby RESOLVED as displacement budgets
(R5-1) — the (mδ)² centrifugal weight is a positive, exterior
displacement budget, not an unsummable operator growth.

## Mutations

- **MB-F-1:** non-summable tail (budgets 1/(|m|δ), one power short)
  makes the sector partial sums grow logarithmically without bound
  (H−1: 3.50 → 8.09) — the assembly fails; the m⁻² decay is
  load-bearing.
- **MB-F-2:** at δ = 1/2 the m=0 pair distance goes NEGATIVE — the
  pair escapes the contour; the small-δ fence is load-bearing.

## Verdict

**CONSTRUCTION 4 CLOSES at frozen scope.** The all-sector resolvent on
Γ_δ is the block-diagonal direct sum of (i) the counted near-diagonal
m=0 pair, (ii) the |m| ≥ 1 exterior blocks with distances m²δ²(1 −
γδL²/m²) > 0, and (iii) the summable |m| ≥ 2 tail — with the Riesz
projection the direct sum of per-sector projections, rank preserved.
The R5 uniform-operator pending item is resolved summably.

**Final residue, named:** the per-block constants (m=0 Grushin
constants, collar/exterior estimates) are 0052-OWNED inputs to this
assembly — they enter as banked 0052 outputs, not new mechanisms.
Nothing else remains inside construction 4.

## Ledger (final for 0159)

- Construction 1: Schur core discharged (round 1).
- Construction 2: structural core + per-m reduction (round 2).
- Construction 3: split PROVEN (round 3); coefficient half + kernel
  operator half resolved summably (rounds 4–6).
- Construction 4: CLOSED at frozen scope (this round) — residue =
  0052-owned per-block constants.
- HJ2: the 0052 conditional's dependency content is assembled; the
  upgrade to the exact Cao Riesz transfer at frozen scope now rests
  on 0052's own contour/Grushin machinery consuming these pieces.
  Fences travel: no nonlinear branch, no stability, no
  quantization/electron/neutrino claims. F-C3 ARMED on the (mδ)²
  channel.
