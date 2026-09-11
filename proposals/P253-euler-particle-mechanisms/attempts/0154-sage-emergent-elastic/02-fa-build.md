# 02 — F-A build: emergent shear modulus of the random frozen tangle

Charter: shepherd 0154 F-A (0153-adjacent shear work is a different
surface: imposed background ON the carrier; this derives the medium FROM
the family). Verdict vocabulary per 01-verdicts-per-family.md. Receipts:
receipts/run_fa.py + run.log (14 assertions, self-counted, exit 0).

## Model assumptions (stated, all model-level)

1. F-A family: statistically isotropic ensemble of closed vortex
   filaments, single circulation quantum Γ, line-length density L₀
   (length per volume), core radius a, outer cutoff ℓ (inter-filament
   spacing). This is the S4 vacuum family taken as exact Euler data
   (frozen-in: 0147 I-Sing is the exactness license — d/dt Lk = d/dt
   Slk = 0; reconnection outside the representation).
2. Energy density (standard cut-off form): ε₀ = K with
   K = (Γ²/4π)·L₀·ln(ℓ/a); mutual terms absorbed into the O(1) constant
   of the log — STATED, not derived (cutoff convention).
3. In-window response is affine: for deformation rates ω ≫ 1/τ (τ =
   tangle relaxation time) the frozen relative shape cannot reorganize,
   so the response to a coarse volume-preserving F is the affine
   push-forward of every material line element.

## Derivation

Material line elements stretch by |Fn|; with F = RU (polar), |Fn| = |Un|
— objectivity automatic; U = (I + 2ε)^{1/2} with ε = sym(displacement
gradient) to O(ε). Expansion:

  |Fn| = 1 + n·εn + ½[n·ε²n − (n·εn)²] + O(ε³).

Isotropic angular averages (receipt R1a–c): ⟨n·εn⟩ = tr ε/3,
⟨n·ε²n⟩ = tr(ε²)/3, ⟨(n·εn)²⟩ = [(tr ε)² + 2 tr(ε²)]/15. Hence

  ⟨|Fn|⟩ = 1 + (tr ε)/3 + tr(ε²)/10 − (tr ε)²/30 + O(ε³),        (R2)

and the energy density ε(E) = K·⟨|Fn|⟩ gives the quadratic form

  W(ε) = K·[ tr(ε²)/10 − (tr ε)²/30 ].

Matching W = μ·tr(ε²) + (λ/2)(tr ε)² (receipt R3):

  **μ_aff = K/10 = Γ² L₀ ln(ℓ/a) / (40π) > 0,   λ = −K/15.**

λ < 0 is the incompressible artifact (volume locked; the (tr ε)² term is
pressure-like); the isochoric response is governed by μ alone. Pure-shear
sanity R7: W = 2μ ε₁₂². Objectivity receipts: R4 (rotation spot check),
R5 (certified form is ω-free; the naive gradient form carries skew
dependence — mutation MA-1 detects the error); average audit MA-2
(corrupted fourth moment changes the form — detected).

## FB verdicts (per 01; drift rules)

- FB-1 (F-A): formula half ESTABLISHED (μ_aff > 0, receipt R3/R8);
  window half PRICED — τ > 0 argued by Kelvin-transit scaling
  (τ ~ ℓ/c_Kelvin; smoothing shortens Λ at fixed topology, so topology
  pins Lk, NOT length — the honest limit of the pinning mechanism), not
  derived. FB-1 verdict PENDING window; no static μ(0) claim.
- FB-2 objectivity: ESTABLISHED (R4, R5).
- FB-3 bookkeeping: ESTABLISHED at the stated order — expansion exact
  through O(ε²) with explicit O(ε³) remainder; log-argument drift
  ℓ(ε)/ℓ absorbed as an O(ε²) constant into K (stated convention, no
  unaccounted sink).
- FB-4 inheritance (gate 1): trivially expected for F-A — ⟨H⟩_bg = 0 as
  an ensemble property of the isotropic family (cited model property,
  the S4 pilot's own statistics); elastic waves of the medium carry no
  helicity flux; conservation inherited trivially. Re-verdict clause
  stands (01): if the S4-measured background helicity level proves
  substantive, FB-4 re-opens as substantive.
- FB-5 propagation (gate 2): in-window shear waves ω = c_s k,
  c_s = √(μ_aff/ρ), real frequencies, second order (R8). ESTABLISHED
  in-window, conditional on the same window as FB-1.
- FB-6 scale separation: HOLD-conditional (ℓ ≫ a pricing needs carrier
  numbers; parameter ratio stated, not derived).

## What F-A contributes to LANE-1

The F-A half of LANE-1 (μ > 0 in-window) is now receipt-backed at
paper level. LANE-1 fires only with the CONTRAST half (03) — neither
half alone licenses anything (01).
