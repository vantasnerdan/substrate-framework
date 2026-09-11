# 02 — F-A build: emergent shear modulus of the random frozen tangle

Charter: shepherd 0154 F-A (0153-adjacent shear work is a different
surface: imposed background ON the carrier; this derives the medium FROM
the family). Verdict vocabulary per 01-verdicts-per-family.md. Receipts:
receipts/run_fa.py + run.log (19 assertions, self-counted, exit 0; F-A
dues paid — see FB-1b, FB-3, CONTRAST below).

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
— objectivity automatic. STRAIN IDENTIFICATION (Biot, un-glossed per
drift dues): the expansion is exact in the rotation-free Biot strain
ε_B = U − I; to the certified order ε_B = sym(F − I) + O(ε²)-rotation
corrections (a |W|²_F/6-class piece), which enter |Fn| only additively —
receipts RB/RB2 verify the rotation piece separates with no deviatoric
strain dependence, so the deviatoric coefficients below are unaffected
and the form works purely in symmetric-ε space (R4 conjugation). With
that identification, expansion:

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
  through O(ε²) with explicit O(ε³) remainder; SOFTENED per drift dues:
  the O(ε²) log-argument drift is NOT a silent constant — it is a PRICED
  cutoff-slaving sensitivity (next bullet). No unaccounted sink beyond
  that priced channel.
- FB-1b cutoff-slaving sensitivity (PRICED, drift dues): if the outer
  cutoff slaves to line density (ℓ ~ L₀^{−1/2}), K drifts at O(ε²) WITH
  the stretch and feeds the quadratic form at O(1)-relative: receipt RS
  gives μ_eff = (K/10)(1 − 1/(2·ln₀)) — same order, same sign, number
  moves by an O(1) factor (~2× for ln₀ ≈ 1). μ's SIGN and FORM are
  robust under either slaving; μ's NUMBER carries O(1) outer-cutoff
  sensitivity — stated as the priced caveat it is.
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
half alone licenses anything (01). CONTRAST receipt REPAIRED per drift
dues (R6b was vacuous-as-coded): the wave-gas zero is now DERIVED —
W_wave = f(det F) with the EXACT incompressibility constraint det F ≡ 1
freezing the state-function argument (R6b-1: the only O(ε²) route is the
det-deficit, worth −(f′/2)·tr(ε²); R6b-2: the constraint annihilates it;
R6b-M: relaxing the constraint revives the response — the assumption is
load-bearing and labeled); R6c-v2 compares the two DERIVED forms
(difference = the memory term on isochoric strains — real discrimination,
not internal consistency).
