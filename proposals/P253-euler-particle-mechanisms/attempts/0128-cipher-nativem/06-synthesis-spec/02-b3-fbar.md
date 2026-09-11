# B3 F-bar freeze (pre-compute; D-08 form) — fluctuation spectrum

Construction: S4 tangle pilot extended — fixed unit-circle carrier, S=256
tangle seeds (Λ≈0.67, Rb=3); observable = |Lk| per seed ("charge-measurement
analogue" ensemble). Spectrum test: Fano factor F = Var(|Lk|)/mean(|Lk|)
(Poisson ⇒ F≈1; Skellam/difference-of-Poisson ⇒ F≈1 on counts).
- B3-HOLD (tangle component viable): F ∈ [0.5, 2.0] (Poisson-like) AND
  histogram mode at 0 with exponential-ish tail (visual check, printed
  deciles) AND seed-doubling (128→256) moves F by <20% (convergence).
- B3-KILL (tangle-vacuum dead): F < 1/3 (sub-Poissonian narrow — spec's own
  "Gaussian-narrow below CV/3" trigger) at converged seeds.
- Gray: F ∈ (1/3, 1/2) ∪ (2, 3] → UNRESOLVED + named refinement (more seeds,
  Λ-leg). Tolerance anchored to S4-measured CV (same ensemble family).
Honesty scope (pre-computed): spectrum supports the NOISE component only —
says nothing about charge quantization (missing-1/5 stand). Run reports
(F, deciles, convergence, verdict), nothing else.
