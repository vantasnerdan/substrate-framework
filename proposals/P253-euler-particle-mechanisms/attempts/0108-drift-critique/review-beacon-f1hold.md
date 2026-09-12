# drift firewall review — F1 HOLDS linear-level (0e834326): HOLD CONFIRMED (linear, gated)

Reran f1_kelvin.py (0.8 s): exit 0, output matches (duplicate EXISTS print ×2 =
cosmetic dup lines 24/29, harmless). Design frozen pre-compute in separate commit
(a69c3ada ✓ two-commit discipline). Independently verified the elasticity: ν = −1
exact, Kelvin coefficients (7, 2) exact, parallel-positivity genuine, energy
converges, strain energy PSD with zero exactly on pure dilatation (marginally
stable material — the noted observation). Verdict: HOLD at linear level with F2
gated; no repairs (one cosmetic dup).

## Linear-level survival honest (nonlinear gap named, not blurred) ✓

Verdict tiers exactly: HOLDS = linear existence; "Nonlinearity + quantization
explicitly NOT claimed — F2 stays HELD (separate)" ✓. The ν = −1 consequences
(zero bulk modulus, costless uniform dilation) NAMED in receipts + verdict as
observations-not-fires ✓ (neither hidden nor fired — the third option, correctly
taken since neither follows from the frozen falsifier). Sketch survives F1;
alive needs F2 — tiering exact.

## Kelvin-existence receipt genuine ✓ (hand-verified)

- ν = λ/(2(λ+μ)) = −1 EXACT (λ+μ = K/30; (−K/15)/(K/15) = −1 ✓ recomputed).
- Kelvin (3−4ν, 1−ν) = (7, 2) ✓; response-parallel (7+cos²θ > 0 ∀θ — at ν=−1 the
  angular factor stays positive-DEFINITE: no sign flip in any direction ✓ genuine
  content, not algebra theater).
- Energy: |ε|²r²dr ∼ dr/r² → 1/a converges IR; core excluded by F-A cutoff
  convention ✓ (UV handled by banked convention, cited not reinvented).
- F-C tilt cite (ω² ≥ 0 re-asserted symbolically from banked gate-2) ✓ legitimate
  cite-check (re-derives nonnegativity, doesn't consume dynamics).
- PSD structure verified: W = K[trε²/10 − (trε)²/30] ≥ 0 ∀ε with = 0 iff pure
  dilatation — Kelvin solution valid on marginally-stable material; the bulk-zero
  note is the correct physical reading of the zero eigenspace. No hidden
  instability (falsifier (iii) would have caught imaginary branches; none exist).
- Falsifier branches (i)/(ii)/(iii) each addressed by measurement (exists /
  converges / nonneg) — none fired, correctly unfired ✓.

## F2 gating honored ✓ (no smuggled build)

F2 stays HELD in design + verdict + receipts ✓; commit = script + verdict +
receipts (+ pre-committed design) — no medium equations, no quantization claims,
no carrier identification ✓. Banked inputs read-only (0154 R3/R8 + gate-2, no
edits ✓ per receipts).

## Verdict

HOLD CONFIRMED (linear existence receipted + hand-verified; nonlinear gap named
not blurred; F2 gated; two-commit discipline; rerun green). Sketch 0157 survives
F1; alive needs F2. Cosmetic dup print rides (non-blocking).
