# 04 — FBDYN D4: PSD window + verdicts (bankable round 4, lane-complete)

Scope: 00-fbdyn-scope.md (frozen pre-compute, 92825c5e), round D4:
total 2nd-variation PSD at long wavelength for p ∈ [0, 1) —
coupling-induced negative modes assessed as kill-or-named-window-
shrink; verdicts per the 01 F-B column paths extended with FB-D
items. Works off the D3 dispersion objects (ω²_± = [K_n k² +
2C_c k²(k̂ᵀMk̂) + K_p λ_±(T)]/I_n; receipts/run_fbd3 banked).
Receipts: receipts/run_fbd4.py + run_fbd4.log (7 assertions —
5 identity + 2 mutations — exit 0). Discipline carried: no
measurement; speeds carry κ_n; model-level numbers labeled; the
flexo-analog acknowledgment (drift R1 540d2c3b) travels.

## Derived structure (each step receipt-backed)

- **RD4-1 (STRUCTURAL MARGIN):** K_n = 10·C_c EXACTLY — the ratio
  (8π/3)/(8/15) = 5 is M₄-INDEPENDENT (one declared kernel feeds
  both constants) — so the k² coefficient is 2C_c(k̂ᵀMk̂ + 5):
  PSD condition (a): **k̂ᵀMk̂ ≥ −5**, a STRUCTURAL margin, not a
  fitted constant.
- **RD4-2 (WINDOW EQUIVALENCE):** ω²_± is AFFINE in k² (d/d(k²)
  constant) with value K_p λ_±/I_n at k → 0 — therefore
  **PSD for all k ≥ 0 ⟺ λ_±(T) ≥ 0 AND k̂ᵀMk̂ ≥ −5**: (a) controls
  the large-k side, the gap the k → 0 side. Exact, both sides
  receipted.
- **RD4-3 (BUCKLING-ANALOG SOFT MODE — the window's mechanism):**
  uniaxial ε = diag(ε⊥, ε⊥, ε∥), k̂ = x̂. TENSION side (ε∥ = t > 0):
  λ_− = −t and ω²_− has the EXACT root k_c² = K_p t/(2C_c(7t + 5)),
  negative inside (ω²_−(k_c/2) = −3K_p t/(4I_n)), positive beyond
  (ω²_−(√2 k_c) = +K_p t/I_n) — sign flip receipted. COMPRESSION
  side (ε∥ = −t): λ_− = 0, ω²_− = 2C_c(5−7t)k²/I_n ≥ 0 given (a) —
  no soft mode. Mechanism NAMED: axial tension beyond isotropy
  destabilizes transverse tilt at long wavelength (Euler-buckling
  analog on the medium side).
- **RD4-4 (COMPRESSION BOUND from (a)):** isotropic ε = −cI gives
  k̂ᵀMk̂ = −25c exactly → (a) ⟺ c ≤ 1/5: the k² coefficient
  survives 20% isotropic compression — outside any linear-tier
  strain window; boundary and both sides receipted on rationals.
  RB10/RB10b acoustoelastic O(p) caveat travels with this bound.
- **RD4-5 (STATIC CONSISTENCY + p-chain):** at ε = 0: ω² =
  K_n k²/I_n ≥ 0 on p ∈ (0, 1), equality ONLY at p = 0 (MB-D1-2
  echo); W_coup = 0 at zero gradient (D2 RD2-5) — the dynamical
  branch adds NO static constraint; **P3-B-dyn** (scope §5's
  registered per-family item) lands CONSISTENT with the static
  branch.

## Mutations (negative controls — the kills are live)

- **MB-D4-1:** a wrong-C_c build (M₄ → M₄/2) breaks K_n = 10C_c —
  a mistyped or fitted coupling constant moves the margin off the
  receipted 5 and is detected; the margin is not adjustable.
- **MB-D4-2:** a windowless all-ε PSD claim is caught by the
  RD4-3 tension counterexample — ω²_−(k_c/2) = −3K_p t/(4I_n) < 0
  EXACTLY: an exact negative mode outside the named window, inside
  any windowless claim. The PSD verdict is a WINDOW verdict.

## D4 verdicts (PROPOSED — drift review requested on the full FB-D set)

- **FB-D window (the kill-(i)/window-shrink question):** NO KILL
  FIRED. The admissible window is NAMED: **{p ∈ (0, 1)} × {ε :
  λ_−(T) ≥ 0 and k̂ᵀMk̂ ≥ −5 (all unit k̂)}** — inside it the total
  2nd variation is PSD at all k ≥ 0. Outside the window the model
  PREDICTS the buckling-analog soft mode with exact k_c formula —
  a named window-shrink WITH mechanism and a falsifiable
  instability prediction, not a model failure and not a hidden
  negative. Kill-(i) (negative mode within the declared window)
  does not fire; the ε∥ > ε⊥ region is where the falsifiable
  prediction lives, and any "PSD everywhere" claim is receipted
  false (MB-D4-2).
- **What-licenses-what (FB-D column):** D1 licenses the stiffness
  (K_n > 0, p²-scaling, k¹ exclusion); D2 licenses the coupling
  (count-2, joint objectivity, exact reduction); D3 licenses the
  dispersion FORMULA (real second-order, direction-dependence,
  polarization-degenerate readout) — kill-(iii)'s kinetics half is
  supplied, closure awaits FIREABILITY (independent-p route,
  named-open, owner-level); D4 licenses the PSD WINDOW statement
  above at model level. No lane-level conjunction fires; no
  carrier/electron claim; speeds carry κ_n; model-level numbers
  ride the declared profile, labeled per round.
- **P3-B-dyn (per-family item, lands here):** the dynamical-branch
  spectrum = the D3 formula + the RD4-2 window, family-wide within
  the declared model class; falsifier FB-D-waves with the frozen
  hygiene (INCONCLUSIVE within priced orders) + fireability
  status; consistent with the static branch (RD4-5).
- **Fences restated:** model-level family build; no carrier/
  electron claim (the RC5-analog readout is a coupling-constant
  handle on the medium side); p fixed (p-field dynamics out of
  scope, named limit); achiral declared; no measurement; no LANE-1
  membership. Scope deltas: none this round (all four rounds ran
  inside the frozen scope; D1's dated deltas remain the only ones).

## Next

The FB-D lane is COMPLETE at frozen scope (D1–D4 all banked, each
with receipts + mutations; D1/D2 drift-reviewed CONDITIONAL PASS
with repairs paid; D3 CONDITIONAL PASS with R1 paid this commit;
D4 drift review now requested on the full claim set). Remaining
named-opens, owner-level: (1) FB-D-waves FIREABILITY route
(independent-p measurement protocol) — the falsifier stays
UNFIREABLE-named until then; (2) the F-B build decision (SYN
D-bucket D4 gate, receipt-round dues per 06 §6) — the dynamical
director completes the priced new-mechanism object that decision
prices.
