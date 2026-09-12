# 03 — FBDYN D3: director-wave dispersion DERIVED (bankable round 3)

Scope: 00-fbdyn-scope.md (frozen pre-compute, 92825c5e), round D3:
quadratic expansion about uniform (n̂₀ = ẑ, rest) in the PRINCIPAL
frame (ε₁₃ = ε₂₃ = 0 — the frame where the pre-stress linear term
vanishes identically; generic-frame linear term 2Kp(u ε₁₃ + v ε₂₃)
receipted as the frame's motivation). Frame-covariance of the
formula is carried by the RD2-4 joint-rotation receipt; principal-
axis diagonalization is a kinematic frame choice, not new
statistics — STOP-1 clear. Receipts: receipts/run_fbd3.py +
run_fbd3.log (8 assertions — 5 identity + 3 mutations — exit 0).
Discipline (D1/D2 repairs carried): this round claims the ω²(k)
FORMULA, its realness/second-order structure, and the
direction-dependence; the PSD window verdict on p ∈ [0, 1) is D4's
object. No spectrum/PSD measurement is made or claimed. The
flexo-analog acknowledgment (drift R1 540d2c3b) TRAVELS: at
NON-uniform ε the divergence-silent O(ε)(∇n̂) family activates
O((∇ε)n̂) bulk pieces + anchoring-like boundary terms, outside this
uniform-ε formula.

## Derived structure (each step receipt-backed)

- **RD3-1 (constraint/projection, EXACT):** on the unit
  parametrization n̂ = (tu, tv, √(1−t²(u²+v²))): the pre-stress term
  is EXACTLY quadratic — W_pre − W_pre(0) = Kp·t²[(ε₁₁−ε₃₃)u² +
  2ε₁₂uv + (ε₂₂−ε₃₃)v²] with NO remainder at any order (n̂₃² = 1 −
  t²(u²+v²) exactly); the coupling gradient block reduces as
  ∂n̂·∂n̂ = t²(∂u∂u + ∂v∂v) + O(t⁴) with the odd-power coefficient
  EXACTLY zero (n̂ → −n̂ evenness, structural).
- **RD3-2 (realness + second-order kinetics):** Euler–Lagrange from
  the quadratic Lagrangian gives In ω²a = Sa with
  **S(k̂) = [K_n + 2C_c(k̂ᵀMk̂)]k²·I₂ + K_p·T** — REAL SYMMETRIC on
  symbolic data; the eigenvalue discriminant (S₁₁−S₂₂)² + 4S₁₂² is
  a sum of squares of reals → ω² REAL structurally; In = ρξ³κ_n > 0
  → genuine second-order oscillator (κ_n priced footnote travels).
  The E-L residual against the closed form is EXACTLY divisible by
  (1 − |k̂|²): agreement on the unit-k̂ cone (only unit directions
  are kinematic).
- **RD3-3 (dispersion formula, two routes):**
  **ω²_±(k, k̂) = [K_n k² + 2C_c k²(k̂ᵀMk̂) + K_p λ_±(T)]/I_n**
  with M = 4ε + 7(tr ε)I, T = [[ε₁₁−ε₃₃, ε₁₂],[ε₁₂, ε₂₂−ε₃₃]],
  λ_± = (trT ± √discT)/2, discT = (T₁₁−T₂₂)² + 4T₁₂²; K_n = 64πKp²ξ²
  (D1), C_c = (4π/15)Kp²ξ²M₄ (D2, M₄ = 24). Annihilates
  det(In ω² I − S) symbolically and EQUALS the direct eigenvalues of
  S on concrete anisotropic data at k̂ = x̂, ẑ, (x̂+ẑ)/√2. Direction
  dependence NONZERO between x̂ and ẑ — the falsifiable anisotropy
  lives (uniaxial ε: the k² coefficient varies as
  4(ε⊥ + (ε∥−ε⊥)cos²θ) + 7(2ε⊥+ε∥); the polarization split is
  Kp(ε⊥−ε∥) via λ_±).
- **RD3-4 (coupling readout, polarization-degenerate):**
  S − S|_{C_c=0} = 2C_c k²(k̂ᵀMk̂)·I₂ EXACTLY — the coupling shifts
  BOTH polarizations identically (k̂-dependent only). The F-C RC5
  analog: a measurable handle on the coupling constant, medium side
  (no carrier claim; MA-4 pattern).
- **RD3-5 (p-chain + inertia footnote):** ω²_±(p = 0) = 0
  IDENTICALLY — the unpolarized medium carries no director waves
  (dispersion-level echo of MB-D1-2). Every speed carries
  κ_n^(−1/2) via In = ρξ³κ_n (priced, not derived; scope 2.5).

## Mutations (negative controls — the kills are live)

- **MB-D3-1:** the receipted ω² is EXACTLY even in k; a parity-odd
  k-linear insertion has odd part exactly 2g₁k ≠ 0 — any k¹ term
  breaks the receipted k-parity (D1-RD1-5 + D2-MB-D2-2 chain),
  detected.
- **MB-D3-2:** an antisymmetric stiffness insertion breaks S = Sᵀ
  (2A ≠ 0) — the complex-ω generator is caught by the exact
  realness receipt.
- **MB-D3-3:** isotropic ε (ε = cI) kills the direction dependence
  EXACTLY (k̂ᵀMk̂ → 25c, d/dθ ≡ 0; T → 0) — the receipted anisotropy
  is carried by ε's anisotropy alone (uniaxial: d/dθ = −8(ε∥−ε⊥)
  sinθcosθ ≠ 0): signal, not artifact.

## Repairs paid (ledger)

Drift review: D3 verdict (CONDITIONAL PASS, 4acc20c1).
**R1 REQUIRED** — the RD3-2 cone-divisibility check as coded was
VACUOUS: `cancel(res/cone)·cone − res` is 0==0 for ANY residual
(the cancel round-trip); drift demonstrated it passes for
non-divisible input. The mathematical content itself HOLDS: drift
hand-ran sp.Poly.div in (kx, ky, kz) — both remainders EXACTLY
ZERO. Paid: RD3-2 re-typed to the structure tier (symmetry,
realness, kinetics); divisibility moved to **RD3-2b** — genuine
quotient/remainder form: sp.div of each residual by (1 − |k̂|²) in
(kx, ky, kz) leaves remainder EXACTLY ZERO, quotients kept and
verified by substitution back. New **MB-D3-4**: injecting a
non-divisible kx term (drift's counterexample class) leaves
remainder exactly kx ≠ 0 — the repaired form detects what the
vacuous form passed. run_fbd3 rerun: 10 assertions (6 identity +
4 mutations), exit 0. A7 fold: faithful, no repair.

## D3 verdicts (PROPOSED — drift review requested)

- **FB-D-waves (formula tier):** the director-wave dispersion is
  DERIVED: real, second-order, k² gradient structure + receipted
  pre-stress gap + polarization-degenerate coupling readout, with
  the explicit direction-dependence formula above. NO kill fired:
  realness holds structurally; K_n > 0 banked at D1. **Kill-(iii)
  kinetics half is NOW SUPPLIED** (real second-order ω² ∝ [k²-form]
  at small k); its closure still awaits the fireability condition
  (p-route independent of the wave measurement) — named, not
  claimed; stays open for D4/owner.
- **D4 window input (stated, not verdicted):** ω²_± ≥ 0 requires
  (a) k²-coefficient ≥ 0: K_n + 2C_c(k̂ᵀMk̂) ≥ 0, and (b) gap ≥ 0:
  K_p λ_± ≥ 0 at k → 0. Negative modes under admissible ε = the
  kill-(i)/window-shrink question — D4's object with its own
  frozen design.
- STOP rules: clear (principal frame = kinematics; coupling nonzero
  without objectivity-breaking terms — STOP-2 n/a; K_n > 0 banked —
  STOP-3 n/a).
- Fences unchanged: model-level family build; no carrier/electron
  claim (the RC5-analog readout is a coupling-constant handle on
  the medium side); p fixed; achiral declared; no measurement; no
  LANE-1 membership. Model-level status stated per round (K_n, C_c
  numbers ride the declared profile; formulas are exact in the
  declared model class).

## Next

D4 window + verdicts per scope 3: total 2nd-variation PSD at long
wavelength for p ∈ [0, 1) — the (a)/(b) conditions above scanned
over the admissible ε window; coupling-induced negative modes =
kill or named window-shrink with mechanism; verdicts per the 01
F-B column paths extended with FB-D items; what-licenses-what
stated; fences restated; drift review requested on the full claim
set. P3-B-dyn per-family item lands with D4 (registered scope 5).
