# 01 — FBDYN D1: director stiffness K_n DERIVED (bankable round 1)

Scope: 00-fbdyn-scope.md (frozen pre-compute, 92825c5e). Round D1 per
that plan: the director stiffness derived via the frozen
gradient-expansion procedure. Receipts: receipts/run_fbd1.py +
run_fbd1.log (8 assertions — 5 identity + 3 mutations — exit 0;
negative controls included). This round bankable alone: D2–D4 proceed
on this K_n; a D2 zero-coupling finding would leave the
free-director propagation claim intact (scope 3).

## Derived structure (each step receipt-backed)

- **RD1-1 (expansion identity):** [n(x)·n(x+r)]² = 1 − M(a,b)
  (r·∇a)² + O(4) with M = 1 + a²/(1−a²−b²) the sphere's INDUCED
  metric (positive definite; = 1 exactly at a = b = 0 — the D3
  linearization point about uniform n̂₀). No linear-in-r term exists
  (the unit constraint kills it pointwise). Bending reduces pair
  alignment QUADRATICALLY, positive coefficient — the k² structure's
  source. Sympy-verified on the exact-unit parametrization
  (split-sqrt form; the combined sqrt's Abs artifacts noted and
  avoided).
- **RD1-2 (isotropy):** ∫dΩ ŝᵢŝⱼ = (4π/3)δᵢⱼ EXACT → the leading
  stiffness is isotropic, (K_n)ᵢⱼ = K_n δᵢⱼ. Direction-dependence of
  the stiffness is a priced higher-order effect (stated, not claimed
  — scope 3 D1).
- **RD1-3 (coefficient formula):** with the declared kernel anchor
  g(r) = (K/ξ³) f(r/ξ) and w_bend = (K_n/2)|∇n̂⊥|²,
  **K_n = (8π/3)·K·p²·M₄·ξ²**, M₄ = ∫₀^∞ s⁴f(s)ds. Declared profile
  f = e^{−s}: M₄ = 24 exactly → K_n = 64π K p² ξ². Gaussian:
  M₄ = 3√π/8. Sign and the p²ξ² structure are profile-independent;
  the NUMBER is model-level (F-A class honesty, stated).
- **RD1-4 (p-scaling):** 𝒦(p) = p² — K_n(0) = 0 IDENTICALLY (the
  frozen requirement is met by derivation, not assumption): an
  unpolarized medium has no director to bend. K_n strictly increasing
  on (0, 1]; positivity holds for ANY f ≥ 0 with M₄ > 0.
- **RD1-5 (k-structure):** the only O(k) scalar is the twist
  n̂·∇×n̂ — receipted PARITY-ODD on a concrete axial tilt field
  (T = α + O(α³) flips sign exactly) — while the receipted |∇n̂|² is
  parity-even. The achirality declaration (scope 2.1) therefore
  excludes EVERY linear-in-gradient term structurally:
  **ω² ∝ k² at small k is derived at D1** (falsifier kill (iii)
  closed structurally before D3).

## Mutations (negative controls — the kills are live)

- **MB-D1-1:** antialigned correlation (G < 0) flips K_n < 0 — kill
  (ii) of FB-D-waves TRIGGERS. Positivity is a real claim riding the
  correlation sign, detectable, not structural.
- **MB-D1-2:** director waves at p = 0 are FORBIDDEN — the stiffness
  is exactly zero there; any propagation claim in the unpolarized
  medium contradicts the receipted scaling and is caught by the
  identity.
- **MB-D1-3:** a linear-in-gradient coefficient rides the parity-odd
  twist scalar — W changes under inversion unless the coefficient is
  zero (antisymmetric part = 2c₁T exactly, = 2c₁α at linear order).
  The receipted W has c₁ = 0.

## D1 verdicts (PROPOSED — drift review requested)

- **FB-D-stiffness:** the dynamical director has a DERIVED stiffness
  K_n = (8π/3)K p² M₄ ξ² > 0 on p ∈ (0, 1], with 𝒦(0) = 0 derived
  and the k² structure structural. NO kill fired: sign positive
  under the declared model, no negative mode, ω² ∝ k² derived.
- Scope deltas (dated, per the freeze's amendment rule): (i) the
  gradient expansion carries the profile's FOURTH moment M₄ = ∫s⁴f
  (the scope §2.3 "∫s²f" wording is read as the generic
  profile-moment placeholder — the frozen items (p-scaling, sign,
  k², coefficient-formula) are unaffected; no falsifier touches the
  moment's power); (ii) the expansion identity holds with the
  sphere's induced metric M(a,b), reducing to δᵢⱼ at the D3
  linearization point (scope's "transverse gradient metric" reading
  confirmed at that point). No other deltas.
- Fences unchanged: model-level family build; no carrier/electron
  claim; p fixed; achiral declared; no measurement; no LANE-1
  membership. Inertia κ_n pricing footnote and RB10/RB10b caveat
  travel with every citation.

## Next

D2 coupling round: enumerate the symmetry-allowed (ε, n̂, ∇n̂) scalars
(scope 2.4's rule), extract coefficients, receipt joint-rotation
covariance of the TOTAL form + uniform-n̂ reduction to F-B-lite.
