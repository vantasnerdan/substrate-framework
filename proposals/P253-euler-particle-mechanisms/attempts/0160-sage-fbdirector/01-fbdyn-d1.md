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
- **RD1-4 (p-scaling):** 𝒦(p) = p² — K_n(0) = 0 IDENTICALLY, derived
  FROM the round-declared pair premise U ∝ p² (the [n·n′]² pair
  structure — DECLARED AT D1, not pre-frozen in 00; delta (iii
  below): an unpolarized medium has no director to bend. K_n
  strictly increasing on (0, 1]; positivity holds for ANY f ≥ 0
  with M₄ > 0.
- **RD1-5 (k-structure):** the only O(k) scalar is the twist
  n̂·∇×n̂ — receipted PARITY-ODD on a concrete axial tilt field
  (T = α + O(α³) flips sign exactly) — while the receipted |∇n̂|² is
  parity-even. The achirality declaration (scope 2.1) therefore
  excludes EVERY linear-in-gradient term structurally: W's gradient
  structure is k² at D1. (R3, per drift review c5b7b854: kill (iii)
  is thereby NARROWED, not closed — its closure awaits D3 kinetics
  (I_n realness) plus independent-p measurement fireability.)

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
  from the declared pair premise and W's k² gradient structure
  structural. NO kill fired: sign positive under the declared model;
  no negative mode FROM THE STIFFNESS TERM at D1 order — the
  spectrum/PSD verdicts are D3/D4 objects (pre-stress and coupling
  contributions live there); the ω² statement and kill-(iii)
  closure await D3 kinetics + fireability. (R3 rewording per
  c5b7b854 — D1 claims the stiffness, not the spectrum.)
- Scope deltas (DATED per the freeze's amendment rule — stamped
  2026-09-12, content commits 23691391/this repair):
  (i) 2026-09-12T03:55Z, 23691391 — the gradient expansion carries
  the profile's FOURTH moment M₄ = ∫s⁴f (the scope §2.3 "∫s²f"
  wording is read as the generic profile-moment placeholder — the
  frozen items (p-scaling, sign, k², coefficient-formula) are
  unaffected; no falsifier touches the moment's power);
  (ii) 2026-09-12T03:55Z, 23691391 — the expansion identity holds
  with the sphere's induced metric M(a,b), reducing to δᵢⱼ at the
  D3 linearization point (scope's "transverse gradient metric"
  reading confirmed at that point);
  (iii) 2026-09-12T04:55Z, this repair (R1) — the pair premise
  U ∝ p² ([n·n′]² structure) was DECLARED AT D1, not frozen
  pre-compute; 𝒦(0) = 0 is derived FROM it. Tier: round-level
  model declaration, same status as the kernel profile. No other
  deltas.
- Fences unchanged: model-level family build; no carrier/electron
  claim; p fixed; achiral declared; no measurement; no LANE-1
  membership. Inertia κ_n pricing footnote and RB10/RB10b caveat

## Repairs paid (ledger)

Drift review: attempts/0108-drift-critique/review-sage-fbd1.md
(CONDITIONAL PASS, c5b7b854). Paid this commit, wording/ledger level
only — no formula touched: **R1** pair premise U ∝ p² labeled as a
D1 declaration (RD1-4 bullet + delta (iii)); **R2** deltas (i)/(ii)
dated per the freeze's amendment rule + delta (iii added);
**R3** spectrum overreach RETRACTED — D1 claims the stiffness, not
the spectrum: "no negative mode FROM THE STIFFNESS TERM at D1
order", kill (iii) NARROWED not closed (RD1-5 bullet, verdict block,
and the receipt's RD1-5 label reworded; run_fbd1 rerun green,
unchanged 8 assertions). Bankable content (K_n formula, p-scaling,
k¹ exclusion) stands per the review.

## Next

D2 coupling round: enumerate the symmetry-allowed (ε, n̂, ∇n̂) scalars
(scope 2.4's rule), extract coefficients, receipt joint-rotation
covariance of the TOTAL form + uniform-n̂ reduction to F-B-lite.
