# P247 attempt 0003 — gate W2: candidate-B enumeration and selection

Preregistered scope: `manifest.yaml` `gates.W2_candidate_B_enumeration`
(commits a8ae525/a82028e). Property checks: `w2_property_checks.py` →
`w2-property-checks.json` (all pass). This document is the paper
deliverable: the enumeration, the frozen-criteria scoring, and the single
selection for attempt 0004.

## 1. Sector and mechanism (from accepted structure, not new numerics)

The conditional P239 action (`src/substrate_framework/m5_covariant_action.py`,
projector-σ kinetic + spectral-Cartan curvature + pinned spectrum potential)
has an exact internal boost symmetry: similarity transformations of the order
parameter change the timelike eigenline u at fixed spectrum, and the spectrum
potential is blind to them. C-M5S-001 records the consequence: exactly 3
propagating species (the timelike-eigenvector boost orbit and partners), all
exactly massless. The 4 stiff static directions carry vanishing quadratic
kinetics — they are static-spectrum directions, not propagating ones.

This explains the attempt-0001/0002 structure mechanistically:

- **Isolation no-go mechanism**: the de-boxed branch's curvature class C[c(R)]
  grows without saturation because the clock/tail sector is exactly massless
  (algebraic tails, no localization scale). The isolation failure is the
  massless-boost-sector failure, not an accident of the radial ansatz.
- **Stability burden**: the q-channel negative mode lives in the stiff static
  sector, which carries NO quadratic kinetics (C-M5S-001). A static saddle
  with vanishing kinetic metric is not a dynamical (exponential) instability;
  dynamical stability of the clock is a statement about the propagating boost
  sector. This is why candidate-B terms that vanish on the static 3×3 sector
  can still deliver a fully stable isolated clock: they leave the static
  q-mode untouched (a static degeneracy) while fixing the dynamics.

## 2. Sector-vanishing term class (enumeration through O(X⁴) ⊕ O(∂²))

Terms admissible for candidate B must be local, Lorentz-covariant (covariant
under the action's own structure: spacetime η, internal metric h(P_t)), and
vanish identically on the complete arbitrary static 3×3 sector (aligned branch
u = ξ = e₀, ∂_t = 0, arbitrary static spatial block M_ij). Every admissible
term therefore contains at least one factor of the boost-sector data
(u − alignment, ∂u, M_{0μ}-block, τ-scalar) — enumerated:

| # | Term | Module primitive / status | Static-sector vanishing | Bounded below? | New continuous parameters | Asymptotic effect | Verdict path |
|---|------|---------------------------|-------------------------|----------------|---------------------------|-------------------|--------------|
| B1 | Boost-alignment mass `s_m = (m²/2)[(uᵀηξ)² − 1]` | new; exact property checks P1 (this attempt) | yes (exact: = 0 on sector) | yes (= sinh²χ · m²/2 ≥ 0) | 0 — mass fixed by convention m² = Λ² (C-M5S-002 selected scale) | masses the boost orbit → Yukawa clock tail e^{−Λr}; cures the named isolation mechanism | SELECTED |
| B2 | Candidate-H τ-current `−κ/2 (∂τ)²`, τ = −uᵀMu | `scalar_current_lagrangian_density` | yes (exact: τ ≡ 0 and ∂τ ≡ 0 on sector, check P2) | yes (Hamiltonian primitive positive) | 2 (κ, α) | kinetic-only; no mass; no cure of isolation alone | companion, not lead |
| B3 | Boost-Skyrme densities ‖[∂u, ∂u]‖-type | new (curvature-family) | yes (∂u = 0 on sector) | yes | ≥ 1 | short-distance only; massless (leading power 4, check P3) | rejected: fails natural-asymptotics |
| B4 | Candidate-K axis-lock `ζ(Tr(Y²) − Tr(P_N Y)²)` | `auxiliary_clock_axis_lock_potential` | not exact on the full arbitrary sector (locks a chosen axis against the whole field, not only the boost block) | yes | 1 | neither mass nor tail cure established | rejected: fails sector-vanishing (pre-numerics) |
| B5 | o₀-block potential monomials `X_{0μ}X^{0μ}`, cubics | new | yes on the block sector | NO (Lorentz contraction indefinite) | 1 | quadratic candidate but unbounded | rejected: fails boundedness |
| B6 | exponential matter factor `exp(2α(τ−g))` prefactor | `exponential_matter_factor` | yes (τ ≡ 0 on sector) | yes | 1 (α) | rescales existing asymptotics; does not create a scale | rejected: fails natural-asymptotics alone; retained as multiplier option |

## 3. Scoring under the frozen proposal criteria (pre-comparator, structural)

Frozen criteria (proposal.yaml `selection_criteria`), applied to the only
lead survivor B1:

1. **Bounded Hamiltonian**: s_m ≥ 0 exactly (equals (m²/2)sinh²χ); the clock
   isorotation kinetic is untouched. PASS.
2. **Local Lorentz covariance**: uᵀηξ is a full contraction; the term is
   covariant under the residual symmetry preserving the alignment ξ. The
   internal boost symmetry is explicitly broken — the assumption cost, named
   in §5. PASS with declared cost.
3. **Exact 3×3/static recovery**: s_m ≡ 0 on the entire aligned static
   sector (machine-checked, P1). Every accepted static claim survives
   unchanged. PASS (exact).
4. **Natural finite-energy asymptotics**: the massless boost orbit acquires
   the gap m = Λ; the clock tail becomes Yukawa-localized; the E(R) growth
   mechanism (massless algebraic tail) is removed at its named cause. PASS.
5. **Assumption and parameter economy**: zero new continuous parameters —
   m² = Λ² reuses the C-M5S-002-selected splitting scale (a stated
   convention, not a fit to any comparator). PASS.
6. **Topology and symmetry fit**: the clock isorotation symmetry (spatial
   internal rotations) is exactly preserved — machine-checked (P1 residual
   0) — so the committed ω = 1/(2I), E_J = E_stat + J²/(4I) clock structure
   carries over unchanged. PASS.
7. **Implementability**: in the reduced radial ansatz the boost angle is ONE
   additional radial profile χ(r) (rapidity along the director direction);
   the existing x-space modal machinery extends by one channel. No new
   solver class. PASS.
8. **Predictive reach**: new observables — clock-sector mass gap, Yukawa
   tail length, modified two-clock exchange (e^{−Λd} against the accepted
   confined d^{−1.696} boxed reference as a structural contrast, not a fit).
   PASS.

**Selection: B1, the boost-alignment mass, as the candidate-B minimal
extension.** B2 is retained as a companion term only if attempt 0004 needs
τ-dynamics for the two-clock sector; it is not part of the minimal delta.

## 4. Q-channel and asymptotic signatures (W2 gate requirement)

- **q-channel signature**: ZERO. s_m vanishes identically on the static
  sector and its quadratic expansion contains only boost angles χ, so the
  static 3×3 Hessian (including the q-channel negative mode) is unchanged.
  The stability burden transfers to dynamical (pencil) stability of the
  rotating clock, which is exactly gate W1's machinery applied to the
  extended model in attempt 0004.
- **Asymptotic signature**: boost-sector gap m = Λ; Yukawa localization of
  the clock tail; the curvature-class growth mechanism (massless tail) is
  named and targeted. The static de-boxed family's no-go stands (it is a
  statement about the unextended class).

## 5. Declared assumption cost and consumer impact

- **Named assumption cost**: explicit breaking of the internal boost
  symmetry (the alignment ξ is a fixed background covector, exactly as the
  conditional action already pins the vacuum eigenline g; the breaking
  pattern mirrors the accepted alignment selection, and does not touch
  spacetime covariance).
- **Downstream consumers**: C-M5S-001 (massless census) and C-M5S-003
  (z_i = 0 collapse) are accepted statements about the UNEXTENDED action and
  keep their scope; the extended class is a candidate-B proposal model whose
  mass-gap consequence (z_i = m²/Λ² ≠ 0) changes the induced-coupling
  evaluation for that class only. No accepted claim is reinterpreted. The
  C-M5S-008 boxed pair interaction (confined class) likewise keeps its scope;
  the isolated two-clock law is the new object.
- **Non-goals carried**: no fitted constants (m² = Λ² is a convention); no
  reinterpretation of the static no-go; no PR before #178 completes.

## 6. Design burden handed to attempt 0004

1. Extend the radial ansatz by the boost profile χ(r) (rapidity along the
   director), reusing the Chebyshev modal machinery; the s_m density in the
   reduced variables is (Λ²/2)sinh²χ with kinetic coupling through the
   projector-current metric.
2. Re-derive the fixed-J reduction: inertia I unchanged to leading order in
   χ (isorotation invariance); the clock tail now sources the massive boost
   field — verify finite E and localization.
3. Restability: order-24 pencil spectrum of the rotating background (W1
   machinery) — the clock is stable iff the boost-sector pencil is
   tachyon-free and the static q-degeneracy stays kinetic-free.
4. Two-clock asymptotics: massive exchange → derive the pairing from the
   linearized massive field (small-ratio-numerics asymptotic-matching
   prescription; carry tail amplitudes in log space).
