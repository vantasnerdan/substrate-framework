# 07 — F-B-LITE build: polarized-tangle elasticity (static director)

Charter: shepherd F-B-lite, owner option A (pricing 150c795d PASS).
SCOPE: static n̂ ONLY — n̂ is a fixed model parameter. The DYNAMICAL
director (own stiffness, ε·∇n̂-coupling) is OUT OF SCOPE: separate priced
lane, unchartered here. Registered dues RIDE this round (all paid — §5).
Receipts: receipts/run_fb.py + run_fb.log (13 assertions self-counted —
10 identity + 3 mutations, exit 0). Verdicts below are PROPOSED and go
to drift firewall per standing discipline; 01's F-B row carries them as
landed-pending-review.

## Model assumptions (stated, all model-level; travels from F-A class)

1. Polarized tangle: F-A's family with a polarization axis n̂ (unit) and
   scalar order parameter p ∈ [0, 1] in the segment-direction averages:
   ⟨nᵢnⱼ⟩ = (1−p)δᵢⱼ/3 + p n̂ᵢn̂ⱼ; ⟨4-fold⟩ = (1−p)·iso4 + p n̂⊗n̂⊗n̂⊗n̂.
2. ACHIRALITY DECLARED (FB-4 requirement): the averages carry no twist
   structure. A chiral polarization would make FB-4 substantive and is
   NOT built here.
3. Static director (scope boundary): n̂ fixed; no director stiffness, no
   ε·∇n̂ coupling.
4. Cutoff convention K = (Γ²/4π)L₀ln(ℓ/a): stated, travels from F-A
   (including its slaving-sensitivity pricing — RS).
5. p is a FREE parameter: what sets it (tangle self-consistency) is
   unpriced dynamics — named limit, not hidden.

## Derived structure (each step receipt-backed)

- RB1 (continuity): p = 0 recovers the F-A form exactly — F-B-lite is
  the p-extension of the established family, not a new construction.
- RB2 + RB3 (DUE-1 paid): the quadratic form is transversely isotropic
  WITH the relation C₁₂ = C₁₃ (cross coefficients extracted and compared:
  both K(p−1)/15, nonzero) — 4 distinct constants, one relation short of
  general TI. The relation is a falsifiable statistics-level signature.
- RB4: exact traceless 2nd-variation spectrum {K(1−p)/10, K(1−p)/5,
  K(3p+2)/10} — linear in p ⇒ analytic PSD for ALL p ∈ [0,1]; rational
  p-grid confirmation.
- RB5: μ⊥ = K(1−p)/10 (transverse shear), μ∥ = K(3p+2)/20 (axial-plane
  shear). **FB-anisotropy falsifier formula banked:
  μ∥/μ⊥ = (3p+2)/(2(1−p)).**
- RB6 (DUE-3 paid, CORRECTS the pricing doc): on traceless strains the
  normal part decomposes as
      W|traceless-diag = K·p·e₃₃  +  K(1−p)/10·(e₁₁²+e₂₂²+e₃₃²).
  The pricing doc's "K(1−p)/30" was WRONG (receipt exposed it; corrected
  in 06 inline). NEW physics the receipt surfaced: the anisotropic
  LINEAR pre-stress K·p·(n̂·εn̂ − tr ε/3) — vanishes at p = 0, i.e. the
  polarized reference state carries an axis-coupled pre-stress at this
  order (bookkeeping: shifts reference stress, does not affect moduli —
  moduli are 2nd-variation objects).
- RB7 (DUE-2 paid): μ∥(p=1) = K/4 exactly (real assert; repairs the
  scratch's truncated print + tautology).
- RB8: p = 1 spectrum {0, 0, K/2}: the degeneracy is a sliding mode
  (aligned filaments free to shear transversely to quadratic order) —
  zeros, never negatives.
- RB9 (objectivity, FB-2): joint covariance W(RεRᵀ, Rn̂) = W(ε, n̂) —
  exact, rational rotation. MB2: rotating the strain WITHOUT the
  director breaks covariance — the axis pins the frame; objectivity
  requires the joint rule.
- RB10 (FB-5, MATERIAL statement — PRE-ACOUSTOELASTIC): direction-
  dependent MATERIAL shear stiffness from the same W:
  W(shear wave, k at angle θ from the axis)
  = K[(1−p)sin²θ/20 + (3p+2)cos²θ/40]. CAVEAT (drift 749179cd,
  REQUIRED, added): the polarized reference carries the RB6 pre-stress
  σ₀₃₃ = +Kp — tensile, uniform, self-equilibrated — whose geometric
  (incremental-moduli) coupling STIFFENS transverse shear waves: wave
  stiffness = material + prestress correction. RB10b receipts the
  correction's SIGN (tensile ⇒ stiffens upward; no shear linear terms —
  axis-diagonal only); its SIZE is O(p) and UNPRICED (incremental
  moduli = a new receipt round if chartered). Signs/nonnegativity
  ROBUST; quantitative wave speeds shift O(p).

## FB verdicts (per 01, F-B column; PROPOSED pending drift)

- FB-1 (anisotropic static-shear branch): FORMULA half ESTABLISHED at
  stated model order — μ⊥ > 0 for p < 1 with the p-window EXPLICIT;
  p = 1 is the sliding DEGENERACY (stated, not hidden; zeros never
  negatives). Window half PRICED (travels from F-A: ω ≫ 1/τ, τ priced
  not derived). NO static μ(0) claim (viscoelastic honesty bar).
- FB-2 (objectivity): ESTABLISHED (RB9 + MB2 — the mutation shows the
  content: joint rotation is load-bearing once an axis exists).
- FB-4 (helicity inheritance): trivial-expected AT THE DECLARED ACHIRAL
  TIER — ⟨H⟩ = 0 needs achirality once the ensemble is polarized;
  declared in assumptions with the re-verdict clause standing (chiral
  branch ⇒ substantive ⇒ gate-1 machinery).
- FB-5 (propagation): direction-dependent MATERIAL shear dispersion
  RECEIPTED (RB10) — second order, real, nonnegative for p ∈ [0, 1)
  (degenerate at 1); conditional on the F-A window as FB-1. CAVEAT
  (drift 749179cd): qualitative content SAFE (tension stiffens upward —
  signs robust); quantitative speeds carry O(p) acoustoelastic
  corrections (geometric terms unpriced) — comparisons against measured
  speeds must budget them.
- FB-6 (scale separation): HOLD-conditional, travels unchanged.
- SYN P3 (per-family amendment): the polarized branch's fluctuation
  SPECTRUM differs from the random-pilot prediction — the amendment
  obligation is OWED and stays with SYN's owner. F-B-lite licenses NO
  lane-level statement that mixes F-A and F-B claims before that
  amendment.

## Falsifier FB-anisotropy (pre-registered; STOP FROZEN PRE-COMPUTE)

- Formula (banked RB5): measured shear-stiffness ratio must satisfy
  μ∥/μ⊥ = (3p+2)/(2(1−p)) with p determined INDEPENDENTLY of the ratio.
- KILL conditions: (i) measured ratio inconsistent with the formula
  beyond stated noise AND beyond the priced-acoustoelastic budget, with
  independently-fixed p; (ii) any measured NEGATIVE stiffness in the
  shear sector (model predicts PSD — a negative mode is a model kill,
  not a window note).
- HYGIENE (drift 749179cd, REQUIRED): a MISS within
  geometric-correction size (the unpriced O(p) acoustoelastic budget of
  RB10/RB10b) is INCONCLUSIVE — model-extendable, NOT a kill. Recorded
  alongside the ratio formula before any measurement charter. The
  material-ratio formula itself is unaffected (RB5 stands as the
  material statement).
- STOP rule (frozen before any measurement exists): the falsifier is
  FIREABLE only when p is measured by a route that does not use the
  ratio it predicts. If no such independent-p route can be constructed,
  declare FB-anisotropy UNFIREABLE (decoration, per SYN-spec discipline)
  and downgrade FB-1's falsifier status accordingly — no silent
  fire-and-miss.
- Measurement host: F-C-hosted readout analog (direction-dependent tilt/
  shear spectroscopy of the crystal/polarized medium); design owed —
  NOT part of this round (no measurements were made or claimed here).

## What F-B-lite does NOT claim

No dynamical director (out of scope; separate priced lane). No p
self-consistency. No isotropic LANE-1 membership; no lane conjunction
fires. No carrier/electron claim. No measurement. Chirality excluded by
declaration, not derived.

## Family state (post-caveat, shepherd routing)

F-A, F-C, F-B-lite ALL BUILT (static-director tier); F-B dynamical
director = SEPARATE PRICED LANE (unchartered). SYN P3 amendment remains
OWED (unchanged). RB10/RB10b caveat travels with every citation of the
anisotropy formula.

## §5 Registered dues — ALL PAID this round

1. C₁₂=C₁₃ extraction line → RB3 (extract-and-compare executed; the
   "or soften" alternative not needed — the relation is receipted).
2. pricing-6 print repair → RB7 in run_fb.py AND the scratch repaired
   in place (real assert μ∥(1)=K/4; tautology removed).
3. Normal-part receipt → RB6, which CORRECTED the value: K(1−p)/10·Σeᵢ²
   (quadratic) + K·p·e₃₃ (linear pre-stress), NOT the doc's K(1−p)/30.
   06 carries the inline correction.
