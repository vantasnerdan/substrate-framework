# 05 — F-B frame-price analysis (EXPLORATORY PRICING, no F-B claims)

Charter: shepherd tasking 2026-09-11 — "draft the pricing (ordered-vs-
random family choice, frame cost, what F-B would need). Analysis only;
build decision stays owner-level." This doc is the analysis. The
companion scratch `scratch/fb-pricing-scratch.py` is EXPLORATORY EVIDENCE
for the cost estimates below — NOT an F-B receipt, no verdicts, no 01
row changes. 01's F-B row ("no computation") is superseded for this
tasking only; owner call restores or lifts it.

## 0. What F-B is (per 00-sketch, unchanged)

Polarized tangle / vortex lattice (drift breaker #7 "S4b ordered-vacuum",
constructive version): anisotropic elasticity; frame price NAMED
(preferred axis) — owner choice: declare anisotropic / domain-average to
isotropy / kill. Where F-A statistically averages over segment
orientations (isotropy by statistics, no preferred frame), F-B declares
a polarization axis n̂ with scalar order parameter p.

## 1. The frame price, decomposed (P1–P3)

**P1 — preferred-axis closure (the actual frame cost).** F-A's frame
was free: filaments are material (frozen-in), isotropy came by
statistics. F-B must say what n̂ IS:
- *Static director* (n̂ a fixed model parameter): closure is model-level,
  same status as F-A's cutoff convention. CHEAP.
- *Dynamical director* (n̂ a field with its own stiffness, coupling
  ε·∇n̂-type): NO derivation exists in this charter — that is a new
  mechanism lane (the analog of F-C's tilt sector for the polarization).
  This is the single biggest price separation in the owner call.

**P2 — modulus algebra (priced by scratch, exploratory).** Cost is
F-A-class: the polarized angular averages ⟨ninj⟩ = (1−p)δij/3 +
p n̂in̂j, ⟨4-fold⟩ = (1−p)·iso4 + p n̂⊗n̂⊗n̂⊗n̂ reduce the quadratic form
to a transversely isotropic form. Scratch results (axis z, exact
symbolic):
- Model-order signature: the dipole-statistics form carries 4 distinct
  constants with C12 = C13 — one relation short of general TI
  elasticity; falsifiable as stated.
- p → 0 recovers the F-A form EXACTLY (family bridge continuous).
- Shear structure: μ⊥ = K(1−p)/10 (transverse, e12), μ∥ = K(3p+2)/20
  (axial-plane, e13/e23). Anisotropy ratio μ∥/μ⊥ = (3p+2)/(2(1−p)).
- **Stability window: the traceless 2nd-variation eigenvalues are
  K(1−p)/10, K(1−p)/5, K(3p+2)/10 — PSD for ALL p ∈ [0,1]** (charpoly
  factored in scratch). No internal instability kill at dipole-statistics
  order; the window CLOSES at p = 1: μ⊥ = 0 — a degeneracy (sliding mode
  of perfectly aligned filaments), not a negative mode.
- Normal part: normal = K(1−p)/30 on traceless strains (all p).
So P2 ≈ ONE receipt round (anisotropic R1a–c analog + 2nd-variation +
objectivity-covariance check + wrong-4th-moment mutation), reusing the
F-A machinery wholesale. The interesting physics (anisotropy ratio, p→1
sliding degeneracy) comes out of the SAME receipt.

**P3 — family liabilities the frame choice commits to.**
- *FB-4 helicity*: F-A's trivial-expected used isotropy. A polarized
  ensemble keeps ⟨H⟩ = 0 only if the polarization is ACHIRAL — model-
  level declaration available, but if the owner ever wants chiral
  (twisted) polarization, FB-4 goes substantive and gate-1 machinery
  becomes load-bearing. Priced as a stated assumption for F-B-lite.
- *SYN P3 amendment* (owner-visible, 00-sketch): an ordered branch
  changes the fluctuation spectrum vs the random-pilot prediction; SYN
  P3 must be re-derived PER FAMILY. Either owner option A or B creates
  an SYN-owner obligation. Never a silent carry.
- *Domain-average option only*: adds domain-wall physics (wall energy,
  defect cost) — NOT priced here; strictly more expensive than A; and
  the resulting isotropic claim inherits double-averaging criticism
  (firewall will flag it).

## 2. Ordered-vs-random: the family choice, stated plainly

These are not competing hosts for LANE-1. F-A is built and
receipt-backed; F-B would EXTEND the medium claim to ordered vacua
(p ≠ 0), with F-A as the p = 0 boundary. What F-B buys that F-A cannot:
- A distinct falsifiable prediction at paper level: direction-dependent
  shear stiffness with ratio (3p+2)/(2(1−p)) — measurable structure, not
  a second copy of μ = K/10.
- The constructive answer to drift breaker #7 (ordered vacuum) inside
  the dipole-statistics model class.
What it cannot buy (honest limits): p is a FREE parameter in the
scratch — what SETS the polarization (tangle self-consistency) is
unpriced dynamics; no mutual-interaction renormalization; achirality
assumed; the dipole-statistics model class is the same level as F-A
(model-level, stated).

## 3. Owner options (costs; decision NOT made here)

| Option | Cost | Licenses |
|---|---|---|
| A: F-B-lite (static director, achiral, p free) | ≈ one receipt round (P2) + stated assumptions (P1-static, achirality) | anisotropic μ>0-type claim conditional on stated model assumptions; NOT isotropic LANE-1 without SYN P3 amendment |
| B: domain-average to isotropy | A + domain-wall pricing (unpriced) | weaker isotropic claim; high firewall risk (double averaging) |
| C: kill F-B | zero | F-A sole medium family; F-C unaffected (standalone) |
| D: defer, pricing banked (this doc) | zero now | preserves A later; scratch stays exploratory |

## 4. What F-B-lite would need to pass the 01 gates (F-B column)

- FB-1: μ-type > 0 with the p-window EXPLICIT (p < 1; p = 1 degeneracy
  stated as the sliding mode, not hidden) — 2nd-variation receipt.
- FB-2: objectivity of W(ε; n̂) — joint-rotation covariance receipt
  (rotates (ε, n̂) together; form invariant). F-A's RB additive-rotation
  receipt transfers at the strain-identification level.
- FB-4: achirality declaration + re-verdict clause (substantive iff
  chiral branch ever opened).
- FB-5: anisotropic shear-wave dispersion ω(k, direction) derivable from
  the same W — receipt-worthy, no new machinery.
- FB-6: travels from F-A (HOLD-conditional, cutoff pricing).
- SYN P3: amendment scheduled with SYN's owner BEFORE any LANE-level
  statement mixes F-A and F-B claims.

## 5. Pricing verdict (analysis conclusion)

The frame price is FRONT-LOADED in P1 (what n̂ is) and P3 (SYN P3 +
achirality), not in P2: the anisotropic modulus algebra is cheap,
continuous with F-A, and internally stable for all p — with a genuine
new prediction (shear anisotropy, p→1 sliding degeneracy). If the owner
wants the ordered-vacuum answer, Option A is proportionate (one receipt
round); Options B and C are priced in the table; nothing here builds F-B
or touches 01's F-B verdict row.
