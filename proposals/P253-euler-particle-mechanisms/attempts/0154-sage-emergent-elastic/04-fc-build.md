# 04 — F-C build: tilt (rotational) coupling — first gate-licensed missing-5 candidate

Charter: shepherd F-C (family order F-A → F-C; F-A formula PASS
b1301a99). Gate passage is MEASURED (derived predicates), not assumed.
Verdict paths per 01-verdicts-per-family.md. Receipts: run_fc.py +
run_fc.log (10 assertions self-counted, exit 0).

## Family (model-level honesty stated up front)

Axially-stacked ring crystal: thin vortex rings (radius R, core a,
circulation Γ) stacked at spacing ℓ along their common axis, one
rotational DOF per ring — rigid TILT of the ring plane (pendulum), not
bend. Existence-grade anchor, not membership claim: single thin rings are
classical exact-asymptotic solutions, and EXACT multi-ring Euler
constructions exist in the campaign's verified source set (S3 primary:
Dávila–del Pino–Musso–Wei Theorem 1, consumed at 0110 scope — no-swirl,
similar-circulation, reduced-law dynamics). The LATTICE model itself is a
mean-field model built from these ingredients; it is not asserted to be
an exact solution family. This is the same model-level status as F-A,
stated per discipline.

## Derived structure (each step receipt-backed)

1. Tilt stiffness (RC1, RC2): leading-multipole ring-ring interaction,
   small tilt, axial stacking — quadratic form
   (Γ²m²/4πℓ³)(θ₁θ₂cos Δφ + θ₁² + θ₂²) with m = πR². Two axial
   neighbors: **J₀ = Γ²πR⁴/(2ℓ³) > 0** — the aligned chain is
   tilt-stable. Mutation MA-3: dropping the −3(m·r̂)² multipole term
   halves-and-FLIPS the stiffness (negative) — the dipole structure is
   load-bearing for stability.
2. Vikulin strain–rotation coupling (RC3): under coarse strain,
   **J(ε) = J₀(1 + 4η_in − 3ε_zz)** (in-plane stretch enlarges the ring
   moment; axial stretch weakens neighbor coupling). This is the
   elastic-rotational coupling CONSTANT derived, not assumed — the
   Vikulin-type structure (strain field modulating a rotational pendulum)
   emerges from the family.
3. Rotation waves (RC4): coupled-pendulum chain dispersion
   **ω²(k) = (4J₀/I)·sin²(kℓ/2) ≥ 0** — second-order, real,
   nonnegative. GATE 2 passage: the emergent sector carries its own
   kinetic structure (micro-rotation inertia I), exactly the extension
   0147's M5 certificate demands — here arising mechanically from the
   family's inertia, not added by hand.
4. Carrier coupling — the missing-5 candidate (RC5, RC7, MA-4): a
   carrier filament threading the ring n times interacts through the
   tilted-disk solid angle: U_c(θ) = n·ΓΓ_c·Ω(θ)/(4π) with
   Ω(θ) = 2π(1 − d/√(R²+d²)) + Ω₂θ² + O(θ⁴),
   **Ω₂ = −πR²d(d²−2R²)/(2(R²+d²)^{5/2})**. Pendulum shift:
   **δ(ω²) = n·ΓΓ_c·Ω₂/(2πI)** — LINEAR in the integer n (RC7).
   Medium-side charge readout: magnitude AND sign (Ω₂ flips at
   d = √2·R, RC5c — relative geometry readable too). MA-4: Ω has NO
   linear-in-θ term — no force template is smuggled (Magnus family
   dead-check: the coupling enters as pendulum-potential shift on the
   MEDIUM side; the dead family was carrier-side force templates; no
   variant is used).

## 0147-gate passage (as measured)

- **Gate 1 (CS current + helicity preserved)**: rigid tilt and affine
  strain keep every ring PLANAR; planar curves have zero writhe
  (Călugăreanu–White class — cited), so Slk = framing (fixed) and
  per-ring helicity Γ²Slk is exactly invariant under the emergent
  dynamics; the derived W(ε, θ) contains no helicity-flux term.
  Measured-as: algebraic planarity statement + absence in the derived
  Lagrangian. (Cited topology flagged per discipline.)
- **Gate 2 (degeneracy deliberately broken)**: RC4 — the tilt DOF adds
  real second-order kinetics absent from the bare action (0147
  I-Noether2's certificate detected the requirement; the family supplies
  it). Measured-as: derived dispersion.
- **Gate 3 (integer class)**: RC7 — the readout is proportional to the
  integer threading number n; no continuous-fractional charge enters.
  Measured-as: linearity structure of the derived coupling.

## FB verdicts (per 01, F-C paths; drift rules)

- FB-1 (elastic branch): the crystal is ANISOTROPIC (in-plane vs axial
  differ); the tilt sector is derived here in full; the full anisotropic
  elastic tensor follows by the F-A machinery restricted to the ring
  plane — PRICED (footnote travels), not re-derived this round.
- FB-2 objectivity: the dipole form is rotation-invariant by
  construction (scalar from m̂₁·m̂₂, (m̂·r̂) products) — stated at
  construction level.
- FB-3 bookkeeping: every expansion carries its stated order (O(θ³),
  O(ε²), O(θ⁴)); no sink.
- FB-4 inheritance: SUBSTANTIVE for F-C (tilt modes touch per-ring
  helicity bookkeeping) and satisfied at model level: rigid tilt keeps
  Slk fixed (planarity), strain preserves the framing class. Re-verdict
  clause stands.
- FB-5: RC4 — in the crystal's elastic window (same viscoelastic
  honesty bar as F-A: ω ≫ 1/τ; τ priced, not derived — F-A footnote).
- FB-6: HOLD-conditional (ℓ vs a vs carrier scale; parameters stated).

## Missing-5 status after F-C

First candidate coupling that passes the 0147 gate AS MEASURED
(gate-1/2/3 each tied to a derived predicate). Pre-registered falsifier
**FB-C1** (D-08 inequality form): measured tilt-frequency shifts across
threadings must satisfy δ(ω²)(n)/δ(ω²)(1) = n + margin·(tolerance), one
geometric constant; KILL: non-proportional response or non-integer-
effective scaling. The carrier's electron-sector identity remains
SYN-conditional; nothing here claims it — this builds and gates the
OBSERVABLE STRUCTURE only.

## F-A footnotes traveling (per charter)

Cutoff/self-energy convention stated-not-derived (classical ring energy
cited); window τ priced-not-derived; inertia I scaling-level (ρℓ³R²·κ,
κ priced); FB-6 HOLD-conditional; no static claims beyond derived
orders; anisotropic elastic tensor priced.
