# 00 — F-B DYNAMICAL DIRECTOR: SCOPE FREEZE (frozen pre-compute)

Charter: shepherd 2026-09-12 (owner delegated): the priced
new-mechanism lane — dynamical director (own stiffness + ε·∇n̂-type
coupling), the object F-B-lite's scope boundary excluded (07 §3) and
the frame price named as the real cost (06 §1 P1: "the analog of F-C's
tilt sector for the polarization; NO derivation exists"). This doc
freezes scope, derivation plan, falsifiers, and STOP rules BEFORE any
compute. Any later change lands as a scoped delta with its own dated
entry — the frozen content below stays auditable.

## 1. Object

Within the F-A/F-B-lite dipole-statistics model class, DERIVE (not
assume) the dynamical polarization sector:

- (i) director stiffness K_n: the energy cost of spatially bending the
  polarization axis, from a gradient expansion of the local-average
  interaction — the analog of F-C's derived J₀;
- (ii) the lowest-order strain–director-gradient coupling (ε·∇n̂ type):
  enumerate the symmetry-allowed terms and derive their coefficients —
  the analog of F-C's derived J(ε);
- (iii) director-wave dispersion: real second-order kinetics for the
  transverse director modes — the analog of F-C's RC4 (gate-2-type
  predicate for the polarization sector);
- (iv) the coupling's back-effect on the F-B-lite static claims: PSD
  window of the total form at long wavelength;
- (v) falsifier(s) with kill conditions and STOP rules (this doc).

## 2. Frozen design decisions (made now, before compute)

1. **Field and symmetry class.** n̂(x, t): unit director, n̂ → −n̂
   symmetric (nematic-type) — consistent with F-B-lite's n̂⊗n̂
   structure and the FB-4 achirality declaration. All W terms even in
   n̂; no twist/chiral invariants (achiral declared; chiral branch
   stays closed with its re-verdict clause).
2. **Fixed p.** The director dynamics moves n̂ at FIXED scalar p.
   p-field dynamics / p self-consistency stays OUT of scope (named
   limit, travels from 07 §2.5 unchanged).
3. **Correlation model (the stiffness's declared input).** Segment
   orientations correlated over a range ξ with declared scaled profile
   f(r/ξ) (model-level, same status as every constant in this family).
   The gradient-expansion PROCEDURE is frozen: expand the interaction
   of neighboring correlation volumes about the local average, keep
   the O(ξ²|∇n̂|²) term. The derived items are: the p-scaling of K_n,
   its sign, the k² structure, and the coefficient's formula in terms
   of the profile moment ∫ s²f(s)ds — the NUMBER is model-level and
   labeled so (F-A class honesty).
4. **Coupling enumeration rule.** W_coupling = all rotational scalars
   at leading order built from (ε, n̂, ∇n̂) satisfying: objectivity
   (joint rotation of ε, n̂, ∇n̂), evenness in n̂, achirality, and
   reduction to zero at uniform n̂. The COUNT of allowed independent
   terms is itself receipt content (a falsifiable model signature).
5. **Kinetics.** Director inertia per unit volume I_n: scaling-level
   ρξ³·κ_n with κ_n priced-not-derived (F-C's inertia footnote
   pattern); the dispersion's REALNESS and k² structure are the
   claims; absolute speeds carry the priced κ_n.
6. **Kinematic constraint.** Unit-n̂ constraint handled at quadratic
   order: δn̂ ⊥ n̂₀ (transverse modes only); the constraint's
   projection is part of the D3 receipt.

## 3. Round plan (bankable intermediates; each round = doc + receipts + commit)

- **D1 stiffness.** Derive K_n = K·𝒦(p)·(moment of f)·ξ²: receipt the
  p-scaling (requirement: 𝒦(0) = 0 — no polarization, no director, no
  stiffness), sign > 0 on the live window, k² structure; the p = 1
  behavior derived and named, not hidden. Mutations: flipped
  correlation-sign flips K_n (detected); isotropic f (no polarization)
  yields exactly zero stiffness.
- **D2 coupling.** Enumerate allowed terms under §2.4's rule; extract
  coefficients; receipts: joint-rotation covariance of the TOTAL form
  (RB9 analog extended to ∇n̂), evenness, uniform-n̂ reduction to the
  F-B-lite form exactly (RB1 analog: the dynamical theory is the
  p-extension's field-theoretic completion, continuous with the built
  family). Mutation: a term breaking objectivity breaks covariance
  (detected).
- **D3 dispersion + readout.** Quadratic expansion about uniform
  (n̂₀, rest): director-wave dispersion ω²(k) = (K_n k² + ΔW_prestress)/I_n
  with ΔW_prestress the RB6/RB10b pre-stress contribution (sign
  receipted there: stiffens upward); require REAL second-order for
  p ∈ [0, 1); derive the direction-dependence formula (anisotropy of
  director-wave speed — the falsifiable structure); derive the
  coupling's k-dependent shift δ(ω²)(k) as the readout of the coupling
  constant (F-C RC5 analog: a measurable handle on g, not a new
  mechanism claim).
- **D4 window + verdicts.** Total 2nd-variation PSD at long wavelength
  for p ∈ [0, 1): coupling-induced negative modes = kill or named
  window-shrink with mechanism; verdicts per the 01 F-B column paths
  extended with FB-D items; what-licenses-what stated; fences restated.
  Drift review requested on the full claim set.

Bankable: each round lands alone (D1's K_n is meaningful even if D2's
enumeration finds zero allowed couplings — that finding itself would be
banked content: the director sector would then be FREE, an even
stronger propagation claim).

## 4. Falsifiers (FROZEN pre-compute)

**FB-D-waves** (the lane's falsifiable prediction): director waves
exist with ω²(k) ∝ k² at small k, real, with direction-dependent speed
given by the D3 formula, p fixed independently.

- KILL conditions (any one):
  (i) any NEGATIVE mode in the derived spectrum within p ∈ [0, 1) at
  declared order (the family's PSD spine — RB4/RB8 pattern);
  (ii) derived K_n ≤ 0 on the live window under the declared
  correlation model (a bent director lowering energy = model-internal
  instability);
  (iii) measured director-wave response with ω² NOT ∝ k² at small k,
  beyond the priced untracked orders, with p independently fixed.
- HYGIENE (travels from FB-anisotropy discipline): a MISS within the
  priced untracked orders (ξ/λ separation; inertia κ_n; the
  multipole/quadrupole footnote; acoustoelastic O(p) budget) is
  INCONCLUSIVE — model-extendable, not a kill.
- FIREABILITY: as with FB-anisotropy, the falsifier is FIREABLE only
  with a p-route independent of the wave measurement; if none can be
  constructed, declare FB-D-waves UNFIREABLE and downgrade the
  propagation claim's falsifier status accordingly — no silent
  fire-and-miss.

**STOP rules (frozen now):**

1. If any round's derivation requires positional co-variation
   statistics beyond the declared direction-statistics model class
   (the boundary 08-syn-p3 already named), STOP that round, bank what
   is derived, name the missing construction. No silent model
   extension.
2. If the coupling enumeration (D2) needs an objectivity-breaking term
   to produce a nonzero coupling, STOP: the lane's coupling claim dies
   with the named mechanism (symmetry), and D3/D4 proceed only for the
   free-director sector.
3. If D1's K_n comes out negative, the lane is KILLED at D1 (condition
   4.i/ii) — bank the mechanism, do not proceed to D3/D4 against a
   dead stiffness.

## 5. Registered obligations (travel with the lane)

- SYN P3 per-family: the STATIC-branch amendment (P3-B) landed for
  F-B-lite; a DYNAMICAL director adds fluctuation modes — the
  dynamical-branch spectrum is a NEW per-family item (P3-B-dyn),
  registered NOW as owed-with-D4, not silently covered by the landed
  amendment.
- Inertia κ_n pricing footnote travels with every speed claim.
- Multipole/quadrupole footnote travels (from 04) wherever K-type
  constants enter.
- RB10/RB10b caveat (acoustoelastic O(p)) travels with every citation
  of static anisotropy formulas inside the dynamical theory.

## 6. Fences (unchanged, restated)

No carrier/electron claim (carrier readout ANALOGY in D3 is a coupling
readout on the medium side; the dead force-template family is not
used — MA-4 pattern stated per round). No measurement. No LANE-1
membership; no lane conjunction fires. Achirality by declaration.
Model-level status stated per round. Drift firewall reviews every
claim round; scope deltas are dated entries, never silent edits.
