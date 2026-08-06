# P225 — sine-Gordon spectral induced gravity

**Campaign issue:** vantasnerdan/substrate-framework#1 — *Can the framework
construct a coherent path from sine-Gordon dynamics to general relativity?*

**Answer: yes.** This proposal contains one end-to-end, framework-native
construction. Every link is either an accepted claim of release v0.159.0, a
named standard import with explicit role, or a declared premise carried in the
debt ledger. The previously free dimensionless coordinate of C-GRV-001 is now
a computable pure number of the upstream sine-Gordon model.

## The chain

1. **Microscopic dynamics (accepted).** The physical sine-Gordon medium of
   C-MED-003: coefficients (λ, T, μ) with scales c = √(T/λ),
   ω₀ = √(μ/λ), ℓ = √(T/μ), E_scale = √(Tμ), J_scale = √(λT). The excitation
   spectrum is the accepted tower: the kink–antikink pair (normalized
   mass 8) and the C-SG-007 breather lattice E_n = 16·sin(nh/16), n·h < 8π —
   whose n = 1 level is the phonon (the Dashen–Hasslacher–Neveu
   identification), so the light quantum is never double-counted.

2. **Emergent geometry (accepted kinematics + declared premise).** The
   geometry felt by excitations is the medium's Gordon geometry (C-GOR-001,
   mostly-plus; C-GOR-002 composition). That this effective metric is the
   spacetime metric is the campaign's declared emergence premise.

3. **The missing link — the coupling (new, C-SIG-001).** One-loop vacuum
   polarization of the accepted tower on that geometry induces an
   Einstein–Hilbert term through the Seeley–DeWitt a₂ coefficient (imported
   standard result, role explicit):

   > (1/G)_s = (1/6 − ξ)·m_s²·log(Λ²/m_s²) / (2π) per species.

   The tower edge 8π is the spectrum's own ultraviolet boundary — the
   heaviest admissible state — so every logarithm is positive and finite and
   **no external cutoff is introduced**. With the quantum premise
   J_scale = ℏ (declared, following C-SCL-001's conditional precedent), the
   induced coupling takes exactly the C-GRV-001 monomial form

   > **G = q·ℓ²c³/ℏ,  q = 1/S,  S = Σ_s mult_s·(1/6−ξ)·E_s²·log((8π/E_s)²)/(2π)**

   At the canonical point (h = 1, ξ = 0, kink multiplicity 2):

   > **S = 120.058077992,  q = 0.00832930209048**

   — a pure number of the upstream model. C-GRV-001 proved dimensions alone
   cannot fix q; P074/AS3 left it a free parameter. Here it is fixed by the
   field content the framework already accepts.

4. **Gravitational dynamics (accepted).** With κ = 8πG/c⁴ (named convention),
   the C-STG-001/002 Einstein–scalar sector supplies the dynamics; the vacuum
   exterior of its spherical reduction is Schwarzschild (f = 1 − 2m/x).

5. **Static weak-field prediction (C-SIG-002).** Newtonian potential
   Φ = −GM/r and light-deflection angle 4GM/(c²b) (standard import, explicit
   role), evaluated with the **derived** G — not a declared one.

6. **Radiative prediction (C-SIG-003).** The accepted quadrupole machinery
   (C-GW-001/004) applied to the exact breather longitudinal moment
   (C-SG-009/010): cycle-averaged power

   > ⟨P⟩ = (2G/15c⁵)·[(E_scale/c²)·ℓ²·ω₀³]²·⟨μ‴²⟩,

   with ⟨μ‴²⟩ = 379.46463806875 (C-SG-010 bracket midpoint). The conditional
   input "A = 2G" of C-GW-004 is discharged by the derived coupling. At the
   unit point (λ = T = μ = ℏ = 1): ⟨P⟩ = 0.4208238630.

## Why these candidates

Candidate B (optical-dilaton direct) is structurally rejected: 1+1 gravity
has no TT sector, so no radiative prediction is possible (attempt 0001).
Candidate C (elastic collective-coordinate gravity) is structurally rejected:
scalar elasticity yields no spin-2 TT structure (attempt 0002). Candidate A
was selected on the frozen criteria — it is the only registered concept that
reaches both a static and a radiative prediction while composing the accepted
sectors without redefinition (attempt 0003).

## Verification

- `verify.py`: 21 named checks, all passing, exit 0. Includes mutation
  sensitivity (ξ → 1/6 or −1/3, kink removal, h mutation), an independent
  float-summation route agreeing to 1e-10 relative, and exact symbolic
  identities for the C-GRV-001 monomial form and the q-reciprocity.
- Package tests: `tests/test_sg_spectral_induction.py` (13) and
  `tests/test_sg_gravity_chain.py` (10), all passing.
- Full repository replay: `scripts/validate.sh` (1929 baseline tests + the 23
  new tests + governance validation).

## Debt ledger (declared premises, carried visibly)

1. Quantum premise J_scale = ℏ — declared, per C-SCL-001's conditional
   precedent. h remains the declared C-SG-007 increment; nothing derives it.
2. Emergence premise — the Gordon geometry is taken as the spacetime metric;
   C-GOR-001 supplies kinematics, and this identification is the campaign's
   declared step beyond it.
3. Semiclassical premise — each tower species contributes to the loop with
   its declared multiplicity (kink–antikink pair = 2, exposed for mutation).
4. Regulator choice — the tower edge 8π is the cutoff; the logarithm's
   additive scheme constant is fixed by that identification (no floating
   scheme parameter).
5. Convention — κ = 8πG/c⁴ named as the C-STG-001 convention factor.
6. Confrontation with observed constants is out of scope (comparators were
   blinded throughout); the model's free upstream scales are not identified
   with any observed system.

Independent-review findings (all low severity, none blocking): the C-SG-010
bracket bounds are literals mirroring the accepted claim's prose (no
importable C-SG-010 artifact exists; a revision would not propagate
automatically); `verify.py`'s exact symbolic arithmetic is slow under load
(minutes); the open-domain level count was made exactly `ceil(8*pi/h) - 1`
in response to the review.

## Scope

Claims C-SIG-001..003 are **proposed** (verification: symbolic/formal +
numeric evidence; review: unaudited pending independent review;
compatibility: compatible extension; epistemic: proposed). The accepted
registry is untouched; adjudication and release promotion are the
maintainer's gates.
