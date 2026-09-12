# 01 — Identity ledger (S3-EXACT)

Every identity is an EQUALITY (D-08 boundary class: exact analytic identities
stay equalities). Each carries: statement with domain/quantifiers, derivation,
receipt line in `receipts/s3exact-sympy/run.log` (assertion-backed; the file
exits nonzero on any failure), and the mutation that must fail the identity
property. Domain throughout: Ω_δ(t) per 03-prefire F2 spec (tubular excision,
weak Lin, core-adjacent variations excluded, distributional core terms kept).

## Conventions

u = velocity, ω = ∇×u vorticity, Π = p + |u|²/2 Bernoulli head, Γ_n filament
circulations. Clebsch ansatz on Ω_δ: u = ∇φ + α∇β (sign-free: all identities
below are invariant under φ ↦ −φ). Bare action (scope doc): S = ∫dt∫_Ω
[½|u|² + α(∂_tβ + u·∇β) − p(∇·u)]; pressure-eliminated u is the Clebsch
field. No EM symbols, constants, or laws enter any derivation.

## I-CS (kinematic CS-current identities) — ESTABLISHED

1. ω = ∇×u = ∇α×∇β on Ω_δ exactly (curl of gradient vanishes twice).
   [receipt I-CS-a; mutation M1: symmetrizing the ansatz
   u = ∇φ + α∇β + β∇α kills ω identically — the identity reads the
   ansatz's antisymmetric structure, not its shape]
2. ∇·(∇α×∇β) = 0 identically: w ≡ ∇α×∇β is a conserved current as an
   algebraic fact, independent of dynamics. [receipt I-CS-b; mutation M2:
   the dot-structure ∇α·∇β is NOT conserved — the checker distinguishes
   current from non-current]
3. u·ω = det J(φ,α,β) exactly: the α∇β term drops because ∇β·(∇α×∇β) = 0
   (cipher R-a, re-derived here as a full determinant identity, stronger
   than the polynomial-degree check). [receipt I-H-a]

Content: the CS current and the helicity density are KINEMATIC objects of
the Clebsch representation — no equation of motion is used. They hold at
every δ (topological/algebraic form statements; prefire δ-gate satisfied
by construction).

## I-Helicity (exact flux identity) — ESTABLISHED

Statement: for any smooth Euler field on Ω_δ (incompressible, vorticity
form ∂_tu = −ω×u − ∇Π),

  d/dt (u·ω) = −∇·[ Π ω + (ω×u)×u ],

pointwise, as an algebraic identity in the fields. Integrated over Ω_δ with
data compactly supported in Ω_δ (weak-Lin discipline) or u·n = 0 walls, the
right side is a boundary flux; for flow-transported tube walls Ω_δ(t) the
flux through the walls is carried by the core bookkeeping (§03, I-Sing) and
the net is exact: dH/dt = 0 with H = ∫_Ω u·ω.

Derivation (two steps, both receipt-backed):
1. Contract the vorticity-form Euler equation with ω and with u:
   ∂_tu·ω = −(ω×u)·ω − ∇Π·ω = −∇Π·ω, using (ω×u)·ω = 0 [receipt M4a —
   the closure mechanism, asserted directly];
   u·∂_tω = u·∇×∂_tu = −u·∇×(ω×u), using ∇×∇Π = 0.
2. Expand with the two generic vector identities [receipts I-H-b, I-H-c]:
   u·∇×(ω×u) = ∇·((ω×u)×u) + (ω×u)·(∇×u) = ∇·((ω×u)×u) + 0, using
   ∇×u = ω. And ∇Π·ω = ∇·(Πω) − Π∇·ω = ∇·(Πω), using ∇·ω = 0.
   Sum: d/dt(u·ω) = −∇·[Πω + (ω×u)×u]. [full identity receipt: I-Helicity
   — the substitution of the Euler evolution into d/dt(u·ω) is performed
   symbolically and the flux form asserted as an equality]

Mutations: M3 (non-solenoidal surrogate curl u + ∇σ) breaks the pure-flux
form — ∇Π·σ-terms survive: the checker requires ∇·ω = 0, not smoothness.
M4b (density u·v with v independent of ∇×u) breaks the (ω×u)·(∇×u) = 0
step: the closure is special to the helicity density, not to any
divergence-form identity. The checker FAILS on both, as asserted.

Noether reading (framing, literature-standard, not load-bearing here): H is
the charge of particle relabeling (SDiff of Lagrangian labels) for ideal
flow; the direct flux derivation above is the same statement without
variational apparatus, and is what the receipt certifies.

## I-Noether2 (A-sector degeneracy) — ESTABLISHED

Statement: for the gauged kinetic functional T = ½|∇φ + α∇β − καA|²
(cipher bridge Reading 1, sign conventions hers), the variational Hessian
w.r.t. the derivative variables vanishes identically,

  ∂²T / ∂(∂_iA_j) ∂(∂_kA_l) = 0   for ALL i,j,k,l ∈ {x,y,z},

(81/81 entries; receipt I-Noether2), and the Euler–Lagrange source ∂T/∂A_j
is algebraic in the fields — no derivative-of-A term exists to integrate by
parts. Hence the A-equation is order zero: A is auxiliary (a constraint),
and NO propagating gauge dynamics is derivable from the bare action at any
order. This is the Noether-2 form of the bridge's Reading-1 kill: an
independent route (degeneracy of the kinetic form) to the same conclusion
(auxiliary-A), cross-confirming cipher's direct expansion.

Mutation M5: adding (ε/2)|∇×A|² — exactly the kinetic term a propagating
completion must add — produces nonzero Hessian entries (12/81), i.e. the
certificate detects precisely the required extension. The absence of F² is
therefore not an artifact of the expansion but a structural fact about the
bare action's symmetry algebra: the extension must BREAK this degeneracy by
construction, which is the honest content of "new physics required"
(0043-door consistent; SYN missing-construction #2 unchanged).

## I-Sing (core terms through the identity) — ESTABLISHED (topological half cited)

With ω = Σ_n Γ_n t̂_n δ²_δ(x − X_n(s)) + smooth on Ω_δ (G2's printed form;
Călugăreanu–White–Fuller class cited, not recomputed):

  H = H_bulk + Σ_n Γ_n² Slk_n + 2 Σ_{i<j} Lk_ij Γ_i Γ_j.

Kelvin's theorem on each tube (circulation Γ_n around a material loop is
exact for the bare action; derivable from I-Helicity's curl or standard
Kelvin argument) makes the flow an isotopy of the tubes: the filaments are
frozen in. Linking and self-linking numbers are isotopy invariants of
closed curves (standard topology; cited). Therefore

  d/dt Lk_ij = 0,  d/dt Slk_n = 0,  d/dt H = 0  (each an EXACT equality)

for closed transport without reconnection. Reconnection is outside the F2
representation (prefire boundary — the domain rule excludes it); the
statement is δ-independent because it is topological (prefire δ-gate).
[The two algebraic ingredients of the frozen-in argument — I-Helicity's
flux form and ∇·ω = 0 — are receipt-backed; the isotopy-invariance step is
cited mathematics, flagged as such per the firewall input/output discipline.]

## I-Decouple (classification: no charge–motion identity in the bare action) — NOT-DERIVABLE class result

The symmetry group of the bare action is generated by: time/space
translations, rotations, Galilean boosts, the relabeling group SDiff of
Lagrangian labels, and the Clebsch gauge shifts. The corresponding
conserved-current algebra is exactly:

- energy ∫½|u|², momentum (vorticity form ∫u×ω up to flux), angular
  momentum — built monomially from (u, ω, p);
- relabeling Casimirs — built from (α, β, ω) Jacobian-type integrals
  (∫f(α,β)-moments and H = ∫u·ω).

Scale covariance closes the enumeration: the bare Euler equations admit
the scaling u_λ(x,t) = λ u(λx, λ²t) (Π_λ = λ³Π), under which the action
scales homogeneously (S → λ⁻³S) — yielding virial-type identities rather
than a new conserved charge; the virial functionals are built from
(x, u, Π) and carry no Lk-dependence, while H itself is scale-invariant
on this branch (H_λ = λ⁰H). Scale invariance therefore supplies no
charge–motion coupling either, and with it the enumeration is closed.


This is an absence claim INSIDE the bare action only — licensed by the
completeness of the symmetry enumeration for this action (the enumeration
is itself auditable: each generator's charge is printed above), and it
makes no claim inside extensions. Every empirical force/observable hunt of
the missing-5 lane (M1 constitutive tag, M2 scattering law, M3 recoil
ratio, M4 precession scaling) hunted exactly such a coupling INSIDE the
bare action; their measured kills (M1: drift-CONFIRMED indistinct+incoherent;
M3: drift-CONFIRMED kill 09:33Z, M4: VOID confirmed 09:41Z, M2: stays down —
the lane concluded inside the bare action; S1/per-m: drift-CONFIRMED
structural orthogonality) now have a structural
explanation: the hunted coupling is not merely hard to see — it is absent
from the identity algebra. The convergence of measured kills with the
classification is recorded as convergence, not claimed as a new kill.

Extensions reopen the question exactly when they change the symmetry
group/constraint set: adding an (∂A)-kinetic term (M5 detects it), a
Kelvin-gauge coupling (0139's priced construction), or a compactness
source (SYN hole #1). Any such extension inherits I-CS/I-Helicity as
consistency constraints — the extended theory must still carry the
conserved CS current and helicity, now coupled. That inheritance condition
is the paper's hand-off to B4-C1 (BF/CS direction) and to any
Kelvin-gauge construction: it is checkable before building.
