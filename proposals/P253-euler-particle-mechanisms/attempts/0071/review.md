# P253/0071 independent review of final 0068

## Frozen transaction

This is one fixed-boundary independent review of root-owned attempt `0068`.
The reviewer authored or implemented none of the target action, derivation,
API, tests, verifier, precision repair, or receipts. The corrected frozen
README has SHA-256
`70dd4e2324173548bdca2b73a3f1b61b6c24998cc9b794c476615139c4d646a6`.
Central activation exited exactly `0`; its command, stdout, empty stderr, and
exit hashes are
`0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
Every target/API/test hash matched the preactivation freeze.

Attempt `0070`, `euler_gauge_emergence.py`, and its tests were not opened or
used. No production numerics or unchanged target test was run. One bounded
0068 domain/evidence-role correction was requested and passed its
correction-only check. Its receipt has SHA-256
`265714b42208c3d2532ec4a916387e2da71cfd2c00089e012acc818b188dfdf5`.

## Source applicability

Accepted 0042/0045 supplies the local transported material tag and the exact
causal retained-history construction for a single Euler velocity field. It
does not supply a gauge field, Maxwell constraint, Lorentz force, mixed-state
projection, charge normalization, or a charged stationary carrier. Final 0068
correctly rederives the enlarged current, action, PDE, and constraint rows
rather than inheriting them.

Maxwell/U(1) is explicitly introduced as a comparator foundation with positive
inputs `epsilon`, `mu` and free coupling `g`. It is not an accepted Euler
claim. The displayed action is therefore the load-bearing source for the
Maxwell signs, material backreaction, energy exchange, and Noether identities.
Attempts 0064--0067 are route context only and do not become authority here.

The pre-review precision record correctly separates free-field and material
coupling terms, Maxwell-subsystem speed from the elliptic incompressible
pressure sector, fixed tag leaf from varying physical charge, static response
from a carrier branch, and translation symmetry from the missing cokernel
calculation. Its one remaining domain defect is the untyped use of
mixed-state `L2` orthogonality, addressed below.

## Joint action, current, and constraints

For a transported signed tag,

    chi_t+u dot grad chi=0,      div u=0,
    rho_q=g chi,                 J=g chi u,

direct expansion gives

    partial_t rho_q+div J=0.

Thus total charge `Q=g integral chi` is conserved when the spatial boundary
term vanishes. Under
`A -> A+grad lambda`, `phi -> phi-partial_t lambda`, the coupling variation is
the spacetime boundary term minus
`integral lambda(partial_t rho_q+div J)`, so gauge invariance holds for the
declared compact or vanishing-boundary gauge functions.

With

    E=-grad phi-partial_t A,     B=curl A,
    L_field=integral [epsilon|E|^2/2-|B|^2/(2mu)],
    L_coup=integral [J dot A-rho_q phi],

independent potential variations give

    epsilon div E=rho_q,
    (1/mu)curl B-epsilon E_t=J,
    div B=0,                     B_t=-curl E.

The signs agree with the positive electric energy and with current work.
Variation of the volume-preserving material map, with the coupling included
exactly once, gives

    rho_m(u_t+u dot grad u)
      =-grad p+rho_q(E+u cross B).

This is a coherent charged incompressible-fluid extension, not bare Euler.
Dotting the material equation with `u` and combining it with the actual
Poynting identity cancels `integral J dot E`. Under the stated integrability
and stress-flux decay, the resulting energy and translation momentum are

    integral [rho_m|u|^2/2+epsilon|E|^2/2+|B|^2/(2mu)],
    integral [rho_m u+epsilon E cross B].

The field momentum coefficient and density signs are correct. These balances
require solved smooth equations and convergent boundary fluxes; they are not
created by the symbolic oracle.

The Eulerian vector field is

    u_t=-P_L(u dot grad u)+(g/rho_m)P_L[chi(E+u cross B)],
    chi_t=-u dot grad chi,
    E_t=(epsilon mu)^(-1)curl B-(g/epsilon)chi u,
    B_t=-curl E.

Its linear constraint subspace is

    K_g={div u=0, div B=0, epsilon div E-g chi=0}.

Taking the corresponding divergences and using continuity proves that the
vector field is tangent to `K_g`. The nonzero-charge whole-space field has a
`1/r^2` electric tail, which is finite-energy and compatible with field
`H^s`; its scalar potential behaves as `1/r` and need not belong to
inhomogeneous `L2`. Periodic net charge would instead violate the Gauss zero
mode unless compensated. Final 0068 consistently uses the whole-space sector
and compact-support/finite-first-moment hypotheses for its Coulomb asymptotic.

## Local PDE and retained history

For `s>5/2`, `H^s(R^3)` is an algebra and gives Lipschitz velocity. The Euler
transport estimate, Maxwell symmetric-hyperbolic energy, Hodge-bounded Lorentz
products, and tag transport close on a common interval for compatible
`H^s` data in `K_g`, with the additional integrability/decay used by global
observables. This is a standard local energy-method theorem, not a Banach ODE
and not an all-time result.

Let `Pi` be a bounded finite-rank idempotent with smooth range, bounded on the
mixed `H^s` state and smoothing from `H^(s-1)` to `H^s`, such that `Pi` and
`Q=I-Pi` preserve `K_g`. For a prescribed resolved path

    v in C(H^s) intersection C^1(H^(s-1)) intersection Ran(Pi) intersection K_g

and compatible
`w_0 in Ran(Q) intersection K_g intersection H^s`, the equation

    w_t=Q F(v+w),       w(0)=w_0

has the corresponding local energy-method solution. Requiring
`v_t=Pi F(v+w)` reconstructs the full system exactly. The initial unresolved
fluid, tag, electric, and magnetic components remain load bearing. This is the
correct enlarged causal-history theorem.

Raw `L2` orthogonality is neither used nor physically typed on
`(u,chi,E,B)`: the components have different dimensions and `chi` has no
canonical Maxwell/fluid energy coefficient. A valid orthogonality claim would
need an explicit nondimensional weighted direct-product Hilbert inner product.
The exact history result needs only the bounded reducing projection above.
The requested bounded correction therefore preserves the theorem while
removing an unnecessary domain ambiguity.

**Constraint/local-history verdict:** established after the bounded
projection-language repair. The local theorem and history map are analytic;
the recorded algebra checks are only corroborating evidence.

## Static Coulomb calculation and stationary-carrier boundary

For a supplied static smooth compact charge profile,

    -epsilon Delta phi=rho_q,
    phi=G*rho_q/epsilon,       G(x)=1/(4pi|x|).

Integration by parts gives the positive electric energy

    (epsilon/2)integral|grad phi|^2
      =(1/(2epsilon)) double_integral rho_q(x)G(x-y)rho_q(y).

The cross term contains no extra factor `1/2`. For separated compact profiles
with finite first moments and charges `q_1,q_2`, expansion of `G` gives

    E_12=q_1 q_2/(4pi epsilon d)+O(d^-2).

The sign and translation derivative yield like-sign repulsion and
opposite-sign attraction. Bounded smooth compact profiles have finite
self-energy because the Newton singularity is locally integrable.

This exact static theorem does not solve the material momentum equation. A
stationary charged carrier would also need pressure/Lorentz balance, a
stationary transported tag, the Maxwell constraints, localization, and the
carrier complement. Final 0068 states that distinction explicitly.

**Static-response verdict:** established for supplied profiles. Actual charged
stationary-carrier existence remains blocked by the stated coupled branch
construction.

## Small-coupling continuation, translation, and stabilizer

The proposed route fixes a material tag leaf and
`C_chi=integral chi_0`, while

    Q_g=g C_chi.

This is the correct small-coupling limit; a single nonzero physical `Q` is not
held fixed while `g` tends to zero. Maxwell fields sourced by the tag/current
are `O(g)`, their energy and Lorentz backreaction are `O(g^2)`, and the
candidate Euler correction is correspondingly `O(g^2)`. In a frame moving at
`c_0`, the principal Maxwell operator

    -Delta+(c_0^2/c_EM^2)partial_z^2

is elliptic for `|c_0|<c_EM`. These are exact conditional order and symbol
statements, not an implicit-function theorem.

The missing branch theorem must quotient translations, fix center/impulse and
neutral rows, identify the stationary nonlinear map and its spaces, invert the
actual KKS/Hessian complement, and prove the translation adjoint pairing.
Translation invariance of the action gives total Noether momentum for solved
evolution; it does not remove the stationary cokernel. The target correctly
keeps that row open.

A locked profile `chi_0=F(I)` is stationary only when
`W_0 dot grad chi_0=0`. Persistence on the joint orbit further requires

    [xi,omega_0]=0  implies  xi dot grad chi_0=0.

This condition says every allowed vorticity stabilizer also preserves the tag.
It is not automatic from passive transport or from an axisymmetric
potential-vorticity analogy. A joint orbit chart and relative
energy--Casimir/electric `H^-1` coercivity remain necessary for persistent
charge localization. These are honestly named missing constructions.

**Branch verdict:** blocked, while preserving an executable conditional route
with exact `Q_g`, field, force, and comoving-symbol scaling.

## Parameter-normalization invariants

Under the constant potential redefinition
`(A,phi) -> a(A,phi)`, the same classical equations are parametrized by

    epsilon -> epsilon/a^2,
    mu -> mu a^2,
    g -> g/a.

Therefore `epsilon mu` and `g^2/epsilon` are invariant, whereas the individual
coefficients depend on the chosen field normalization. Independent rescaling
of the classical tag can be compensated in `g`; fixing `C_chi` (or setting it
to one) removes that redundancy along the proposed branch. Multiplying the
whole action by a constant leaves classical equations unchanged while scaling
symplectic/action periods. Consequently
`c_EM=1/sqrt(epsilon mu)` does not select charge magnitude or quantum action.
No smallest charge, spin magnetic moment, gyromagnetic factor, quantum rule,
or electron parameter is derived.

**Normalization verdict:** established as an exact nonselection ledger, not
as physical parameter selection.

## Oracle audit and finding

The importable API exactly evaluates
`(rho_q,J)=(g chi,g chi u)`, the Maxwell speed, and the leading Coulomb
energy/force, with useful dimension and nonzero-separation checks. Its three
focused tests are valid algebraic regression at that scope.

The standalone oracle's continuity and gauge rows encode the already-derived
identities; its Coulomb derivative and speed checks evaluate real formulas.
However `PASS 5` checks only `work+(-work)=0`, and `PASS 6` only
`-div_current+div_current=0`; they do not independently derive Poynting's
identity or constraint propagation. This does not weaken the analytic action
and divergence proofs above, but the validation evidence role must be stated
truthfully.

| Finding | Direct evidence | Minimum correction | Upgrade path |
|---|---|---|---|
| Mixed-state raw `L2` orthogonality is physically untyped and unused. | `(u,chi,E,B)` has dimensionally unlike components and no declared tag energy weight. | Declare a nondimensional weighted product inner product, or drop `orthogonal` and retain bounded finite-rank idempotent `Pi`, bounded `Q`, smooth range/smoothing, and reduction of `K_g`. | Supply the explicit nondimensionalization and prove orthogonality/boundedness in that Hilbert product if an energy-orthogonal split is later needed. |
| The validation overcredits tautological oracle rows as derivations. | `PASS 5` and `PASS 6` contain cancelling placeholder symbols rather than the action/PDE expressions. | Relabel them as consistency/regression checks and credit the analytic derivation, or replace them with derived symbolic identities. | Encode the full Poynting divergence and differentiated Gauss/magnetic constraints with meaningful sign mutations. |

The bounded correction package fully closes these connected domain and
evidence-role findings. No equation, API result, local-history theorem, or
charged-carrier boundary was weakened.

## Four-axis decision and frontier

- **Verification:** the joint action, local constraint propagation, static
  Coulomb law, conditional branch scaling, and normalization ledger are
  analytically established; API/oracle evidence is regression at its actual
  scope.
- **Review:** established as an explicit constrained U(1) foundation extension
  after one bounded domain/evidence-role correction; charged stationary carrier
  remains blocked.
- **Compatibility:** internally compatible with positive Maxwell energy,
  incompressible local Sobolev evolution, and transported tags. It is a changed
  foundation, not an accepted bare-Euler consequence.
- **Epistemic:** active comparator/conditional extension; no LP4, P4, P5,
  electron, particle, or campaign conclusion.

The strongest supported result is an exact signed material current, joint
Euler--Lorentz/Maxwell action and energy exchange, static Coulomb sign, local
constraint-preserving evolution, and exact retained-state/history
representation on the extended state. The decisive next construction is the
small-`g` stationary continuation on an actual persistent carrier, including
the carrier complement, translation cokernel, stabilizer-preserving tag, and
joint localization coercivity. Parameter selection and P4/P5 remain separate.

## Correction check

The single bounded correction is complete and passed. Receipt SHA-256 is
`265714b42208c3d2532ec4a916387e2da71cfd2c00089e012acc818b188dfdf5`.
It pins these corrected artifacts:

| Artifact | Corrected SHA-256 |
|---|---|
| `0068/derivation.md` | `41bfff94e0c72d0c028564b2ee7c04d231fc3ea3eca75d1d7c556facb86d7816` |
| `0068/validation.md` | `15a83b8888a4737264d3f877ea0030fcbf52e400ce79692c8c587ffcd28c6a58` |
| `0068/sol-high-pre-review-precision.md` | `c41c7236be4f8c8db14655805b28dc353139437139d8f05d41c1356aee2161a6` |
| `0068/artifact-hashes.sha256` | `0e86b493f0361f6806e187b1c8c02ec226a4ffc66ee47ea5a775c48c8f3c7683` |

The check was limited to the projection-domain language, oracle evidence-role
statements, and their direct result edges. `Pi` is now a bounded finite-rank
idempotent with smooth range; both `Pi` and `Q=I-Pi` are bounded on the typed
product `H^s` and `H^(s-1)` spaces and reduce `K_g`; only `Pi` carries the
finite-range smoothing assumption. No mixed-state orthogonality is asserted.
The validation assigns the Poynting and constraint-propagation proofs to the
displayed action/PDE calculation and accurately labels the placeholder rows as
regression. Maxwell-only propagation wording is explicit, while the full
incompressible no-cone statement remains.

The affected stale-language scan, target YAML parse, terminal-newline checks,
and scoped `git diff --check` pass. The API, verifier, and tests were unchanged
and not rerun. No second substantive pass was performed, and no additional
correction is required.
