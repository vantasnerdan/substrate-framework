# P253/0100 — Ertel axial current and same-carrier chiral-lock test

Owner: `root`

Status: **README-only preregistration awaiting central schema activation.**
No source body, derivation, API, verifier, comparator, or numerical work may
open until repository validation is replayed with exit exactly `0`.

## Frozen objective and authority boundary

This attempt continues the P6 current supplier without claiming P6. It asks a
fixed exact question: can one transported material scalar and the Euler
vorticity supply, from the same carrier, a conserved axial current that is
dynamically related to the transported polar current?

Corrected 0068 and independent 0071 supply only the transported scalar polar
current. The exact axisymmetric charged Cao branch is consumed only through
final 0080/0084 at its transported-tag scope. P253/0097 is author-stage
motivation, not accepted authority; every identity used here is rederived.
Active 0095 and 0096, the 0094 action route, and every empirical neutrino
parameter are excluded.

## Route A — exact Ertel material current

For incompressible Euler and a transported true scalar `chi`, define

    q_E=omega dot grad chi.                                (1)

Derive directly

    D_t omega=(omega dot grad)u,
    D_t grad chi=-(grad u)^T grad chi,
    D_t q_E=0.                                             (2)

Thus `(q_E,q_E u)` is a conserved Galilean axial continuity-current tuple.
Freeze its `O(3)` parity, regularity, support/decay, and boundary flux domain.
Compare it with the polar tuple `(chi,chi u)` from the same material tag.

The global charge test is load-bearing. Since `div omega=0`,

    q_E=omega dot grad chi=div(chi omega).                 (2a)

With this parity choice, `q_E` is a pseudoscalar and `q_E u` is an axial
spatial current; a pseudoscalar `chi` reverses that ledger. Therefore
`int_R3 q_E dx=0` whenever `chi omega` is smooth and compactly
supported, or decays strongly enough for the flux at infinity to vanish. On a
bounded domain the integral is exactly the boundary flux of `chi omega`. A
nonzero integrated Ertel charge consequently requires boundary flux, a
puncture or distributional defect, or a non-global/multivalued label whose
topology obstructs the global exact-divergence representation. Nontrivial
topology by itself does not evade Stokes when `chi omega` is globally smooth
and single-valued on a boundaryless domain. The local conservation law alone
does not supply an independent nonzero global charge.

Route A earns an established verdict only for (1)--(2). It does not create a
Lorentz chiral current, weak coupling, spinor projector, or neutrino.

## Route B — constitutive vector/axial lock

Test the strongest local lock

    q_E=lambda chi,                                        (3)

with a typed coefficient `lambda`. Because both sides are transported, (3)
is preserved if imposed initially. The attempt must determine whether Euler
selects it dynamically or whether it is merely an initial constitutive row.

If `lambda` varies, preservation of (3) additionally requires
`D_t lambda=0` wherever `chi` is nonzero. On a closed vorticity line
`dx/dtau=omega(x)`, equation (3) becomes

    d chi/dtau=lambda chi.                                 (4)

For real constant nonzero `lambda`, periodicity forces `chi=0` on that closed
line. For variable real `lambda`, the exact monodromy condition is
`exp(int_loop lambda d tau)=1`, equivalently `int_loop lambda d tau=0`; it does
not force `lambda` to vanish pointwise. For complex locks retain the full
exponential condition with the chosen single-valuedness convention. This is a
route-scoped obstruction for closed vorticity-line carriers, not for open
lines, defects, complex phases, nonlocal locks, or added connections.

## Route C — exact Cao same-carrier test

For the axisymmetric no-swirl Cao carrier, vorticity is toroidal and the
reviewed material charge tag is `chi=F(I)` with poloidal gradient. Therefore

    omega dot grad chi=0.                                  (5)

Audit the cylindrical normalization and the tag-support band, then decide the
Ertel axial-current route on the actual charged carrier. A vanishing (5)
refutes this current supplier only for that same-carrier tag geometry. It
activates a nonaxisymmetric tagged carrier, a multi-label current, or a
nonlocal/defect current rather than narrowing P6.

## Route D — nonaxisymmetric and multi-label continuation

Construct or identify an actual persistent carrier on which `q_E` is nonzero
and the polar/axial current pair is tied to one interaction action. Competing
subroutes are:

1. a nonaxisymmetric Euler carrier with open or nontrivially linked vorticity
   lines and a smooth transported tag;
2. a boundary-flux or puncture/distributional-defect realization, or a
   non-global/multivalued label whose topology obstructs a global exact
   divergence, supporting nonzero integrated Ertel flux; this includes a
   two-label advected flux two-form only when its globality, compactness, and
   defect coefficients are derived; and
3. an internal material connection whose holonomy replaces the prohibited
   real exponential lock on closed lines.

The positive target requires finite energy, an open nonlinear persistence
neighborhood, an electron-linked interaction current, and the shared P4
spin/action representation. A prescribed tag profile or helicity sign alone
does not satisfy it.

## Verification boundary

The exact oracle will derive the cancellation in (2), the divergence identity
and global-flux consequence in (2a), the constant and variable closed-line
monodromy in (4), and the toroidal-poloidal orthogonality in (5). Importable
APIs may expose those identities and their parity/domain ledgers. They cannot
prove a persistent nonaxisymmetric carrier or a physical chiral weak current.

No production numerics are planned. Each route receives one verdict and no
route verdict propagates to LP6, a neutrino, or the parent campaign.
