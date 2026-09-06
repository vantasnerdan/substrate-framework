# Exact transported-current U(1) extension

## 1. The reviewed material tag gives a conserved current

Let `chi` be a smooth signed, compactly supported or integrable material tag
and let `div u=0`. The reviewed 0042 transport equation is

    partial_t chi+u dot grad chi=0.

For a constant coupling `g`, define

    rho_q=g chi,              J=g chi u.                 (1)

Then

    partial_t rho_q+div J
      =g(partial_t chi+u dot grad chi+chi div u)=0,       (2)

and `Q=integral rho_q dx` is conserved under the stated decay. A signed
initial tag gives exact opposite charges. The normalization `g` is a new
constant; material transport selects neither its magnitude nor a smallest
nonzero tag integral. Its individual value depends on the chosen tag and
gauge-potential normalization.

With `E=-grad phi-partial_t A`, `B=curl A`, separate the free field
Lagrangian from its coupling:

    L_field=integral [epsilon |E|^2/2-|B|^2/(2 mu)] dx,
    L_gauge=L_field+integral [J dot A-rho_q phi] dx.     (3)

For smooth `lambda` compactly supported in spacetime, or with vanishing
endpoint and spatial boundary terms, under
`A -> A+grad lambda`, `phi -> phi-partial_t lambda`, the coupling
changes, up to spacetime boundary terms, by

    -integral lambda(partial_t rho_q+div J) dx dt=0.      (4)

Thus the transported tag supplies the precise current identity needed for
gauge invariance; it does not supply the gauge field or coefficient values.

Variation of `L_gauge` gives

    div(epsilon E)=rho_q,
    (1/mu)curl B-epsilon partial_t E=J,
    div B=0,             partial_t B=-curl E.            (5)

The homogeneous principal speed is

    c_EM=1/sqrt(epsilon mu).                             (6)

## 2. Gauss constraint gives the electric Coulomb sign

For a static charge with `A=B=0`, (5) gives

    -epsilon Delta phi=rho_q,
    phi=(1/epsilon)G*rho_q,   G(x)=1/(4*pi|x|).          (7)

The field energy is

    E_E=(epsilon/2)integral |grad phi|^2
       =(1/(2 epsilon)) integral rho_q(x)G(x-y)rho_q(y) dxdy.  (8)

Two separated compactly supported smooth charge form factors with finite first
moments and total charges `q_1,q_2` have

    E_12=q_1 q_2/(4*pi*epsilon*d)+O(d^-2),               (9)

so like signs repel and opposite signs attract under `-grad_d E_12`. Smooth
compact bounded charge density has finite self energy because the Newton
kernel is locally integrable. This correct sign comes from the constrained
Gauss variable and positive electric energy, unlike eliminating a healthy
propagating scalar with a linear source.

**Route A verdict:** established exactly for the declared U(1) extension.
It supplies signed conserved charge, gauge invariance, Coulomb sign and finite
propagation for the Maxwell subsystem, with `g,epsilon,mu` explicit new inputs.

## 3. One material-map action gives backreaction

Let `eta(a,t)` be volume preserving, `chi_0(a)` the signed material tag, and
`chi(eta(a,t),t)=chi_0(a)`. The joined action is

    A_joint=integral dt [
       (rho_m/2) integral |partial_t eta|^2 da
       +L_field
       +g integral chi_0(a){partial_t eta dot A(eta)-phi(eta)} da]. (10)

The last integral is precisely the Eulerian coupling in (3). Divergence-free
variations of `eta` and arbitrary gauge variations give

    rho_m(partial_t u+u dot grad u)
       =-grad p+rho_q(E+u cross B),                      (11)

together with (1)--(5). Thus the same coupling produces source, Lorentz force
and work. It is a modified incompressible charged-fluid theory, not unforced
Euler.

For smooth decaying solutions,

    E_tot=integral [rho_m |u|^2/2+epsilon |E|^2/2
                         +|B|^2/(2 mu)] dx              (12)

is conserved: the fluid work is `integral J dot E` and Poynting's identity
gives its negative. Translation invariance similarly conserves

    P_tot=integral [rho_m u+epsilon E cross B] dx        (13)

when these integrals and the Maxwell stress boundary flux converge.

**Route B verdict:** established at the exact smooth variational and balance
scope. It changes the substrate equations and introduces the gauge state.

## 4. Local state and retained history

In Eulerian form the coupled evolution is

    u_t=-P_L(u dot grad u)+(g/rho_m)P_L[chi(E+u cross B)],
    chi_t=-u dot grad chi,
    E_t=(epsilon mu)^-1 curl B-(g/epsilon)chi u,
    B_t=-curl E.                                        (14)

For `s>5/2`, compatible initial data in `H^s` with the finite-energy/decay
conditions, `div u=div B=0`, and `div(epsilon E)=g chi` have the standard
quasilinear transport plus symmetric-hyperbolic energy structure. The
commutator energy estimate closes locally because `H^s` is an algebra and
`grad u` is bounded. The constraints propagate by (2) and (14). This gives a
common smooth local interval depending on the `H^s` norms; it is not an
all-time theorem.

Let

    K={div u=0, div B=0, div(epsilon E)-g chi=0}.        (14a)

For any bounded finite-rank idempotent `Pi` that reduces `K`, with smooth
range, with both `Pi` and `Q=I-Pi` bounded on the declared product `H^s` and
`H^(s-1)` spaces, and with the finite-rank smoothing bound needed from
`H^(s-1)` to `H^s`, write the product state
`Y=(u,chi,E,B)=v+w`, `w=QY`. Given
`v in [C(H^s) intersection C^1(H^(s-1))] intersection Ran(Pi) intersection K`
and `w_0 in Ran(Q) intersection K intersection H^s`, the unresolved equation

    w_t=Q F(v+w),       w(0)=w_0                         (15)

has the corresponding local energy-method solution on a common interval, and
substitution into `v_t=Pi F(v+w)` is an exact causal history map. The initial
unresolved fluid, tag and gauge fields remain load bearing. Equation (15) is
the enlarged analogue of 0042; the old velocity-only response cannot simply
be reused.

The static Green calculation in section 2 is an interaction law for supplied
charge form factors, not a stationary charged-carrier theorem. For a generic
transported `chi`, the self-force `rho_q E` need not be a gradient and its
deviatoric part can deform or disperse the Euler carrier. A particle
application needs a separate small-`g` continuation with `chi` locked to a
carrier material invariant, the tag leaf fixed while `Q_g=g integral chi`
varies along the branch, translations quotiented, and the `O(g^2)` Maxwell
self-field controlled in the carrier KKS/Hessian complement. Admissible
variations at each fixed `g` stay on that branch member's fixed `Q_g` level.

**Route C verdict:** established at local smooth retained-state scope for the
declared extended equations, not inherited as a bare-Euler theorem. The
Maxwell subsystem has characteristic speed (6); the coupled incompressible
system still contains the elliptic pressure projection and therefore does not
have a strict finite propagation cone for every component.

## 5. Particle ledger

The extension supplies exact classical signed charge, current, Coulomb and
magnetic Lorentz coupling, and a finite field speed. It does not supply:

- a mechanism selecting `g`, `epsilon`, `mu`, or a smallest charge;
- an intrinsic spin magnetic moment or a gyromagnetic relation;
- the action quantum, Hilbert probabilities, exchange sign, or reset;
- a persistent carrier or its scale; or
- a weak chiral current, neutrino mixing, or oscillations.

Under `(A,phi)->a(A,phi)`, the same theory is parametrized by
`epsilon->epsilon/a^2`, `mu->mu*a^2`, and `g->g/a`; hence `epsilon*mu` and
`g^2/epsilon` are invariant while the three displayed coefficients separately
require a field/tag normalization convention. Multiplying a classical action
by an overall constant does not change (5),
(11), or (14) while it rescales symplectic/action periods. Independently,
`epsilon*mu` fixes the wave speed. The gauge speed therefore does not select
the quantum action normalization. Route D is established as an exact
input/output ledger and refutes any claim that this minimal extension alone is
an electron.

The positive result is reviewable if the owner later chooses a foundation
extension: it is the smallest current/action system found in this campaign,
not a global minimality theorem,
that gives the electric Coulomb sign. Until that choice, it remains a
comparator and cannot earn LP3, LP4 or P5 for the frozen Euler objective.

## 6. Failure-derived charged-carrier continuation

The arbitrary-tag defect identifies a concrete joined route with 0066.  In the
translating frame of a persistent carrier put

    W_0=u_0-c_0 e_z,
    chi_0=F(I),                  W_0 dot grad chi_0=0,    (16)

where `I` is an actual smooth material/streamline label and `F` is supported
strictly inside the regular core. Freeze this tag leaf and its integral
`C_chi=integral chi_0` while varying `g`; then the physical charge level is
`Q_g=g C_chi`. Equivalently normalize `C_chi=1` and use the small physical
charge `Q_g=g` as branch parameter. Seek

    u_g=u_0+O(g^2),       E_g,B_g=O(g),                  (17)

with admissible variations restricted to the `Q_g` level at each fixed `g`,
not with one nonzero `Q` held constant across the limit `g->0`. In a standard
gauge the comoving Maxwell principal operator is

    -Delta+(c_0^2/c_EM^2) partial_z^2,                  (18)

which is elliptic exactly when `|c_0|<c_EM`. Solve (18) first with the smooth
current `g chi_g u_g`. Its Lorentz self-force enters the augmented Euler
equation at order `g^2`.

The remaining construction is a same-leaf Lyapunov--Schmidt/IFT at the
corresponding fixed `Q_g` level:
quotient translations, fix center/impulse and every carrier neutral row, and
invert the actual carrier KKS/Hessian complement. Translation invariance and
zero integrated Maxwell self-force make the Noether row a candidate
solvability identity, but the actual cokernel pairing on the chosen KKS/graph
domain must still be derived. A positive 0066 coercive or
isolated-mode complement would therefore turn the 0068 comparator into an
actual persistent charged carrier. At present that complement is not earned,
so (16)--(18) are an executable conditional route rather than a branch
theorem. Charge/action quantization and magnetic moment remain separate.

Existence alone would not prove persistence of the charge profile. At `g=0`,
`chi` is passive and has rearrangement/filamentation directions invisible to
the Euler kinetic energy. The joint orbit must therefore satisfy

    [xi,omega_0]=0  implies  xi dot grad chi_0=0.         (19)

A relation `chi_0=F(I)` supplies (19) only after `I` is proved to be the same
global carrier invariant preserved by every allowed stabilizer; identifying
`I` with potential vorticity is restricted to the axisymmetric no-swirl
setting and is not a full three-dimensional constitutive law. After the orbit
chart, a joint relative energy--Casimir/electric `H^(-1)` coercivity estimate
must control the topology in which charge localization is claimed. Until that
estimate is proved, all-time charged-profile persistence remains open even if
the small-`g` branch exists.
