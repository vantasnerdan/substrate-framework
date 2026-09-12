# Neutral chiral-current and flavor-propagation suppliers

## 1. Passive material tags do not oscillate

Let `u` be a smooth incompressible Euler velocity and let
`chi=(chi_1,...,chi_N)` satisfy

    D_t chi := partial_t chi + u dot grad chi = 0.          (1)

For a material trajectory `Xdot=u(t,X)`, the chain rule gives

    d/dt chi(t,X(t))=0.                                    (2)

Consequently every time-independent local observable `F(chi)` and every
constant mixture `U chi` is constant on that trajectory. In particular, a
finite multiplet of passive tags cannot autonomously generate flavor
oscillation.

If one declares a space-time dependent unitary frame `y=U(t,x)chi`, then

    D_t y=(D_t U)U^* y=-i H_U y,
    H_U=i(D_t U)U^*.                                      (3)

Equation (3) is an identity for a chosen frame. Along any single material
curve its connection is pure gauge. It becomes physical only if an
independently derived interaction fixes the flavor frame at preparation and
readout and the same action supplies nontrivial curvature or holonomy. A
hand-chosen `U` merely rewrites the constant vector in (2).

This proves Route N0 is refuted as a passive-tag oscillation mechanism. It
does not address genuine Euler modes, retained nonlocal variables, or a new
internal connection.

## 2. Exact vector and axial continuity currents

For any transported scalar `s`, incompressibility gives

    partial_t s + div(s u)=0.                              (4)

Thus `(s,su)` is a Galilean continuity-current tuple. When `s` is a true
scalar, its spatial current is polar. When `s` is a pseudoscalar, `su` is
axial. Neither tuple is a Lorentz four-current before the shared P4
representation is constructed.

Euler also has a source-derived axial current. Put

    omega=curl u,   h=u dot omega,
    B=p+|u|^2/2.                                         (5)

The equations imply

    partial_t u=u cross omega-grad B,
    partial_t omega=curl(u cross omega).                  (6)

Using `div omega=0` and
`div(A cross u)=u dot curl A-A dot omega`, one obtains

    partial_t h
      =-div(B omega)+div((u cross omega) cross u)
      =-div[h u+(p-|u|^2/2)omega].                        (7)

Hence

    j_H=(h, h u+(p-|u|^2/2)omega)                         (8)

is a conserved axial Galilean current on any domain where the displayed
flux is meaningful. Under spatial reflection `h` is a pseudoscalar and the
flux in (8) is axial. It is not a passive convective tag current: it contains
the nonlocal Euler pressure and vorticity flux.

A formal combination `j_V-kappa_W j_H` is conserved whenever both terms are,
but Euler does not fix the dimensionful coefficient `kappa_W`, identify the
combination with a one-handed Lorentz spinor current, or link its coupling to
the electron current. Full `O(3)` covariance also forbids identifying a
nonzero scalar density with a pseudoscalar density without a parity-breaking
structure. A helical carrier can select geometric handedness, but does not
by itself supply a fermionic chiral projector. Route N1 therefore establishes
the exact axial supplier (8) and refutes only the zeroth-order passive
convective scalar/pseudoscalar-tag route to a derived chiral current. A
complete first-derivative invariant-tensor and continuity classification,
including gradient, cross-product, and superpotential currents, remains
blocked rather than refuted.

## 3. Exact two-state propagation algebra

Let an actual normalized two-mode state at fixed momentum `p` evolve under a
Hermitian propagation matrix

    H_prop(p)=U(theta) diag(nu_1(p),nu_2(p)) U(theta)^*,   (9)

where `U(theta)` is the real two-dimensional rotation. If the interaction
basis state is initially `(1,0)`, direct exponentiation gives

    P_(1->2)(p,t)
       =sin^2(2 theta) sin^2(Delta nu(p)t/2),
    Delta nu=nu_2-nu_1.                                  (10)

For a normalized momentum distribution `w(p)>=0`, the incoherent detector
average is

    Pbar_(1->2)(t)
       =sin^2(2 theta)/2 [1-Re C(t)],
    C(t)=integral w(p) exp[-i Delta nu(p)t] dp.            (11)

The group-velocity difference is

    Delta v_g(p)=grad_p Delta nu(p).                      (12)

For a packet of width `Delta p`, coherence therefore requires control of
`t Delta p |Delta v_g|` and of higher dispersion derivatives. Equation (10)
alone contains no localization or coherence theorem.

If the shared P4 action is `S_*` and the physical branches have the
relativistic form

    nu_j(p)=sqrt(p^2 c_*^2+m_j^2 c_*^4)/S_*,             (13)

then at high momentum

    Delta nu(p)=Delta(m^2)c_*^3/(2p S_*)+O(m^4/p^3),     (14)

and for travel time `t=L/c_*` and `E approximately p c_*`,

    Delta Phi
      =Delta(m^2)c_*^3 L/(2 E S_*)+higher order.          (15)

Thus the familiar `Delta m^2 L/(2E)` structure is a consequence only after
the same physical action and relativistic dispersion are derived. Equations
(9)--(15) are an exact conditional propagation ledger. They do not derive
`theta`, the masses, `S_*`, the propagation modes, or the current that fixes
the interaction basis.

## 4. Test of the current Cao suppliers

The independently reviewed Cao spectrum supplies distinct simple
positive-Krein modes at fixed carrier parameters. Distinct radial/Sturm
labels and rational-ray crossings establish real Euler mode structure, but
they do not yet provide:

1. a neutral persistent full-three-dimensional P2 carrier neighborhood;
2. a shared P4 spin-one-half/action representation;
3. momentum-dependent freely propagating branches `nu_j(p)`;
4. an interaction current that fixes a flavor basis and is linked to the
   electron current; or
5. a same-action noncommuting propagation/current pair with derived mixing.

The fixed-frequency Schwinger--Hopf algebra can reproduce (10) after choosing
two axes, but prepared beating is not the mass/mixing mechanism required by
P6. Route N2 is blocked at the five named same-carrier constructions, rather
than refuted.

## 5. Connection and extension routes

A material-frame connection generated only by `U(t,x)` in (3) is a gauge
choice. A positive Route N3 construction must give a base space, an
Euler-derived connection with nontrivial curvature or defect holonomy, a
finite-energy core, a retained-state action, and an interaction that fixes
the endpoint frame. Classical continuous holonomy does not itself supply
spin one-half, fermionic exchange, or a universal action. No such compact
same-carrier connection is constructed here, so N3 is blocked at that object.

A chiral `SU(2)_L` connection with a chosen left-handed representation can
make a coherent weak-current and mixing framework, but its gauge action,
coupling, representation, anomaly/current conditions, mass data, and action
normalization are added foundation inputs. The existing Euler plus transported
`U(1)` action does not derive them. Route N4 is blocked at construction of the
joint chiral action, representation, and same-carrier map; its explicit input
ledger is the route's positive proposal result.

## 6. Strongest continuation

The exact positive atom is the conserved Euler axial current (8); the exact
failure is passive-tag oscillation (2). The next constructive supplier must
join (rather than name separately) one persistent neutral P2 carrier, the
shared P4 action/spin state, a momentum-dependent propagation matrix, and an
electron-linked interaction current. Its first exposing test is whether the
two matrices fail to commute in a covariantly fixed physical basis while the
packet satisfies (11)--(12). The current 0088 response theorem may become an
input only after independent review and only at its accepted scope.

None of these results establishes a neutrino, P6, or a no-go for nonlocal,
topological, mode-based, or separately approved foundational routes.
