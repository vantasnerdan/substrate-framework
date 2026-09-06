# Global linearized Euler domain and the Gavrilov return wave packet

## 1. Correct global operator, interface, and action sign

Let `u_*` be a stationary divergence-free Euler field and let `P` be the
whole-space Leray projector.  Directly linearising

    partial_t u+u dot grad u+grad p=0,       div u=0

gives the global velocity operator

    partial_t v=G_*v
               =-P[(u_* dot grad)v+(v dot grad)u_*],                 (1)
    div v=0.

This formula contains the nonlocal pressure exactly; there is no local
replacement for `P`.  For the Gavrilov carrier selected below,
`u_* in C_c^infinity(R^3)`.  Let `C_sigma` be the smooth compactly
supported divergence-free core in

    H_sigma={v in L^2(R^3;R^3):div v=0}.                             (2)

The operator `-P(u_* dot grad)` on `C_sigma` is skew.  Denote its natural
skew-adjoint projected-transport closure by `A_0`; `H^1 intersect H_sigma` is
a convenient core, not an asserted description of the maximal domain in the
flat exterior.  The shear `B_*v=-P(v dot grad u_*)` is bounded by
`||grad u_*||_infinity`.  Thus

    G_*=A_0+B_*,       D(G_*)=D(A_0),

and the bounded perturbation theorem gives a strongly continuous group
`S_*(t)` with

    ||S_*(t)||_(H_sigma -> H_sigma)
       <= exp(|t| ||sym grad u_*||_infinity).                         (3)

The sharper coefficient in (3) follows directly from

    (1/2)d/dt ||v||_2^2
      =-int v dot (sym grad u_*)v dx.                                (4)

Thus the full open-space initial-value problem exists in the physical kinetic
norm.  Equation (3) is an existence/continuity bound, not LP2 stability.

Let

    B_orb w=-ad_w^* m_*,       div w=0,
    X_DA=closure_(rho0^(1/2)L2){B_orb w:w smooth and decaying}.       (5)

Equivariance of Euler under volume-preserving relabelling gives
`S_*(t) Ran(B_orb) subset closure Ran(B_orb)`; hence `X_DA` is a closed
invariant subspace.  The physical operator for this attempt is the part of
`G_*` in `X_DA`.  Its graph domain is

    D_DA={v in D(G_*) intersect X_DA:G_*v in X_DA},
    ||v||_Y^2=rho0||v||_2^2+T_0^2 rho0||G_*v||_2^2.                  (6)

This is the graph closure of the smooth orbit-image core.  Euclidean neutral
directions are handled by the quotient distance `inf_z ||v-z||_2`, with `z`
in the finite-dimensional symmetry tangent, rather than by asserting that an
arbitrary orthogonal complement is invariant.  In the axisymmetric packet
below only axial translation occurs; its oscillatory pairing can be removed
by subtracting the exact neutral translation tangent.

In the axisymmetric variables `(eta,chi)` this kinetic norm is exactly the
0019/0014 norm

    rho0||v||_2^2
      =2pi rho0 int[eta K eta+chi^2/r^2]dnu.                         (7)

The corrected orbit convention is `i_X Omega=dH`, `dot q=J grad H`.  On the
stabilizer quotient the quadratic action is therefore

    S^(2)=int[-Omega(q,dot q)/2-Q(q,q)/2]dt.                         (8)

The sign in (8) gives `dot q=J Lq`; the opposite kinetic sign gives
`dot q=-J Lq`.  Since (1) is the direct Euler linearisation, it is also an
independent convention check on (8).

For comparison, the Cao--Zhan generator is not globally smooth.  In its
moving frame let `D={psi_e>0}`, `Gamma=partial D`, with outward normal `n`.
The active generator has

    xi=(psi_e/e) 1_D,
    zeta=[alpha/e^2+psi_e/(e^2r^2)]1_D,
    [zeta]_Gamma=alpha/e^2.                                         (9)

A normal interface displacement `h` belongs to the physical tangent through

    eta=eta_bulk+(alpha/e^2)h delta_Gamma,                           (10)

while `chi` has no surface delta because `xi` is continuous and zero on
`Gamma`.  The complete weak axisymmetric linearisation is

    partial_t chi+u_m dot grad chi+v_m dot grad xi=0,
    partial_t eta+u_m dot grad eta+v_m dot grad zeta
       =(2/r^4)partial_z(xi chi),                                    (11)
    v_m=(-partial_z Keta/r,partial_r Keta/r).

The singular part of (11) is equivalently the kinematic Rankine--Hugoniot row

    (partial_t+u_tau partial_s)h
       =v dot n+h partial_n(u_m dot n)                               (12)

in a signed-distance normal graph.  Thus the Cao domain contains the linked
bulk/interface state `(eta_bulk,chi,h)` with (10)--(12); it is not the raw
interior pair alone.  In Poisson notation this is precisely the missing

    delta[J(q)dA(q)]=J_*d^2A_* delta q+(delta J[delta q])dA_* .       (13)

The second term vanishes on a smooth critical core but not by an interior
identity across `Gamma`.  The primary Gavrilov route avoids this interface:
its cutoff is flat, so `u_*`, its vorticity, and every coefficient in (1)
vanish smoothly through the support boundary.

## 2. A concrete regular annulus in the same Gavrilov family

Use Gavrilov's dimensional radius `R>0`, local coordinates

    x=(r-R)/R,       y=z/R,

and analytic function `a=alpha(x,y)`.  The primary source gives

    a=2(x^2+y^2)+3x^3+3xy^2+O((|x|+|y|)^4),                         (14)

with a strict nondegenerate minimum at `(1,0)` in the original `(r/R,z/R)`
coordinates.  It defines

    p=R^4 a/4,
    b=R^3 sqrt(H(a))/4,
    u=(p_z e_r-p_r e_z+b e_theta)/r,                                (15)

and the source series implies

    H(a)=4a+O(a^2).                                                  (16)

Choose once and for all a smooth `W:R-> [0,1]` such that

    supp W subset (1,2),       W=1 on [5/4,7/4],                    (17)

with flat extension by zero.  For a dimensionless `0<sigma<<1`, set

    omega_sigma(p)=W(a/sigma),
    u_sigma=omega_sigma(p)u.                                        (18)

This is one explicit member of the same exact cutoff family used in 0010.
It is smooth and compactly supported, and its swirl is positive on its
support, so its axial angular momentum remains finite and strictly positive.

By the Morse lemma, for all sufficiently small `sigma` every level
`a=c in [5sigma/4,7sigma/4]` is a smooth simple closed curve around the
minimum.  Hence

    A_sigma={5sigma/4 <= a <= 7sigma/4}                              (19)

is a concrete compact regular invariant annulus.  It lies in the plateau
`omega_sigma=1`; no wall or internal boundary condition is imposed at its
edges.  The exact field outside `A_sigma`, its flat cutoff collars, and the
whole-space Leray kernel remain present in (1).

For any pressure value `p_0` in this annulus, its meridional circuit has the
exact quadrature

    T(p_0)=oint_{p=p_0} r ds/|grad p|,                               (20)

and, using travel angle `vartheta=2pi tau/T(p_0)`,

    Phi_t(p_0,vartheta)
       =(p_0,vartheta+Omega(p_0)t),
    Omega(p_0)=2pi/T(p_0),
    dnu=[T(p_0)/(2pi)]dp_0 dvartheta.                               (21)

The last identity follows from coarea and `dt=r ds/|grad p|`; it fixes the
measure rather than assuming a canonical chart.  The exact cotangent return is

    (k_p,k_vartheta)
      ->(k_p-T(p_0)Omega'(p_0)k_vartheta,k_vartheta).                (22)

Thus a twisted annulus has a skew-product return on the cotangent bundle, not
one finite matrix whose powers remain in a fixed wave-vector fibre.

## 3. Solvable small-shell trajectory and full pressure cocycle

Put `X=r-R`, `Z=z`, and `s=sqrt(X^2+Z^2)`.  Equations (14)--(16), uniformly on
the plateau annulus after the rescaling `s=O(R sqrt(sigma))`, give

    p=(R^2/2)s^2+O(R s^3),
    u_X=R Z+O(s^2),
    u_Z=-R X+O(s^2),
    u_theta=(R/sqrt(2))s+O(s^2).                                   (23)

In a local Cartesian azimuthal coordinate `Y=R theta`, the rescaled leading
flow is therefore the exactly solvable helical shear

    u_0=(Omega Z, gamma s,-Omega X),
    Omega=R,       gamma=R/sqrt(2).                                 (24)

Its trajectories have constant `s`, meridional period

    T_0=2pi/Omega,

and azimuthal advance `gamma s T_0`.  In the moving orthonormal frame
`(n,t,e_Y)` (radial, meridional tangent, azimuthal), the deformation and
cotangent returns over one meridional circuit are

    D Phi_T0=[[1,0,0],[0,1,0],[gamma T_0,0,1]],
    k(T_0)=D Phi_T0^(-T)k(0),                                       (25)

or explicitly

    (k_n,k_t,k_Y)->(k_n-gamma T_0 k_Y,k_t,k_Y).                     (26)

This already separates the two propagation regimes: an axisymmetric ray
`k_Y=0` returns to its covector fibre, whereas a nonaxisymmetric ray is sheared
between fibres and must be iterated as the skew product (26).

Now derive pressure rather than deleting it.  For a high-frequency velocity
packet with wave covector `k` and divergence-free amplitude `A`, the leading
substitution into (1), together with differentiation of `k dot A=0`, gives

    dot k=-(grad u_0)^T k,
    dot A=-(grad u_0)A
       +2 k[k dot (grad u_0)A]/|k|^2,
    k dot A=0.                                                       (27)

The second term in the amplitude row is the order-zero Leray pressure.  For
the returning axisymmetric ray `k=e_t` and initial physical polarization
`A(0)=e_n`, equation (27) reduces exactly to

    dot A_n=0,       dot A_Y=-gamma A_n,       A_t=0.                (28)

Therefore the physical one-circuit amplitude return on `(A_n,A_Y)` is

    M_0=[[1,0],[-gamma T_0,1]]
       =[[1,0],[-sqrt(2)pi,1]].                                     (29)

It has determinant one and preserves the canonical two-form on this
two-dimensional polarization plane.  Independently, the orbit Hamilton
equation uses the corrected KKS action sign (8); no identification of this
auxiliary polarization area with the globally normalized KKS form is assumed.
Its two eigenvalues are one but it is not semisimple.
For the explicit accessible polarization `(1,0)`,

    |M_0(1,0)|^2=1+2pi^2.                                           (30)

This is full-pressure algebraic lift-up, not an eigenvalue or an inference
from 0014's raw Jordan matrix.

The same conclusion is visible in the axisymmetric vorticity variables.  If
`kappa=|k|` and `K#=r^2/kappa^2`, the energy variables are equivalent on the
annulus to `y=(eta/kappa,chi)`.  After removal of the common convective phase,
the leading paired entries are

    y1_dot=i[2xi k_z/(r^4kappa)]y2-(dot kappa/kappa)y1,
    y2_dot=i[r^2b_xi/kappa]y1,
    b_xi=(xi_r k_z-xi_z k_r)/r.                                    (31)

Their product is the finite quantity

    -2xi b_xi k_z/(r^2 kappa^2).                                    (32)

Equations (27)--(29) are the physical velocity/Leray rederivation of this
order-zero coupling for the solvable annulus limit.  They retain wave-vector
transport and show why raw differential order alone was insufficient.

## 4. Bridge from the solvable return to an actual compact Euler wave packet

The analytic series give, after scaling the annulus to a fixed one,

    ||grad u_sigma-grad u_0||_(L-infinity(A_sigma))=O(sqrt(sigma)),
    |T_sigma-T_0|=O(sqrt(sigma)).                                   (33)

Let `C_sigma(t)` be the coefficient on the right of (27), restricted to
`k_sigma(t)^perp` and expressed in the physical parallel frame along the
actual circuit.  The actual one-step polarization cocycle is the path-ordered
exponential

    M_sigma=Pexp[int_0^T_sigma C_sigma(t)dt].

The cylindrical curvature correction and the physical covector drift during
one circuit are `O(sqrt(sigma))`.  The latter is the finite-`sigma` remnant of
(22): the final covector need not equal the initial covector, so `M_sigma` is
a map between the corresponding transverse polarization fibres, not a fixed-
fibre Floquet multiplier.  After parallel identification of those fibres,
continuous dependence of (27) and Gronwall give

    ||M_sigma-M_0||=O(sqrt(sigma)).                                 (34)

No fitted constant is used.  In particular, for all sufficiently small
`sigma`, the continuation of the polarization in (30) has

    |M_sigma A_sigma(0)|^2 > 1+pi^2.                                (35)

To bridge (35) to the full PDE and the physical orbit, choose a smooth tube
inside `A_sigma` around the selected circuit and an integer phase
`N vartheta`.  More explicitly, take a real smooth tube amplitude `A(p,
vartheta)` and the real part of
`a_N=N^-1 A exp(iN vartheta)` in the exact coadjoint generator of 0014, with
`g_N=0`.  Then

    q_N=J_*(0,a_N)=({xi,a_N},0)                                     (36)

is exactly dynamically accessible.  Its Biot--Savart velocity has leading
polarization normal to the meridional circuit and wave covector tangent to it,
which is precisely `(k,A)=(e_t,e_n)` in (28).  The exact Biot--Savart
reconstruction is divergence free and finite energy.  The order-minus-one
symbol gives the desired unit leading velocity amplitude; its longitudinal
and coefficient-gradient corrections are `O(N^-1)`.

For completeness, the fixed-time continuum bridge used here is obtained by
inserting

    exp(iN phi)[A_0+N^-1 A_1]

into (1).  The eikonal equation transports `d phi` by the first row of (27),
the transverse order-zero equation is its second row including the Leray
term, and `A_1` cancels the remaining longitudinal order-zero residual.  Since
the tube stays a positive distance from the axis and support boundary and
`|d phi|` is bounded below, one integration by parts in the whole-space
Fourier multiplier for `P` leaves an `L^2` residual `O(N^-1)`, uniformly on
one circuit.  This also shows explicitly where global pressure, rather than a
local pressure closure, enters the estimate.

Let `v_N(t)=S_sigma(t)v_N(0)` be the exact solution of the global operator
(1), not a truncated ray.  Substitution of the WKB transport (27) produces a
residual `R_N` satisfying on the fixed interval `[0,T_sigma]`

    sup_t ||R_N(t)||_2 <= C_sigma/N,
    sup_t ||v_N(t)-v_N^WKB(t)||_2
       <= (C_sigma/N)T_sigma exp(C_sigma T_sigma).                   (37)

This is the Duhamel estimate from the exact semigroup (3).  Choose `sigma`
first so that (35) holds, then choose the tube and `N` so that the squared-norm
error is less than `pi^2/2`.  The resulting exact, finite-energy,
dynamically-accessible axisymmetric linearized Euler solution satisfies

    rho0||v_N(T_sigma)||_2^2
       > (1+pi^2/2)rho0||v_N(0)||_2^2.                              (38)

The axial-translation tangent is an exact neutral accessible solution.  Its
pairing with this oscillatory packet is `O(N^-m)` for every fixed `m`; subtracting
that component enforces the modulation gauge without changing (38) for large
`N`.  Equation (38) is therefore the promised physical quotient-space growing
wave packet with its continuum bridge.  It occurs on the same smooth compact
finite-`j` Gavrilov carrier and retains the global pressure.  It proves robust
one-circuit transient energy amplification.  It does **not** prove exponential
spectral instability or unbounded all-time growth: (34) controls one return,
while twist and the cotangent shear (22)/(26) generally move later packets
between fibres.

## 5. What the actual return decides and what remains

Three routes now have exact scoped verdicts.

1. **Global operator/domain route -- established for Gavrilov.**  Smooth flat
   cutoff gives the whole-space group (1)--(7), with no hidden interface.  The
   Cao comparison instead requires the linked surface state (9)--(13).
2. **Regular-annulus/physical packet route -- established for sufficiently
   small plateau-shell Gavrilov members.**  The exact leading return is the
   symplectic non-semisimple shear (29), and (33)--(38) transfer its one-cycle
   energy amplification to exact accessible solutions of the full linearized
   PDE.
3. **All-time semigroup/phase-mixing route -- still open.**  Unit eigenvalues
   of (29) are not stability, while one-cycle gain is not exponential
   instability.  The exact cotangent return (22)/(26) shows that generic twist
   invalidates repeated powers of one fixed finite matrix.  The next operator
   must therefore be the skew-product/Koopman return including `K`, not a boxed
   eigenvalue list.

The next executable dependency is to estimate repeated returns of the actual
skew product: conjugate the full `S_sigma(t)` by the characteristic shear,
retain the order-zero pressure cocycle, and prove a polynomial/modulated bound
or construct a sequence of accessible packets with gains persisting over an
unbounded number of circuits.  Nonaxisymmetric `k_Y!=0` begins from the exact
cotangent shear (26) and the full amplitude equation (27); a local multiplier
claim is unavailable until that moving-fibre cocycle is solved.  This is a
specific LP2 route gap, not a global Euler or particle impossibility.
