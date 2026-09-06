# Nonlinear same-leaf range calculation

## 1. Scope and dependency status

This calculation is unconditional up to elementary whole-space Euler/Hodge
identities and the existence of the fixed smooth Cao carrier.  Whenever the
finite bending cluster or a common graph domain from 0052/0054 is mentioned,
it is explicitly an **author-stage conditional input pending 0057 review**.
No 0054 implicit-function, commutator, Riesz, or resolvent statement is used to
cross the nonlinear range.

Write the Lie bracket as

    [X,Y]=(X dot grad)Y-(Y dot grad)X.                    (1)

Let `omega_0=curl v_0` be the centered, axisymmetric, no-swirl Cao ring in a
frame translating with speed `c_0`, and let

    W_0=v_0-c_0 e_z-Omega_0 R,                            (2)

where `R=partial_theta` is the physical rotation generator.  The base ring is
axisymmetric, so inserting the prospective pattern speed `Omega_0` in (2)
does not change its steady equation.

## 2. Exact rotating equation and global Hodge row

For any compact-vorticity, finite-energy candidate set

    W=v-c e_z-Omega R,       omega=curl v.                (3)

The rotating/translating vorticity equation is

    curl(W cross omega)=0
      iff [W,omega]=0.                                    (4)

The sign follows from
`curl(W cross omega)=(omega dot grad)W-(W dot grad)omega`
when both fields are divergence free.  Since `omega` is compactly supported,
`W cross omega` is compactly supported even though `Omega R` is affine.  On
simply connected `R^3`, (4) gives

    W cross omega=grad B,                                 (5)

with `B` constant in the exterior component.  Substituting
`omega=curl v` in (5) recovers

    P_L[(v dot grad)v]-c partial_z v-Omega[R,v]=0.         (6)

Conversely, curling (6) gives (4).  A divergence-free, curl-free, decaying
whole-space remainder is zero, so there is no independent `L^2` harmonic
one-form on `R^3`.  In a core/exterior decomposition, however, the exterior
solid-torus circulation is a genuine Hodge row.  It is fixed by the global
Biot--Savart field and the circulation constraint; dropping it in a core-only
solve would break the equivalence.  Thus (4)--(6) are equivalent only with the
global decay and circulation row retained.

## 3. Physical reverser, period, and KKS sign

In Cartesian coordinates let `S=diag(1,-1,-1)`, which implements
`(r,theta,z)->(r,-theta,-z)`, and define

    mathscr R u(x)=-S u(Sx).                              (7)

For `B(u)=P_L[(u dot grad)u]`, orthogonal equivariance and quadraticity give

    B(mathscr R u)=S B(u)(Sx)=-mathscr R B(u).            (8)

Because `S` reverses `z`,

    partial_z(mathscr R u)=-mathscr R(partial_z u).       (9)

Hence the translating Euler vector field
`X_c(u)=-B(u)+c partial_z u` satisfies

    X_c(mathscr R u)=-mathscr R X_c(u).                  (10)

Also `S R S^{-1}=-R`, so the same identity reverses the rotating parameter.
This verifies the README reverser in the physical translating frame rather
than importing an abstract reversible-Hamiltonian theorem.

For the real cluster basis use the supervisor-corrected conventions

    [R,q_c]=ell q_s,
    X_J=-[R,omega]=-a ell q_s+O(a^2),
    sigma=Omega_KKS(q_c,q_s).                            (11)

With `i_(X_J) Omega_KKS=dJ`,

    J'(a)=Omega_KKS(X_J,q_c)
         =(-a ell) Omega_KKS(q_s,q_c)
         =a ell sigma,
    J(a)=J(0)+(ell sigma/2)a^2+O(a^3).                   (12)

Thus a branch, if constructed, has nonzero action when `sigma!=0`.  Its full
SO(2) group period and its field period under `Z_ell` isotropy are respectively

    T_SO2=2 pi/|Omega|,
    T_min=2 pi/(ell |Omega|),                             (13)

with laboratory displacements `c T_SO2` and `c T_min`.  Equations (11)--(13)
are kinematics; they do not establish a branch.

## 4. An exact same-leaf carrier-map equation

Let `g=Exp xi` be a volume-preserving diffeomorphism equal to the identity
outside the collar and set

    omega_g=g_* omega_0,
    v_g=BiotSavart_R3(omega_g).                           (14)

Naturality of the bracket turns (4) into

    [Y,omega_0]=0,
    Y=g_*^{-1}(v_g-c e_z-Omega R).                        (15)

Therefore the exact nonlinear branch problem can be written as the
**carrier-map/centralizer equation**

    v_g-c e_z-Omega R=g_*Y,
    Y in Cent(omega_0),                                  (16)

together with the exterior trace, circulation, impulse, center, and phase
rows.  This is an exact representation change.  `Y` is a commuting relative
velocity, not a freely changed vorticity profile, while (14) preserves the
entire coadjoint leaf.

The centralizer can be calculated for the actual no-swirl carrier.  In the
positive core write

    omega_0=zeta(r,z) partial_theta,       zeta>0.         (17)

For
`Y=Y^r partial_r+Y^theta partial_theta+Y^z partial_z`,

    [Y,zeta partial_theta]
      =-zeta partial_theta Y^r partial_r
       +(Y^r zeta_r+Y^z zeta_z-zeta partial_theta Y^theta)
          partial_theta
       -zeta partial_theta Y^z partial_z.                 (18)

In a nonzero Fourier harmonic `n`, the coefficient matrix in the ordered
basis `(Y^r,Y^theta,Y^z)` is

              [ -i n zeta       0           0       ]
    C_n =     [ zeta_r       -i n zeta    zeta_z     ],   (19)
              [    0             0        -i n zeta ]

and

    det C_n=i n^3 zeta^3 !=0.                             (20)

Thus `Cent(omega_0)` has **no nonaxisymmetric interior component**.  At
`n=0`, (18) instead requires the poloidal part to be tangent to levels of
`zeta`; locally at regular levels its streamfunction is `F(zeta)`, while the
axisymmetric toroidal component is free subject to the physical domain and
boundary rows.  These are licensed same-leaf streamline variables, but they
are axisymmetric.

Using

    d/dt|_0 (Exp(t xi))_* Z=-[xi,Z],                      (21)

the derivative of (16) at the base is

    M(xi,y,c_1,Omega_1)
      =B_R3(-[xi,omega_0])+[xi,Y_0]
         -y-c_1 e_z-Omega_1 R,                           (22)

where `y in Cent(omega_0)`.  Equations (19)--(20) show that in every
nonzero harmonic

    M_n xi_n=B_R3(-[xi_n,omega_0])+[xi_n,Y_0];            (23)

neither `y`, `c_1`, nor `Omega_1` supplies an arbitrary critical trace in
that harmonic.  The centralizer formulation is still useful: it identifies
the exact global matching equation and the allowed profile variables.  It
also refutes the stronger proposed claim that adding an arbitrary same-leaf
streamline function automatically makes every nonaxisymmetric trace map onto.

## 5. What the physical critical operator actually is

Linearizing (4) in vorticity gives, modulo the finite translation/rotation
rows,

    L eta=-[W_0,eta]-[B_R3 eta,omega_0].                  (24)

For frozen constant coefficients and a divergence-free Fourier amplitude
`k dot eta=0`,

    B_R3 eta=i (k cross eta)/|k|^2,

    sigma(L)(k) eta
      =-i(W_0 dot k) eta
       -(omega_0 dot k)(k cross eta)/|k|^2.               (25)

The second term is the full pressure/Hodge coupling.  On `k^perp`,
`(k cross)^2=-|k|^2`, so its two polarizations have local eigenvalues

    -i D(k) minus-or-plus i (omega_0 dot k)/|k|,
    D(k)=W_0 dot k.                                      (26)

Consequently a zero of the **raw transport** divisor is not by itself the
exact scalar critical equation at finite covector.  In a fixed angular
harmonic, the whole-space Hodge term becomes a nonlocal lower-order radial
operator.  It is compact/order-lowering relative to transport in the usual
core scale, but it remains load-bearing in the critical layer and makes the
operator nonnormal.  Compact perturbation also cannot remove an essential
endpoint merely by moving finitely many eigenvalues.

It follows that the README trace
`T_(m,j)N=N(I_*)` is an exact criterion only **after** a physical
block-conjugation or limiting-absorption reduction has proved that (24) is
equivalent to scalar multiplication by `D_(m,j)` plus a regular divisible
remainder.  Neither the Cao source nor the author-stage 0054 graph/Riesz
calculation supplies that nonlinear, embedded-spectrum reduction.  Using an
IFT on the isolated finite cluster here would bypass the very range theorem
being sought.

For reference, once such a scalar reduction is earned, the elementary local
division lemma is exact.  If `s>3/2`, `D(I_*)=0`,
`|D'(I_*)|>=d>0`, and `N` is one derivative smoother, then

    N=D G in H^(s-1)

implies `N(I_*)=0`; conversely `N(I_*)=0` gives

    N(I)=(I-I_*) integral_0^1 N'(I_*+t(I-I_*)) dt,        (27)

and division by the analogous formula for `D` yields
`||G||_(H^(s-1))<=C_d ||N||_(H^s)`.  This proves the
divided-difference step, not the missing conjugation of (24).

## 6. Structural transparency fails as a universal Euler identity

The quadratic vorticity interaction of two divergence-free Fourier inputs is,
up to the common harmless phase,

    N(k_1,eta_1;k_2,eta_2)
      =k cross (v_1 cross eta_2+v_2 cross eta_1),
    v_a=(k_a cross eta_a)/|k_a|^2,
    k=k_1+k_2.                                           (28)

Take

    k_1=(-1,0,0),       eta_1=(0,-1,0),
    k_2=(0,-1,0),       eta_2=(0,0,-1).                  (29)

Then

    v_1=(0,0,1),        v_2=(1,0,0),
    k=(-1,-1,0),        N=(1,-1,0).                      (30)

Both inputs and the output are divergence free.  For the frozen relative
velocity `W=(1,-1,0)`,

    D=W dot k=0,             N!=0.                       (31)

Hence the Euler quadratic tensor is not polynomially divisible by the output
transport factor.  Incompressibility, curl structure, and the quadratic
bracket alone cannot prove every Cao recursive trace vanishes.  This is a
route-scoped refutation of **universal structural transparency**, not a
calculation of the source-specific Cao numerator and not a no-go for the
branch.

## 7. Dynamical accessibility does not by itself close the range

The principal symbol of a same-leaf vorticity tangent is

    delta omega=curl(xi cross omega_0),

    sigma_DA(k)xi
      =i k cross (xi cross omega_0)
      =i[(k dot omega_0)xi-(k dot xi)omega_0].            (32)

For divergence-free displacement `k dot xi=0` and
`k dot omega_0!=0`, (32) is a nonzero scalar on the transverse plane.  Thus
localized same-leaf directions are symbolically rich enough to prescribe a
transverse vorticity amplitude even when `W_0 dot k=0`.

But the range equation applies (24) to that tangent.  At a critical layer the
transport part loses its inverse and the nonlocal term in (25) remains.  A
symbolically onto DA chart therefore does not provide a bounded right inverse
for `L`, nor does it identify the adjoint trace of the full nonnormal operator.
Together with (20) and (23), this yields the exact route verdict:

- DA tangents are **established symbolically** on the physical coadjoint leaf;
- arbitrary nonaxisymmetric centralizer/profile trace control is **refuted**;
- the source-specific critical trace is **blocked by the missing full-Hodge
  critical-layer reduction and adjoint limiting-absorption calculation**.

## 8. Failure-derived conditional critical-layer candidate

The failed analytic division selects a materially different representation:
solve a resonant layer nonperturbatively rather than demand an analytic
coefficient `N/D`. Suppose first that the missing full-Hodge reduction of
(24) constructs a generalized adjoint resonant polarization, its physical
symplectic measure, and canonical action/phase coordinates. Only under that
hypothesis can one write the generic normal form

    H_rel(I,phi)
      =H_*+(alpha/2)(I-I_*)^2
          +epsilon V_* cos(phi)
          +O((I-I_*)^3+epsilon|I-I_*|+epsilon^2),         (33)

where

    alpha=partial_I D_(m,j)(I_*),                         (34)

and `V_*` would be the generalized-adjoint/limiting-absorption pairing of the
exact whole-space Biot--Savart/carrier-map residual. The conditional leading
normal form is

    dot phi=alpha(I-I_*),
    dot I=epsilon V_* sin(phi),                           (35)

with candidate half-width and small-oscillation frequency

    Delta I=2 sqrt(|epsilon V_*/alpha|),                 (36)
    omega_cl=sqrt(|alpha epsilon V_*|).                  (37)

Equations (33)--(37) do not yet resolve the Euler critical layer. They record
the generic scale that a successful physical reduction must derive, including
its sign and normalization.

A prospective same-leaf realization must also carry the singular scaling. In
local volume coordinates `dI dphi ds`, use a generator of the form

    S_epsilon(I,phi,s)
      =epsilon^a chi((I-I_*)/sqrt(epsilon),s) sin(phi).   (38)

Its `k`th transverse derivative costs `epsilon^(-k/2)` before the amplitude
factor. A collar Bogovskii extension carries corresponding Sobolev losses.
The exponent `a`, support, divergence correction, and the norm of
`g_epsilon-id` must be fixed in the physical tame ledger. If this produces a
smooth volume-preserving flow, then

    omega_epsilon=(g_epsilon)_* omega_0,                 (39)

and the exact residual is

    B_R3((g_epsilon)_*omega_0)-c e_z-Omega R
       -(g_epsilon)_*Y_epsilon.                          (40)

There is a separate topology condition. Equation (16) makes the relative
field smoothly conjugate to `Y_epsilon in Cent(omega_0)`. A smooth
near-identity pushforward cannot create pendulum centers, saddles, or a
separatrix if the chosen axisymmetric centralizer lacks that foliation. A
positive smooth route must construct an admissible `Y_epsilon` with the same
critical-set topology. Otherwise (33) is only a coefficient-space normal
form, or the continuation must justify a singular/weaker-regularity conjugacy.

Thus (38)--(40) define a concrete candidate matching problem rather than an
established Euler inner layer. The remaining theorem must first construct the
full-Hodge limiting-absorption adjoint mode and evaluate `alpha,V_*`; then it
must close the derivative/Bogovskii ledger, centralizer topology, inner/outer
matching, global pressure, and convergence.

## 9. Route and obligation verdicts

1. **Route A, analytic KKS Lyapunov--Schmidt:** blocked at the exact physical
   range step.  Universal quadratic divisibility is refuted by (28)--(31), and
   the exact Cao trace cannot be defined as point evaluation until the full
   nonnormal reduction of (24) is built.
2. **Route B, reversible invariant manifold:** the physical reverser and KKS
   sign are established by (7)--(13), but the coefficient equations inherit
   the same embedded critical range.  Finite-cluster terminology does not
   remove it; the convergent manifold remains unconstructed.
3. **Route C, same-leaf profile/range variable:** the exact carrier-map
   formulation (15)--(16) is established.  Equations (18)--(23) refute direct
   onto-ness of an arbitrary nonaxisymmetric centralizer trace.  The licensed
   axisymmetric centralizer functions remain useful global matching unknowns.
4. **Failure-derived inner critical layer:** the generic pendulum scale and
   exact carrier-map residual are recorded conditionally in (33)--(40). The
   physical adjoint/symplectic reduction, scaled generator and Bogovskii
   ledger, compatible centralizer topology, Cao resonant projection,
   inner/outer tame matching, and convergence remain open.

Accordingly N1 is established only for the unconditional carrier-map,
reverser, action sign, and DA tangent identities; its finite spectral cluster
portion remains conditional on 0057.  N2--N4 and the frozen 0058 success
contract remain open.  There is no exact rotating Cao branch, restoring
neighborhood, particle claim, or quantum inference in this attempt.
