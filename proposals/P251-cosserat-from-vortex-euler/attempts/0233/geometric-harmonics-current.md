# Actual first geometric harmonics and the polar-current return

The axial carrier coefficient is computed in the companion body. This
body keeps the separate geometric coefficient kappa=1/R and the actual
closed-ring observation. The two are combined only with k=n/R, n=-1.

## 1. Actual centered inner coefficient

Write phi0=delta/2+A_seed J0(lambda r). Expansion of the actual
0211 equation -Delta phi+kappa phi_x/(1+kappa x)=f(phi), in its inner
affine region, gives

    (Delta+lambda^2)phi1=partial_x phi0.

The full Green expansion in 0186 (equation (4) of
global-ring-and-tag.md), with the 0211 source factor, has only the
first transverse harmonic at this order. The logarithmic linear term
is the translation/cokernel term already removed by the actual U border.
The invertible radial border and source-center condition exclude an
independent first-order radial amplitude. Recenter on the actual
meridional critical point; the remaining regular m=1 homogeneous
solution is fixed by grad phi1(0)=0. The resulting local coefficient is

    phi1=-A_seed x J2(lambda r)/2
            =-Omega0 x r^2/8+O(Omega0 lambda^2 r^5).      (1)

The quoted Green remainder and the smooth bordered map give the
corresponding O(R^-2 log^2 R) remainder on fixed core norms. This uses
the full ring construction, not arbitrary boundary data for the inner
Helmholtz equation.

In the right-handed coordinates of 0222,

    u1,perp=-Jgrad phi1-x u0,perp,
    W1=lambda phi1-x W.

At leading inner-core order these are

    u1,perp=(3Omega0 xy/4, Omega0(y^2-5x^2)/8),
    W1=-Wc x+O(lambda Omega0 r^3),
    div_perp u1,perp=-u0,x.                             (2)

Thus both the true first field coefficient and its volume divergence
are retained.

## 2. Full pressure and all generated first harmonics

In orthonormal ring coordinates the first geometric Euler equations
for the horizontal velocity v1 have the actual form

    (partial_t+T0)v1+Du0,perp v1=F1-grad p1,
    div_perp v1=-w0,x,
    F1=2W b0 e_x-u1,perp.grad w0,perp-Du1,perp w0,perp. (3)

The axial Euler equation retains grad W1.w0,perp and the connection
W w0,x. Their leading terms cancel because grad W1=-Wc e_x. Hence
there is no spurious large axial forcing left by deleting one of them.
The remaining axial terms, including W'v1, are controlled at their
actual small inner-core order.

The initial horizontal Kelvin force after its exact gradient is removed
is

    Finit,1=lambda[(h-lambda S)grad phi1
             +x(2lambda S-h)grad phi0+2WS e_x].           (4)

Use its curl for the initial zeta1. Every term in (3)-(4) has only
local scalar harmonics 0 and -2, starting from m=-1. The scalar
vorticity and divergence solution is explicit:

    (partial_t+Omega0 partial_theta)zeta1
                           =curl F1-2Omega0 div v1,
    v1=Jgrad Delta^-1 zeta1+grad Delta^-1(div v1).         (5)

For each generated harmonic j=0,-2 the time coefficients are finite
linear combinations of

    e^(-ij Omega0 t),
    I_j^ell(t)=integral_0^t e^[-ij Omega0(t-s)]
                                     s^ell e^(iOmega0 s)ds.

The integral is explicitly the ell-th derivative, divided by i^ell,
with respect to eta of
[e^(i eta t)-e^(-ij Omega0 t)]/[i(eta+jOmega0)] at eta=Omega0.
The denominator is plus or minus Omega0, not an assumed pressure gap.
Lin reconstruction adds its explicit time integrals, including resonant
t e^(-ij Omega0 t) terms when they occur. These formulas specify every
unobserved harmonic as well as the observed rows.

For example the original h=lambda S supplier has, at leading order,

    zeta1,0(0)=2iOmega0(A'+A/r),
    zeta1,-2(0)=-2iOmega0(A'-A/r),
    zeta1,0(t)=zeta1,0(0)e^(iOmega0 t),
    zeta1,-2(t)=zeta1,-2(0)[2e^(2iOmega0 t)-e^(iOmega0 t)]. (6)

There really is an off-resonant 2Omega0 harmonic. Its absence from one
angle row is angular selection, not absence from Euler dynamics.

All Newton inverses in (5) are the full planar coefficients of the
physical ring pressure expansion. Compactness follows from their
computed exterior moments as in 0226, not a wall. The initial potentials
f,H allow polynomial-moment returns on an off-tag interval. At this
first order the polynomial multipliers have degree at most four along
a composed term, there are at most two inverse Laplacians, and the
Cartesian harmonic bound is at most four. Vanishing polynomial moments
through degree sixteen is therefore a sufficient explicit common
budget. Extra unused moments add no scientific assertion. A maximal
independent basis of those actual polynomial rows remains independent
on an off-tag open interval by analytic continuation, so their targets
can be solved there without changing the physical rows or the fixed
inner core jet.

## 3. Physical selection and the actual centroid coefficient

The leading global angle and G test the local m=-1 material density.
The generated scalar 0,-2 harmonics in (5), the first tag correction
chi_phi phi1, and the first Jacobian/Euclidean-test corrections have
zero first geometric contribution to those leading rows by their
explicit angular products. They do not all disappear from the centroid.

For the general corrected lift of the companion body, put
F=lambda f-g. At leading inner-core order the m=0 toroidal Lin row is

    Xi_zeta,1,0(t)=C(r)[e^(iOmega0 t)-1]
                   -iWc lambda t(A_F'+A_F/r)e^(iOmega0 t),
    A_F=F,
    C(r)=(r^2 h'+3rh)/16
          -(Wc/Omega0)(A'+A/r)
          +(Wc lambda/Omega0)(F'+F/r).                   (7)

This follows from the actual leading axial Lin forcing

    -u1,perp.grad h-2Wc Xi0,x+u0,x h.

The m=0 coefficient of its h part is
iOmega0(r^2h'+3rh)/16. Its remaining coefficients are obtained from
Xi0,perp=Jgrad S+i lambda Omega0 t Jgrad F. The small terms already
retained in the full equations supply the stated higher inner-core
remainder; they are not asserted to vanish exactly.

The measured linear-time centroid coefficient is proportional to
integral r chi(F'+F/r)=-B(F), and vanishes by the explicitly solved
physical matrix. The other coefficient is

    integral r chi C
       =-B(r^2 h)/16+(Wc/Omega0)B(A)
                             -(Wc lambda/Omega0)B(F).    (8)

It is generally nonzero. Thus the raw actual centroid has a constant
plus an oscillatory correction, which cannot be hidden in a spin name.

There is an explicit compact repair. Add the initial toroidal component
kappa C(r)exp(-i varphi). Its scalar m=0 moment satisfies

    integral r C(r)dr=0,

because every term is a compact radial boundary derivative. Its regular
m=0 inverse potential is therefore compact. Complete the generator
exactly with -ik kappa grad Delta_0^-1 C and the Piola map. The added
initial scalar is constant under the leading free transport; it removes
the constant part in (7). The surviving polar coefficient (8) multiplying
exp(iOmega0 t) is retained as physical data.

The divergence completion is of order kappa^2, but its density moment
contributes at the first relative centroid order because of the real
lever arm R. Dropping it would incorrectly call the toroidal field a
tag-changing solenoidal generator. Its first angle/G/spin cross rows
vanish by the displayed local harmonic selection, and its first phase
and energy cross rows with m=-1 vanish as well. Its small actual Euler
velocity is obtained from the full Kelvin projection, not supplied
as a force.

The exact Euclidean observation formulas remain

    Q=delta integral rho chi(x+i y)z,
    theta=Q/(i D_tag),
    G=integral rho chi x cross Xi,
    S_mech=G_t+2rho integral chi Xi cross u,
    M_tag delta X=integral rho chi Xi.

For the global n=-1 preparation, only the transverse centroid/spin and
the xz,yz symmetric-shape components survive the global azimuthal
integrals. The latter are Q itself. In particular, to first order

    G_+=-iQ+2 integral rho chi y_local h e^(i varphi),     (9)

with the actual volume measure. The extra term is retained in G/theta
and therefore in its matched inertia. The physical matrix in the
companion body sets S_mech-G_t to zero at the retained order; it does
not delete the nonzero term in (9).

## 4. Scope of the actual first matrix

The exact source equations, compact inverses and explicit moment minors
give a controlled first-curvature/first-carrier physical-history
construction on the SAME Euler ring. The whole-space residual is
estimated before applying the actual Euler/Lin energy estimate. The
small inner-core remainder, the O(kappa^2 log^2 R) geometric remainder,
the finite-return conditioning and the weighted physical moment norms
are separate terms. Finite time and a chosen finite derivative order
are fixed before making them smaller than the actual phase/current
margins.

No exact all-time scalar oscillator at finite lambda a, arbitrary-order
row compatibility, common-K continuum, or stationary array is inferred.
The next consumer receives the actual first matrix, its positive
normalization construction, and its retained polar/shape rows.
