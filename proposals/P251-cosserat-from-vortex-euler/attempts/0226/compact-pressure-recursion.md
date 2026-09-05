# A compact full-pressure recursion for the inner-core expansion

This is a proved operator construction that can be used in candidate C.
It removes the specific assumption that a small local coefficient error
automatically gives a small global pressure response. Its input is a
finite polynomial inner-core expansion, not the unproved assertion that
all of the curved ring's observed clock rows have already been canceled.

## 1. Compact inverses and their exact moment cost

For m>=0 let

    D_m f=f''+f'/r-m^2 f/r^2,
    M_j(g)=integral_0^infinity r^(m+1+2j)g(r)dr.

For smooth compact regular harmonic profiles, integration by parts gives

    M_0(D_m f)=0,
    M_j(D_m f)=4j(m+j)M_(j-1)(f), j>=1.                 (1)

The regular Newton inverse of g(r)exp(im theta) has a single exterior
harmonic tail, r^-m M_0(g) for m>0, or log(r)M_0(g) for m=0. Thus it is
compact exactly when M_0(g)=0. The exterior constant in the m=0 case is
fixed by the actual Newton convention; for a zero-mass radial source it
is zero. Consequently q successive regular Newton inverses are compact
when M_0(g),...,M_(q-1)(g) all vanish. Equation (1) gives the successive
moment coefficients, rather than a boundary condition imposed on f.

Equivalently, let V_(N,M) consist of compact smooth planar functions
with finitely many angular harmonics |m|<=M and with every polynomial
moment of degree at most N zero. If N>=M, its actual whole-plane Newton
inverse is compact; it belongs to V_(N-2,M). The two-moment loss follows
by solving Delta Q=P for each polynomial P of degree at most N-2:
there is a polynomial Q of degree at most N, and
integral P Delta^-1 g=integral Q g=0. All boundary terms vanish because
the inverse's exterior harmonics were removed by the preceding moments.

Multiplication by a polynomial of degree d costs at most d moment orders
and increases the finite harmonic bound by at most d. Differentiation
does not create an infinite harmonic set and has the corresponding
finite moment cost by integration by parts. Thus a finite composition
of polynomial-coefficient differentiations and Newton inverses has a
finite, explicitly countable initial moment budget. It maps a sufficiently
moment-flat compact input to a compact output. This is an exact
whole-space inverse construction, not a radial wall or a spectral gap.

## 2. Why the free core propagation preserves the construction

At the uniform inner-core limit use u0=(Omega Jx_perp,W0) and the zero
axial carrier. For a solenoidal horizontal field w_perp=J grad psi,
the actual linear Euler solution is rigid transport of vector components
and coordinates:

    w_perp(t,x)=Rot(Omega t)w_perp(0,Rot(-Omega t)x),
    w_z(t,x)=w_z(0,Rot(-Omega t)x).

Its pressure is p0=2Omega psi(t). Indeed the vector-rotation part and
Du0 w supply 2Omega J w_perp=-2Omega grad psi. These signs use the
same J as 0222. Compactness, finite angular harmonics and polynomial
moment flatness are preserved. On each harmonic the time propagator
is a literal exponential exp(-im Omega t).

Now take any specified finite-order expansion of the actual local
Euler coefficients and tube metric in small parameters. Each retained
coefficient is a polynomial in the transverse coordinates; the single
global ring harmonic contributes the explicit n/R factors. Expanding
the pressure/divergence equations to the same order requires only
finitely many operations of the class in section 1. The inhomogeneous
time equations use the free propagator above and finite time integrals.
Their time coefficients are finite exponential-polynomial functions;
time integration does not change the spatial moment budget.

This need not assume differentiability of the global stationary-ring
construction to arbitrarily high order in 1/R. At each selected actual
ring, write its inner coefficient discrepancy as delta_R V_R, with
the proven finite C^j bound on V_R, and Taylor-expand V_R in the inner
transverse variables to the required finite degree. The coefficients
of that polynomial are the actual derivatives of the chosen field.
The Euler operator depends linearly on the velocity discrepancy. A
finite Duhamel recursion therefore has a delta_R^(q+1) remainder, with
the separate transverse Taylor remainder retained, without taking any
derivative of V_R with respect to R. Metric factors have their explicit
geometric series. All constants and preparation moment costs belong
to this specified finite recursion. This is an operator expansion about
an actual field, not an assumed high-order stationary Green expansion.

Choose that budget before selecting the initial compact profile. The
complete finite recursion then has compact pressures and velocity
coefficients supported in the selected inner preparation region. It
does not excite an unrepresented exterior critical layer at its
retained orders. It makes no statement about a Riesz pole of the full
all-poloidal Euler operator.

The initial radial return and the finite polynomial moments are
compatible with a nonzero tag response. Their kernels are polynomial
radial weights after integration by parts. Choose the tag derivative
weight to be a smooth nonanalytic compact bump on the observed annulus.
It is not in their finite polynomial span. Hence the response functional
does not vanish on their common kernel. The additional frequency-flat
tag rows of density-normalization.md remain linearly independent: on
an open interval where the tag weight is nonzero they reduce to
independent powers of the strictly monotone radial clock, while an
off-tag open interval tests the polynomial return rows. The positive
and negative quadratic returns can still be chosen on disjoint
frequency subintervals in this finite-codimension kernel.

This compatibility statement is for polynomial exterior-pressure
moments and the listed radial tag rows. An arbitrary new curved
observation row is not thereby proved independent of the desired
angle. That distinct question is retained below.

## 3. Full nonlocal Euler error, without an inverse-gap shortcut

The first actual toroidal pressure correction illustrates the recursion.
For the ring coordinates of 0222 and global harmonic n,

    Delta_ring=Delta_perp+(R+x)^-1 partial_x
                              -n^2/(R+x)^2.

The first coefficient of a Newton solve is therefore
Delta_perp p1=-partial_x p0, in addition to the explicitly retained
first source coefficient. If p0=A(r)exp(i theta), its derivative has
only harmonics 0 and 2. The two source profiles are

    g0=-(A'+A/r)/2,   g2=-(A'-A/r)/2.

Their exterior moments satisfy

    integral r g0=0,
    integral r^3 g2=2 integral r^2 A.                    (3)

Thus the same compact m=1 return moment that appears in 0222 removes
both tails of this first curvature correction. For the actual Kelvin
pressure A includes its actual W(s) multiplier; expanding that
multiplier requires the further polynomial moments, not deletion of W.
Higher source and metric coefficients are handled by the explicit
finite budget of section 1. Equation (3) is an executed first step on
the actual toroidal Laplacian, rather than an assumed local inverse.

Suppose the recursion is carried to order q and produces compact fields
w_app,p_app, and an exactly solenoidal compact velocity after its finite
divergence completion. On its support assume the actual background
coefficient remainder has the stated C^(s+1) bound. Compute, in physical
coordinates,

    r=partial_t w_app+u.grad w_app+Du w_app+grad p_app.

If ||r||_(C^r_t H^s)<=epsilon_q and the actual initial velocity differs
by epsilon_q in the same norms, the exact whole-space Euler solution
obeys on a fixed interval

    ||w-w_app||_(C^r_t H^s)<=C_(T,s,r) epsilon_q.          (2)

Projecting the displayed equation gives the actual residual P r, whose
norm is at most ||r||. The standard differentiated Euler energy estimate
then proves (2), with coefficients controlled by the actual global
derivative bounds. A large constant far velocity is transport and does
not create an exponential norm cost. There is no replacement of P by a
local projection and no use of a box Poincare constant. Lin transport
gives the corresponding displacement estimate, including its initial
data and its explicitly computed residual. The same estimate then
controls actual tag, G, spin and shape rows with their own norms.

For the curved geometry one should use the exact Piola/curl completion,
as in 0222, rather than claim that a truncated solenoidal series is
exactly divergence-free. Its finite correction and its residual belong
in epsilon_q. Finite parameter derivatives can be included by
differentiating this same residual before the energy estimate.

## 4. Earned repair and the next actual coefficient calculation

Sections 1-3 prove a full-pressure compact parametrix construction for
any finite polynomial inner-core expansion with its stated remainder.
They are materially stronger than multiplying a local error by an
uncontrolled global inverse norm. Together with the exact nonlinear
root of density-normalization.md, they provide two independent tools:
compact pressure returns and positive quadratic-action centering.

They do not prove that the actual curved ring's complete registered
angle history is one exponential to that order. The recursion can
produce additional exponential-polynomial tag rows. Their independence,
their possible cancellation by actual preparations, and the associated
quadratic forms must be computed at the first order where they occur.
Likewise an actual stationary array and its acoustic centroid/ambient
response are not supplied by the nominal density accounting. Those are
the next constructive obligations; no parent completion is inferred.
