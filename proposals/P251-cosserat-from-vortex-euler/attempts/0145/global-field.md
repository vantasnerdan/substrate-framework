# A global constant-curl toroidal field and bounded periodic approximation

## 1. An explicit entire stationary Euler field

Fix kappa,k,lambda>0 with kappa²+k²=lambda² and a real amplitude A.
In cylindrical coordinates define

    psi=A r J1(kappa r) cos(kz),
    u=(-psi_z/r,lambda psi/r,psi_r/r)
      =A(k J1(kappa r)sin(kz),lambda J1(kappa r)cos(kz),
          kappa J0(kappa r)cos(kz)).                       (1)

Full cylindrical differentiation gives div u=0 and curl u=lambda u.
The physical pressure is p=-rho|u|²/2, so the entire field is stationary
Euler. Smoothness on the axis follows from J1(kappa r)/r being an
entire function of r² and J0 being even; the cylindrical basis creates
no coordinate string. The field is bounded with bounded derivatives,
periodic in z and decaying radially like r^(-1/2). It is not of finite
total energy or positive volume-averaged energy density by itself.

At z=0 and kappa R=j_(0,n), a positive zero of J0, the streamfunction
has a nondegenerate extremum. J1 at that zero is nonzero by Bessel ODE
uniqueness. Its Hessian divided by its central value is

    Hess(psi)/psi(R,0)=diag(-kappa²,-k²).

Choosing A so U=A lambda J1(kappa R)>0, the central circle has actual
speed U. The transverse material first jet, in radial/vertical local
coordinates, is

    xdot=(Uk²/lambda) z, zdot=-(U kappa²/lambda)x,
    Omega=U kappa k/lambda.                            (2)

At kappa=k this is circular and Omega=U lambda/2. The axial angular
velocity has radial derivative -2U/R² at the circle; it is not silently
set to zero. Increasing the radial-zero index gives arbitrarily large
R at fixed kappa,k. Normalize A at each such finite R first. Its global
norm grows with R but is finite and independent of later quadrature
accuracy. No limit claiming a uniform nonvanishing density is taken.

On a fixed cross-sectional neighborhood, as R increases the normalized
field tends in every fixed derivative to

    u_infinity=(Uk/lambda cos(kappa x)sin(kz),
                U cos(kappa x)cos(kz),
                -U kappa/lambda sin(kappa x)cos(kz)),   (3)

in the radial/tangential/vertical frame. This follows either from Bessel
large-argument asymptotics or its ODE with initial values J0(R kappa)=0,
J1'(R kappa)=-J1(R kappa)/(R kappa). The local ODE coefficients differ
from the constant-coefficient ones by O(1/R). Curvature of the frame is
also O(1/R), so it is retained in finite-R comparisons.

The resulting invariant closed tori, actual flux-action twist and
arithmetic-free core/boundary choices are proved separately in twist.md.
These are actual unknotted closed tubes, not arbitrary knot types.

## 2. Exact helicity projection of the Bessel integral

Let f=J0(kappa r)cos(kz). The field(1) has the representation

    u=(A/kappa)[curl curl(f e_z)+lambda curl(f e_z)].    (4)

The [Bessel integral identity](https://dlmf.nist.gov/10.9.E1) writes f
as the average of exp(i q.dot(x)) over the two circles
q=(kappa cos theta,kappa sin theta,+/-k). Their positive normalized
measure has total mass one and |q|=lambda. For a plane wave the vector
coefficient in(4) is

    b(q)=(A/kappa)[lambda² e_z-q q_z+i lambda(q cross e_z)].

Direct algebra gives q.dot(b)=0, i q cross b=lambda b, and
b(-q)=conjugate(b(q)). Thus symmetric quadrature gives real exact
Beltrami fields, not approximate solutions of Euler. Moreover for all
directions on the sphere

    |b(q)| <= sqrt(2)|A|lambda²/kappa,
    ||partial^alpha u_quad||_infinity
       <= sqrt(2)|A|lambda^(|alpha|+2)/kappa.           (5)

The bound follows from the total positive quadrature weight, and does
not grow as the number of quadrature points increases. Smoothness of
the integrand and its x derivatives gives C^m convergence on every
fixed compact set. This is an existence/analytic approximation bound,
not a numerical quadrature accuracy claim based on an unmeasured mesh.

## 3. Periodicity without losing the curl eigenvalue or the global bound

Rational points are dense on the unit sphere: stereographic projection
of pairs of rational numbers supplies them explicitly. Replace each
quadrature direction q/lambda by a nearby rational unit direction,
retaining antipodal pairs, and use the SAME helicity coefficient b(q)
for its new direction. The curl identity remains exact, since the norm
is still lambda. Compact C^m errors tend to zero with direction error;
on a ball of radius L their bounds involve at most an extra factor
1+lambda L from differentiating the phase. Bound(5) is unchanged.

A common denominator D of the finite rational directions gives a
common Cartesian period 2pi D/lambda. Consequently the construction
produces exact periodic stationary Euler fields converging to(1) on
any prescribed finite union of compact sets, with global derivative
bounds independent of the requested approximation accuracy. This is
the missing quantitative feature of this route: an unconstrained local
Runge theorem alone does not provide such global bounds.

## 4. Compatibility with robust knotted tubes in the same field

If arbitrary knotted tubes elsewhere in the same global field are
required, take one periodic EPS field w with that robust topology,
using the primary periodic realization licensed in0116. Spatial
rescaling makes its curl eigenvalue lambda; since its original
eigenvalue is an integer, its wave directions are rational on the
unit sphere. Its smooth global C^m norm M_w is finite. Fix the field
w FIRST; then choose a small amplitude epsilon so epsilon M_w is
below the persistence/observation margin of the toroidal core.

Place a marked knot neighborhood of epsilon w at a sufficiently large
radial distance X. The entire field(1) and its fixed finite derivatives
are arbitrarily small there by radial decay; choose them smaller than
epsilon times the topology-persistence margin of w. Thus

    u_total=u+epsilon w(x-X)

is one exact constant-lambda Euler field with both the original closed
core tube and the prescribed distant knotted tubes. It remains bounded.
Approximate u by the periodic construction above simultaneously near
both neighborhoods. The approximating modes and w have a common
period after taking a common denominator, so the resulting sum is
periodic and retains both neighborhoods by their strict margins.

This does NOT put an optical packet at the core circle inside one of
the distant knotted tubes. The distinction remains explicit if that
stronger same-tube identification is used downstream. It supplies
coexistence and exact global stationary compatibility, not a topology
claim about an uncomputed mode support.

Uniform translation of this periodic field in its common cell, together
with one Haar rotation of the WHOLE field and its marks, gives a
stationary isotropic law of actual stationary Euler fields. Its energy
density is finite and positive and its marked tubes have positive
density. This law need not be Gaussian or translation-ergodic; it is
not silently substituted into C-CST-008's Gaussian hypotheses.

## 5. Controlled finite-time transfer of localized actual Euler data

Here is the exact license supplied by the uniform global bound. Suppose
one has already constructed a C² family of smooth localized isovortical
initial data indexed by carrier q in a compact interval, and its actual
linear Euler solutions v(q,t) on a fixed interval[0,T] about the target
field. Let v and its first two carrier derivatives be continuous into
H^(s+1)(R³), s sufficiently large for the material observations. This
is an explicit upstream hypothesis, not a quasimode asserted by(1).

The whole-space pressure projector has L² operator norm one. Uniform
global C^(s+2) background bounds give uniform H^s linear Euler energy
bounds. With L_u v=-P[(u.dot(grad))v+(v.dot(grad))u], the difference
w=v_N-v obeys

    wdot=L_(u_N) w+(L_(u_N)-L_u)v.

For a large fixed ball, the forcing is bounded by the local background
C^(s+1) approximation error times ||v||_(H^(s+1)), plus the bounded
global coefficient norm times the H^(s+1) tail of v. The same estimate
holds for both carrier derivatives, because u_N and u do not depend
on that preparation parameter. Compactness of[0,T] times the carrier
interval makes these tails uniformly small. Choose the ball first,
then the local approximation accuracy. Duhamel gives arbitrarily
small C²_q C_t H^s error with NO approximation-dependent exponential
growth constant. This retains full nonlocal pressure.

Compact Lin generators with a smooth carrier parameter supply smooth
initial velocity families via the fixed full-space projector. The
initial data and complete KKS/Hessian forms converge by the same local
coefficient and L² projector estimates. Actual finite-time transported
tags converge by the local flow/Lin equations on a ball containing
their trajectories; a nonsingular physical-angle chart is retained
with its stated positive margin.

For response errors E(q), a C² bound gives the differentiated-scale
estimate

    |E(q+d)+E(q-d)-2E(q)| <= d² sup |E''|.           (6)

Thus a positive local second-carrier-jet margin can be preserved at
its own natural scale, not merely at the optical zero-order scale.
One fixes the local construction and its finite positive margins before
choosing the periodic approximation. This theorem does not supply those
local optical margins;0142 is still constructing them. Nor does it
identify a globally autonomous Bloch band, replace a physical mean
observation by a tag, or close the full coupled continuum objective.

For genuinely nonlinear Euler histories, choose a small physical
perturbation amplitude e after the fixed finite T and the background
bounds. The relative H^s energy inequality has the form
Y'<=C_M Y+C_s Y². Taking e sufficiently small closes the bootstrap
through T and the nonlinear history differs from its linear tangent
by C_T e². Uniform background bounds make this comparison uniform
over the periodic approximants. This is local-in-time smooth Euler
existence on the specified observation interval; no global regularity
theorem for arbitrary three-dimensional initial data is assumed.
