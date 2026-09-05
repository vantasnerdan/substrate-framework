# Complete toroidal inverse-curl transfer, including two discrete carrier jets

This executes the pressure step identified in Section5 of
full-poloidal-operator.md. The stationary same-profile background family
and actual stationary tag continuation are constructed separately in
global-ring-and-tag.md; they are not inferred from this kernel estimate.

## 1. Exact kernel and normalization

Put a fixed smooth meridional source disk D inside |(x,z)|<c and embed
it by X_R(x,z,phi)=((R+x)cos phi,(R+x)sin phi,z). Let E(phi) be the
orthogonal rotation taking the cylindrical frame at phi=0 to that at phi.
For a source w(y)exp(i n phi) in cylindrical vector components, the actual
velocity at phi=0 is the full Biot–Savart integral

    B_R,n w(x)=1/(4pi) integral_D integral_-pi^pi
       [E(phi)w(y) cross (X_R(x,0)-X_R(y,phi))]
       /|X_R(x,0)-X_R(y,phi)|³
       exp(i n phi)(R+y_x) dphi dy.                 (1)

Here x and y denote meridional pairs; the zero entries in X specify the
toroidal angle, not their z coordinate. The sign convention in(1) is
the usual omega(y) cross (x-y). All angular arcs and all exterior
velocity values are retained. The norm on the source/observation disk
uses the physical weight1+x/R; it is uniformly equivalent to dx dz.

Take integer n_R with p_R=n_R/R in a fixed compact interval
0<p_min<=p_R<=p_max. Under s=R phi, the phase is exp(i p_R s).
The local straight limit of(1) is

    B_infty,p w(x)=1/(4pi) integral_D integral_R
       [w(y) cross (x_x-y_x,-s,x_z-y_z)]
       /(|x-y|²+s²)^(3/2) exp(i p s) ds dy,         (2)

which is precisely the complete axial Fourier inverse curl, equivalently
curl_p(-Delta_xy+p²)^(-1), with orientation fixed by the same Cartesian
coordinates. It is not a scalar pressure surrogate.

The exact squared distance is

    |X_R(x,0)-X_R(y,phi)|²
      =|x-y|²+4(R+x_x)(R+y_x)sin²(phi/2).          (3)

This identity supplies the estimates below directly, including near
coincident meridional points. The inverse-curl singularity is locally
integrable; no principal-value cancellation is needed for its L2 Schur
bound (unlike a derivative of that kernel).

## 2. Near and far parts, without truncating pressure

Choose an even smooth cutoff chi(s/L), equal to one on |s|<=L and zero
on |s|>=2L, with 1<<L<<R. This splits an exact integral; the far part
is estimated, not removed from the definition.

For |s|<=2L, (3), E(s/R)=I+O(s/R), and the exact source Jacobian give
uniform Schur bounds for the kernel difference. After integrating over
the bounded meridional source disk, the moments needed for j=0,1,2 obey

    || integral near s^j (kernel_R-kernel_infty) ds ||_Schur
                <= C (1+L²)/R, j<=2.             (4)

One way to expose the singularity is to split |s|<=1 and |s|>=1.
In the first part, coordinate changes are uniformly C1-close and their
kernel difference is bounded by C/R times(|x-y|²+s²)^(-1), integrable
over the two source variables and s. In the second part, the rotation
error is C |s|/R times s^(-2). With the largest weight s² its integral
is C L²/R. Denominator and Jacobian errors are smaller:
C/R+C s²/R² relative to the same kernel. This proves(4) without an
L2 pointwise bound at x=y.

For |s|>=L, integration by parts in the true oscillatory phase is
available because p>=p_min. On |s|<=pi R, (3) is comparable to
|x-y|²+s², with constants independent of R. The smooth far kernel
and its s derivatives satisfy the corresponding inverse-power bounds.
After N integrations by parts, its Schur norm with j<=2 moments is
bounded by

    C_N L^(j-1-N),                                (5)

with harmless smaller inverse-R terms. The original angular integrand
is periodic for integer n, so endpoint terms at phi=+-pi cancel.
Equivalently one integrates on the circle with a smooth central cutoff;
there are no artificial endpoint pressure terms. The straight kernel
has vanishing endpoints at infinity after the same integration by parts.
Taking N>=4 suffices for all three rows.

## 3. Actual discrete carrier derivatives

At fixed R, define forward differences with physical carrier step1/R:

    D_R B_R,n=R(B_R,n+1-B_R,n),
    D_R² B_R,n=R²(B_R,n+2-2B_R,n+1+B_R,n).          (6)

Inside the exact angular integral these multiply the kernel by
R^j(exp(i phi)-1)^j. Those factors are genuinely periodic, so(5)
still applies; treating them as nonperiodic derivatives of a fractional
toroidal harmonic would have introduced erroneous endpoint terms.
For |s|<=2L,

    R^j(exp(i s/R)-1)^j=(i s)^j+O(|s|^(j+1)/R), j=1,2.

The extra near Schur error is again C(1+L²)/R. The far factors and their
derivatives have bounds by powers of s that give(5); their periodicity
retains exact cancellation at the far angular join.

Choose L=R^(1/4), N>=4. Schur's test proves the quantitative complete
local operator convergence

    ||D_R^j B_R,n_R - partial_p^j B_infty,p_R||_L2(D)->L2(D)
                      <= C R^(-1/2), j=0,1,2.     (7)

Here D_R^0 is identity. Differentiating the integral for(2) uses the same
oscillatory far interpretation and the massive p_min>0 domain, so its
parameter derivatives are well defined. A symmetric second difference
can be used as well, with the same bound. Formula(7) controls the actual
integer harmonics; it does not assert a continuous physical angle
quasiperiodicity on a closed ring. Taking p_R->p0 adds the ordinary
continuity error in p.

## 4. Consequence for the full compact-vorticity Euler operator

For a no-swirl ring, the exact compact term in the preceding proof is
the sum of the toroidal derivative of vorticity-weighted velocity and
the velocity multiplied by the background vorticity derivative. In the
physical frame, n/R appears rather than a divergent n: the coordinate
coefficient of the toroidal vorticity is its physical magnitude divided
by r. Formula(7) and C1 convergence of the background vorticity therefore
give operator-norm convergence of this compact term, including its
first two discrete carrier differences. The norm bound contains the
fixed smooth-taper derivatives; the taper is chosen first and not sent
to zero simultaneously with R.

Pullback by a C2-close volume-preserving meridional coordinate map adds
the corresponding small coefficient errors. The actual transport
difference has one angular derivative, but it is small relative to the
reference transport resolvent because partial_theta(z-T0)^(-1) is
bounded by Section2 of full-poloidal-operator.md. The resolvent factor
for L0=T0+K0 has the same domain, so that relative bound transfers to
any fixed Euler resolvent contour.

Thus, GIVEN a global stationary same-profile ring family with these
background/coordinate convergences, (13) in the preceding proof follows
for all sufficiently large finite R. Its positive isolated column mode
continues to a true closed-ring Euler pole. The three discrete kernel
rows also provide the pressure part of the carrier-two continuation;
the stationary family, the eigenvalue/eigenvector normalization and the
tag rows must all be differentiated consistently to complete that join.

This is a derived full-kernel transfer implication, not a numerical
localized residual or an assertion that an arbitrary EPS field supplies
the background family. The exact next construction at this stage is the
global stationary same-profile ring continuation and the material tag
moment map, rather than an unspecified missing pressure estimate. The
third companion executes it without changing the pressure operator.
