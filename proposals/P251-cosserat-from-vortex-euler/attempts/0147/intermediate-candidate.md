# Intermediate candidate record — superseded by the executed repair

This is the pre-verification construction record, NOT the final theorem.
Its proposed estimate(13) is insufficient: two carrier derivatives cost
(pL)², so first slow-time matching leaves delta²(pL)² rather than delta².
The TG anisotropic fast harmonics also cannot be bounded by declaring
every time derivative slow. The executed repair in finite-packet.md
uses n=8 and seven slow-time moments on Lundquist, and toroidal-transfer.md
constructs high-jet positive multi-CK geometry instead of silently
dropping the fixed TG anisotropy. The positive full-packet KKS and
distributed-marker calculations here are retained as derivation history.

Central preregistration was validated263/12 before0145 was read. This
file repairs the Fourier-fiber normalization at the new0147 boundary.
It does not edit0142 or identify a registered marker with an absolute
director. The second file supplies the actual toroidal transfer.

## 1. A direct material-generator construction on the circular cell

Use oriented normal coordinates (x1,x2)=(z,x) and longitudinal s=y
for the Taylor--Green limit of0145, kappa=lambda/sqrt(2):

    u1=-(U/sqrt(2))sin(kappa x2)cos(kappa x1),
    u2= (U/sqrt(2))cos(kappa x2)sin(kappa x1),
    W=U cos(kappa x1)cos(kappa x2), Omega=lambda U/2.

This is an exact entire constant-lambda stationary Euler field,
independent of s, with bounded derivatives independent of packet
parameters. Its linear material/Kelvin generator is exactly

    A xi=P(xi cross omega0)-[u,xi], div xi=0.                  (1)

Choose xi initially, set v=P(xi cross omega0), and evolve xi_t=Axi.
Then v is an actual linear Euler solution and xi its actual Lin
displacement. This avoids assuming an invertible global periodic
Lin resolvent. The full-space/fiber P is retained.

For the negative longitudinal carrier -p, put delta=sqrt(lambda/p)
and ell=(2/(lambda p³))^(1/4). The normal first jet is Omega Jx;
the axial quadratic jet is U-U lambda²r²/4, exactly as in0142.
At quartic order the axial difference from Lundquist is

    W_TG-W_L=-U lambda^4 r^4 cos(4theta)/192+O(U lambda^6 r^6).
                                                                  (2)

The cubic normal-velocity difference acts on the scaled packet at
order Omega delta³; (2) contributes to pW at order Omega delta^4.
These are derived background Taylor differences, not a claimed
constant-lambda approximation to an arbitrary variable-factor core.

The all-finite-order0142 construction extends directly to (1).
In the laboratory longitudinal-following frame, the uniform-rotation
principal eigenvalues on helicity/azimuthal components are
(m-2)Omega and (m+2)Omega. The C4 sector containing m=2 also contains
m=-2; both are principal resonances and BOTH are retained. The first
pressure/axial-trap correction on these two resonant sectors is the
signed radial oscillator: its coefficients are respectively
+sqrt(2)(2n+2) and -sqrt(2)(2n+2). The selected n=2 positive value
6sqrt(2) is simple and separated from every other radial value by
at least2sqrt(2). All other C4 angular sectors have a nonzero integer
multiple-of4Omega principal gap. Thus no angular degeneracy is
silently discarded when (2) mixes harmonics.

Incompressibility fixes the longitudinal component via the nonzero
carrier. The full Leray inverse is (p²-Delta_perp)^-1. Its exact
resolvent expansion and weighted remainder are0142(2). At every
finite order, polynomial-Gaussian right-hand sides are solved by
the nonresonant helicity/angular gaps and the signed oscillator
inverse on the complement of the selected mode. Projection onto
that mode fixes the frequency coefficient. Its leading nonzero KKS,
together with the exact Hamiltonian identity for(1), makes this
coefficient real; equivalently one pairs the order equation with
the mode in the Hermitian KKS form, where all previously fixed
coefficients and the Hessian are real/Hermitian. The normalized
coefficient is unique. This constructs a smooth solenoidal xi_N and
real laboratory-frame frequency to arbitrary prescribed finite order.

Define v_N=P(xi_N cross omega0) EXACTLY. Then

    ||A xi_N+i omega_N xi_N|| <= C_N Omega delta^N ||xi_N||,    (3)

including any fixed weighted Sobolev and carrier-derivative norms,
after taking enough recursion orders. There is no residual Kelvin
leaf error. Exact Euler evolution is used after preparation. The
global bounded TG coefficients give a carrier-independent energy
bound for(1); derivatives in the preparation carrier have only the
finite losses already budgeted in0142. No isolated eigenvalue or
nonlinear periodic solution is inferred.

The corrections from the extra resonant angular component start at
order delta², since its forcing is order delta³ and its selected
signed-oscillator gap is order delta. They cannot change the leading
order-delta physical curvature or the leading marker-rank matrices.
They are retained in the exact observation/reconstruction operators.

## 2. A finite-action packet, not an assigned axial length

Let the reference carrier p_* be fixed and positive. A smooth spectral
cutoff is one on |q-p_*|<p_*/8 and zero for |q-p_*|>p_*/4. Write it
as chi_band(q). It is fixed when the carrier p is varied. With L>0,

    a_L(q-p)=L exp[-L²(q-p)²/2] chi_band(q),
    Xi_p(x,s)=(2pi)^-1/2 integral a_L(q-p)xi_N(q,x)e^(-iqs)dq.
                                                                  (4)

The same formula with exact evolved xi(q,t) gives the actual packet
history. It has finite full-space Sobolev norm, is divergence free,
and is C² in p. It is also Schwartz in s and normal position at t=0.
No arbitrary Lz remains. Whole-space pressure is exactly the direct
integral of the fiber pressures, so no end-pressure term is omitted.
Let Xi_1=Re Xi_p, Xi_2=Im Xi_p. Plancherel gives its complete KKS

    beta_packet(p)=integral |a_L(q-p)|² beta_fiber(q)dq>0,       (5)

where beta_fiber is the KKS per unit longitudinal length in the
same complex normalization. For a constant profile this is exactly
sqrt(pi)L beta_fiber. Smooth q dependence gives the full integral,
not an assigned multiple of the center value. Two real Euler
columns preserve (5) at every time. Spectral-cutoff errors and
their carrier derivatives are exponentially small in (p_*L)².

The field is a packet of actual Euler solutions. Frequency variation
across its narrow band is kept in its integral; no common exact
frequency is asserted. Choose p first, then L with pL sufficiently
large (pL>=delta^-3 is ample below), then all geometry/transfer
parameters. Finite optical time |Omega t|<=T0 is fixed.

## 3. Distributed helical material tag and its physical meaning

Fix a small constant c>0 and let ell_z=cL. At reference time define
a nonnegative material fraction

    w0=mu chi(r)chi_z(s0)[1+eps b(r,s0)cos(2theta-p_*s0)].       (6)

Here chi(r) is the fixed radial marker of0142, chi_z is a Gaussian
of width ell_z with a fixed smooth remote cutoff, and b is bounded
on its support. Every choice is frozen at p_*. Take mu<=1/2 and
|eps b|<1. The material label s0 labels actual transported sheets;
it is not reset to Eulerian position during evolution.

On each sheet define its physical central quadrupole in the normal
plane. Its initial reference phase is p_*s0. Rotate by this KNOWN
reference registration and compare with the same sheet transported
by the unperturbed field. The collective angle is the chi_z-weighted
average of these actual sheet-angle variations. It is measurable
from the transported labeled material. It is not the unregistered
quadrupole of the whole long packet and not an absolute core director.

The physical spin is the actual total tagged angular momentum about
its physical centroid, projected along the fixed longitudinal axis.
It includes position variation and the entire transported tag. The
transverse C2 symmetry makes the planar centroid correction vanish
at the linear TG comparison; the exact toroidal transfer retains
the full centroid subtraction when that symmetry is perturbed.
The factor cos(2theta-p_*s0) in(6) makes the total axial spin coherent
with (4), instead of cancelling between successive carrier periods.
It scales as L, exactly as (5). No signed material density, fictitious
extra spin, or cancellation of the surrounding-fluid exchange is used.

We impose the reference-row identity

    integral chi(r)b(r,s0)r³dr=Q_* !=0                        (7)

for every sheet on the tag support. The axial controls constructed
below preserve this pointwise identity. The local quadrupole chart
therefore never divides by a zero. The order-delta small Q_* required
by the physical spin/action matching has no additional1/L factor.

## 4. The finite-packet carrier factor must be derived

First take uncut Gaussian weights; their remote cutoffs are handled
with their differentiated exponentially small tails. The angle's
longitudinal integral involves the product

    exp[-L²(q-p)²/2] exp[-ell_z²(q-p_*)²/2].

Completing the square gives

    q_bar=a p+(1-a)p_*, a=1/(1+c²),
    var_q=1/(L²+ell_z²),
    prefactor=(1+c²)^(-1/2)exp[-D(p-p_*)²/2],
    D=L²ell_z²/(L²+ell_z²).                                  (8)

The radial angle polynomial and phase of0142 are therefore averaged
around q_bar, not p. Taylor's formula under this exact Gaussian
integral, with the finite number of weighted profile derivatives,
gives on the fixed optical window

    gamma_packet=-2Omega-(sqrt(2)/3)Omega sqrt(lambda/q_bar)
      +O_T(Omega[delta²polylog+delta/(pL)²]),

    p_*² partial_p² gamma_packet²|_(p_*)
      =a² sqrt(2)Omega²delta_*+O_T(Omega²[delta_*²polylog
                                                +delta_*/(p_*L)²]).
                                                                  (9)

The large real Gaussian amplitude factor in(8) cancels from the
phase derivative; it remains in the physical mass and parameter
connection. It is not dropped from the action. The finite positive
factor a² is part of the actual physical carrier response. As c is
fixed, this retains a strict fraction of0142's positive margin.
The same time-connection terms as0142(8)--(9) are computed from the
actual packet angle row, and are included in the physical equation.

## 5. Why using only the old radial controls fails, and the repair

There is a concrete new normalization failure if b is independent
of s0. Its spin carries the real prefactor exp[-D(p-p_*)²/2],
whereas Pi=-beta_packet/c carries its inverse. Their second carrier
derivatives differ by2D times the matched initial row. Asking only
radial derivatives to compensate that term makes the old moment
inverse grow as (pL)² and Q_* shrink accordingly. The resulting
reference dephasing is no longer below the curvature scale. Merely
taking L large does not repair this. This is a route-scoped failure,
not a failure of finite-action Euler packets.

Use two bounded even axial controls instead. Their zeroth and second
moments under chi_z(s)exp[-s²/(2L²)] are independent: before cutoff,
the basis {1,s²/L²} has positive determinant2v², where
v=c²/(1+c²)>0 is the normalized variance of s/L. A fixed remote
cutoff preserves this strict margin. Dualize these controls into
phi_0,phi_2 so their normalized moments M0,M2 equal (1,0),(0,1).
They are smooth and bounded on the compact marker support. Set

    b(r,s)=b_0(r)phi_0(s/L)+b_2(r)phi_2(s/L).                  (10)

Choose Q(b_0),Q(b_2) to be the two coefficients of the constant
function1 in this dual basis, times Q_*. This gives (7) EXACTLY.

The necessary radial functions are

    1, e^-x/2 P, e^-x/2 xP, D_r(e^-x/2 P), D_r(e^-x/2 xP),
    P=x²/2-5x+9, D_r=-1/2+(3/2)x partial_x.                    (11)

They are independent: the four polynomial rows after extracting
the exponential have rank4, while1 is independent of exponential
polynomials. Choose five radial points with a nonzero evaluation
minor and narrow smooth bumps. Their fixed exact moment matrix is
invertible. There is no solver-selected or floating-point rank.

At leading narrow-band order b_0 controls the zeroth spin and slow
time rows. b_2 controls their second carrier derivatives after
division by L². The first carrier derivative rows, scaled by p,
use D_r on b_0-b_2; this term follows by differentiating both the
carrier exponential and the first profile-envelope correction.
Thus the eight physical equations are: the two reference conditions
in(7), and the three carrier jets of spin/Pi and of their first slow
time moment. Add the harmless choices D_r P(b_2)=D_r(xP)(b_2)=0.
The resulting10 equations for the two sets of five radial bump
coefficients are block triangular with two copies of the invertible
matrix(11). Nonzero dimensional factors (including L² for the
second derivative and delta for the slow-time row) are divided out
before comparison.

The exact finite-band pressure/profile/physical-action equations,
including the Gaussian factor(8), perturb this dimensionless matrix
by O(delta polylog+(pL)^-1). Its inverse persists after p,L are
selected in that order. Solve the exact10 equations, using the
actual beta_packet and actual initial Euler derivatives. All targets
are derived from Pi, not assigned inertias. The two reference
conditions keep every sheet's quadrupole nonzero. Q_* of order
delta yields bounded coefficients uniformly in L; a common small
scaling enforces positivity in(6). This is the executed representation
repair for the large-D discrepancy, not an assumption that it vanishes.

As in0142, replacing Pi by a fixed eta Pi gives the same construction;
eta=1/2 is available when a coherent standing-pair preparation calls
for that normalization. The choice precedes the marker construction.

## 6. Finite time, error scale and exact action

The exact two Euler columns define a constant positive KKS beta_packet
and an actual collective angle row c(t,p),psi(t,p). Its exact scalar
pullback is

    L=M/2[(theta_dot-(c_t/c)theta)²-psi_t²theta²],
    M=-beta_packet/(psi_t c²)>0.                              (12)

The full packet phase has psi_t<0 on the fixed optical window. The
packet's large real carrier-envelope curvature remains in c and M.
The first slow-time matching and all three carrier jets are imposed
on the ACTUAL total physical tag spin. The same Taylor/Duhamel
argument as0142, now with bounded normalized axial moments and the
two-control inverse, gives a row error

    max_(j<=2) |p^j partial_p^j(S-eta Pi)|/Pi_scale
       <= C_T[delta² polylog+delta/(pL)²]+high-order Euler error.
                                                                  (13)

Time and carrier derivatives of Gaussian integrals are taken before
limits; the large factors L in envelope derivatives are retained in
the normalized matching matrix, not mistaken for profile derivatives.
Finite-order remainder bounds carry corresponding powers of L.
Since p,L,T and marker geometry are finite, first choose the order
of the material quasimode and its small delta sufficiently to make
their weighted Euler remainder smaller than the displayed physical
margin; alternatively use the explicit polynomial dependence in
pL and choose a fixed hierarchy pL=delta^-3 and a correspondingly
higher finite recursion order. The full generator's TG energy bound
is independent of p,L, so this ordering contains no exponential
large-box amplification.

The reference quadrupole dephasing is O(delta³)/Q_*=O(delta²),
with bounded axial-control factors, and not O(L delta²). Thus it
is below the signal in(9). Gaussian spatial cutoffs are taken after
the finite algebraic moments with their differentiated tail errors
below the same margin. The genuine nonlinear Euler disturbance
amplitude is chosen last for finite T as in0145.

This file supplies finite-action C² preparation families and actual
material histories. Toroidal closure, periodic approximation and
the exact large-radius pressure estimate are supplied separately;
it does not assert that a long straight packet already is a torus.
