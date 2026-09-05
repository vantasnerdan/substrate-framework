# Actual finite-time annular packet on the same compact closed ring

This is a finite-time Euler/Lin construction, not a continuation of a
fixed-m pole through the other poloidal bands of a bent ring. Its inputs
are exactly 0198's smooth compact axisymmetric ring and 0200's actual
noncritical straight-column mode and stationary material tag.

## 1. Core scaling and a uniform global coefficient bound

Let epsilon be the thin-ring parameter in 0198. Translate a point on
the source's major circle to the origin, magnify all lengths by
1/epsilon and use the corresponding core time. After 0198's amplitude
normalization the major radius is R=1/epsilon and the core velocity
scale is fixed. Denote the actual velocity in these units by u_R.

The cutoff is a fixed smooth function of the scaled streamfunction,
supported in a fixed compact annulus of positive normal radius. On
that annulus all derivatives of the source's regular branch have their
ordinary Taylor bounds. The excluded central circle and both flat
cutoff edges have neighborhoods where the completed velocity is zero.
Axial symmetry makes the same bounds apply all around the circle.
The tubular metric has bounded inverse on the support for large R.
Consequently, for each fixed finite s,

    sup_(R>=R_s) ||u_R||_(C^s(R3)) <= C_s,                    (1)

and on every fixed Cartesian ball about the chosen circle point,

    ||u_R-u_infinity||_(C^s(B_D)) <= C_(D,s)/R.              (2)

Here u_infinity is the exact constrained annular column of 0200.
The support volume grows like R, but the coefficient norm in (1)
does not. This is the relevant distinction for finite-time Euler
energy estimates. No exponentially growing global norm is hidden in
the order of the geometric approximation.

## 2. A finite-action exact column preparation

Fix a sufficiently large but FINITE carrier p0, a small compact band
I around p0 lying inside 0200's noncritical branch, and a smooth band
cutoff equal to one on a smaller neighborhood of p0. All observation
conditioning constants are fixed at this stage. Put

    A_L(q-p)=L exp[-L^2(q-p)^2/2] chi_band(q),
    Xi_L(p,x,z)=1/sqrt(2pi) integral A_L(q-p)
                                  xi(q,x) exp(i q z)dq.     (3)

The two real phase columns are the real and imaginary parts. The
actual initial velocity is P_R3(Xi_L cross omega_infinity). Evolve
the exact Euler and Lin equations, equivalently superpose the actual
fiber histories with frequencies omega(q). The full-space pressure
commutes with the axial Fourier integral. Every fixed number of p and
time derivatives is finite in H^s and is Schwartz in z; the normal
tails are exponentially decreasing uniformly for q in I. Thus (3)
is finite action in R3, not a monochromatic fiber assigned a length.

Write beta_f(q)<0 for 0200's KKS per unit axial length, with its real
phase convention. Plancherel gives the full forms

    beta_L=integral |A_L(q-p)|^2 beta_f(q)dq <0,
    h_L=integral |A_L(q-p)|^2[-beta_f(q)omega(q)]dq >0.        (4)

Their leading normalization is sqrt(pi)L times the corresponding
fiber density, with all q derivatives retained. The actual phases
are not a single exact oscillator at finite L.

Keep the SAME fixed quiet-cavity tag from 0200, with its finite axial
window and four fixed radial moments. Its angle and full displacement
moment are linear observations of (3). Gaussian approximate-identity
estimates, differentiated under the actual integral, imply at p0

    max_(j<=2, l<=3) |partial_p^j partial_t^l
        [theta_L-theta_f]| <= C_T/L^2,
    max_(j<=2, l<=2) |partial_p^j partial_t^l
        [G_L-Delta theta_L]| <= C_T/L,                      (5)

on every prescribed finite time window. In the first line theta_f
means the actual p-dependent fiber observation with the same initial
real phases. The conservative second bound follows from
J_tag(q)-Delta=O((q-p0)^3), including the twice differentiated row;
no uniform matching outside those carrier jets is assumed. The linked
axial position/momentum similarly obey

    Z_L=(p/m)Delta theta_L+O_(C^2_p C^2_t)(1/L),
    P_L=partial_t Z_L, S_L=partial_t G_L.                    (6)

All constants include the real tag quadrupole and cavity pressure
amplitude. They can be large but are fixed before L is selected.

For the phase chart (theta0,r0), r0=theta_t(0), write
theta_L(0)=c_L b and theta_L,t(0)=d_L a. Both c_L,d_L are positive
and approach c(p),omega(p)c(p). The actual initial forms are

    J_L=-beta_L/(c_L d_L)>0,
    H_L=h_L diag(c_L^-2,d_L^-2).                            (7)

Thus J_L grows like L while the tag-current coefficient Delta remains
finite. The physical overlap is small but positive, of order 1/L;
neither the full action nor the tag mass is renormalized to hide it.

## 3. Exact high-angular-band Kelvin preparation on the ring

On the nonzero axial band, the full inverse curl supplies a smooth
vector potential for (3). Cut off this potential, not its curl, at a
large fixed normal radius B. Its discarded H^s tails, including the
p derivatives, can be made arbitrarily small. Keep the potential's
axial Fourier band. Periodize its axial coordinate around the circle
by sampling q_n=n/R, with the actual 1/R Riemann-sum normalization.
Equivalently use its Poisson sum with period 2pi R.

Push the potential as a one-form through the cylindrical tubular map.
Taking its physical curl is precisely the volume-preserving Piola
push of the displacement: D Phi Xi/det(D Phi). The Jacobian
1+x/R and all vector-frame derivatives are retained. The resulting
initial displacement Xi_R(p) is globally smooth, solenoidal and
supported in a fixed-width toroidal neighborhood. It has only angular
harmonics n with n/R in the fixed positive band I; the conjugate
negative band supplies its real columns. Choose R so every |n|>=N,
where N is a preselected finite integer, say N>=12.

This is not an integer-quantized carrier curve: p remains continuous
in the smooth coefficients A_L(n/R-p), on a fixed background and fixed
set of harmonics. Its derivatives are derivatives of actual initial
data. The actual Kelvin velocity is

    v_R(0)=P_R3(Xi_R cross curl u_R).                        (8)

After that preparation, evolve the full equations

    v_t=-P_R3[u_R dot grad v+v dot grad u_R],
    Xi_t+u_R dot grad Xi-Xi dot grad u_R=v.                  (9)

They preserve the actual Kelvin relation, not a finite-dimensional
mode constraint. Axial symmetry preserves the selected global
azimuthal band under this exact evolution.

## 4. Full pressure convergence, not a local-projector assertion

Both the column and the ring in (9) use the SAME full R3 projector.
For the velocity difference its forcing, apart from the initial-data
error, is exactly

    -P_R3[(u_R-u_infinity) dot grad v_L
                           +v_L dot grad(u_R-u_infinity)].  (10)

Fix s>5/2 and a finite derivative inventory, increasing s for the
time derivatives in (5). The standard differentiated energy identity
for linear Euler has coefficient depending only on the uniform finite
C^(s+2) bounds in (1), not on R or support volume. The norm of P on
every H^s is one. Splitting (10) into B_D and its complement gives

    ||forcing||_Hs <= C_(D,s)/R ||v_L||_(H^(s+1))
                  +C_s ||v_L||_(H^(s+1)(outside B_(D-1))).   (11)

Use a smooth spatial partition to interpret the final local Sobolev
norm. All derivatives of that fixed partition are included. The
second term is small uniformly in the prescribed time interval by
the actual fiber tail estimates; it is not an omitted pressure return.
The Lin difference has the same local/tail splitting and a direct
transport/stretch energy estimate. Parameter derivatives differentiate
only the initial data, so the same bounds apply to them without a new
operator assumption. Finite time derivatives follow from (9), with
the required additional spatial regularity.

Therefore first choose the finite L, then a tail tolerance and B,D,
then R large. The exact ring histories converge in the full mixed
C^2_p C^3_t H^s inventory to the finite column packet, with bound

    E_ring <= C_(T,s)[E_cut(B)+E_tail(D)+C_(D,s)/R
                                      +E_Poisson(R/L)].     (12)

Poisson errors decrease faster than each inverse power of R/L. The
constants may depend on the already fixed L,B,D,p0 and tag, but the
energy amplification is uniform in R. Equation (10) is the nonlocal
pressure license absent from an argument based only on (2).

## 5. Actual stationary material observation and both actions

The fixed tag can be chosen strictly inside the quiet cavity with a
positive normal margin. For all sufficiently large R its whole support
still has u_R=0. It is an actual stationary material mass fraction in
the reference ring, measured in the fixed laboratory transverse plane
at the selected circle point. It is not identified with flowing
annular vorticity, an EPS central streamline or a rotated paint clock.

Ordinary angle, centered axial G and S, the axial centroid/polar
moments and the entire symmetric displacement/shape tensors are taken
from the actual Xi_R,v_R and this tag. The exact moving-tag formula
gives S=G_t here because the actual base velocity on the tag is zero.
Changes of transverse centroid in the perturbed ring are retained;
they alter a centered linear moment by the explicit reference-centroid
terms, and their difference from the straight m=2 zero row is bounded
by (12). The nonzero reference quadrupole makes the angle map smooth.

The full KKS and stationary Euler Hessian also converge under (12).
For example, with A_R Xi=P(Xi cross omega_R)-[u_R,Xi],

    Omega_R(Xi,Eta)=rho integral omega_R dot(Xi cross Eta),
    H_R(Xi,Eta)=-Omega_R(Xi,A_R Eta).                        (13)

These forms include the induced exterior velocity through A_R. They
are the actual conserved forms on the two initial phase columns.
Their reference positive restriction persists for sufficiently small
conditioned error; no positive coadjoint norm from a different field
is substituted.

Separately, let the actual observed history in initial angle/rate
coordinates be theta=F(t)theta0+Gfun(t)r0, with its actual constant
phase coefficient J>0. Put

    W=F Gfun_t-F_t Gfun,
    N=F_t Gfun_tt-F_tt Gfun_t,
    M_obs=J/W, K_obs=J N/W^2.                               (14)

The exact two-column moving-observation action is
`(M_obs theta_t^2-K_obs theta^2)/2`, up to its computed boundary
term. W,N stay positive by (5),(12), so both coefficients are
positive on the prescribed finite interval. They need not be constant.
This moving-embedding action is distinct from the conserved initial
Euler Hessian (13); closeness of both is proved, equality to an exact
autonomous oscillator is not asserted. Their actual carrier derivatives
and the linked polar energy remain in all later averages.

The next file supplies actual common-laboratory Bloch preparations on
a periodic ring array and tracks the density and complete ambient row.
