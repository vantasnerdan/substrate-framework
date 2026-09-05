# Actual stationary cavity quadrupole and full carrier-two current

Use the genuine noncritical m=2 mode constructed in the companion proof,
at a selected sufficiently large FINITE positive axial carrier k0. Its
laboratory frequency omega>0, full KKS coefficient beta<0 and all
parameter derivatives are the actual radial Euler branch. Fix one
reference axial period containing its carrier; nearby k values have
the usual Bloch/fiber interpretation with fixed phase-averaged measure.
This is not an infinite monochromatic mode assigned a finite R3 norm.

## 1. A material observation where the actual base is quiet

The admitted profile vanishes identically for `r<r0-a`. Select a compact
radial tag support strictly in that cavity, and a finite axial interval
inside the reference period. This is a real stationary material region:
the background velocity is ZERO on a neighborhood of the whole support.
It is not the moving annular vortex material, and it is not a claimed
nonzero EPS core at the excluded central circle.

In the cavity the actual regular mode pressure is

    P(r,k)=A(k) I_m(k r), A(k)>0,
    xi=grad[P(r,k) exp(i m theta+i k z)]/omega(k)^2.     (1)

The positivity of A follows from the positive radial ground function
and the regular I_m solution. It cannot vanish: a zero cavity solution
would have zero pressure/displacement Cauchy data at the annulus and
make the entire eigenmode zero. All pressure tails remain in the full
mode and in beta. Formula (1) follows from BOTH Euler and Lin on the
quiet region; it does not prescribe an arbitrary potential flow there.

Choose even smooth `0<=chi_z<=1` of compact axial support so that

    Z(k)=integral chi_z exp(i k z)dz>0

on a fixed neighborhood of k0; support within `|z|<pi/(4 k_max)`
suffices. Put `Z0=integral chi_z`. With smooth `0<=chi_r<=1/2` and
bounded signed radial b_tag supported where chi_r is positive, take

    w0=chi_r(r)chi_z(z)[1+epsilon b_tag(r)cos(m theta)],
    B=integral chi_r b_tag r^(m+1)dr !=0.

Choose epsilon>0 after b_tag, so `|epsilon b_tag|<1`. Then this is a
nonnegative stationary material mass fraction. Its ordinary transverse
quadrupole is `Q0=rho pi epsilon Z0 B`. The measured physical angle is
`theta=delta(arg Q)/m`, with `Q=integral rho w(x+i y)^m` and m=2.
Its response to a rigid axial rotation of the physical tag is exactly
one. No helical scalar, rotating paint clock or Floquet relabeling is
used for this observation.

## 2. Full material angle, spin and displacement rows

Write the two real pressure phase coordinates as
`a_T=-omega b,b_T=omega a`. Direct angular integration of (1) gives

    Gobs(k)=integral chi_r r^m[P_r+mP/r]dr>0,
    Ttag(k)=integral chi_r b_tag r P(r,k)dr,
    delta Q=rho pi m Z(k)Gobs(k)(a+i b)/omega^2,
    theta=c(k)b,
    c(k)=Z(k)Gobs(k)/(omega^2 epsilon Z0 B).             (2)

The term with the painted modulation has no extra hidden contribution
to delta Q: angular orthogonality removes it; the unpainted radial
fraction supplies the displayed Gobs. The actual full phase action is

    M(k)=-beta(k)/(omega(k)c(k)^2)>0,
    L_theta=M(k)[theta_T^2-omega(k)^2 theta^2]/2.        (3)

This uses the mode's full action and the stationary measured clock.
Its mass derivatives are not put to zero.

The actual displacement moment and mechanical spin on the same tag are

    G_z=integral rho w0(r cross xi)_z
          =rho pi epsilon m Z(k)Ttag(k)b/omega^2,
    S_z=G_z,T.                                         (4)

The general material spin is
`integral rho w0[r cross xi_T+2xi cross u0]`; here the second term
vanishes because u0 really is zero on the tag, not because it was
neglected in a flowing region. The reference transverse centroid and
its first perturbation vanish for m=2 by angular orthogonality. Hence
(4) is the actual centered axial spin, not an origin-dependent impulse.

There is also actual axial tag motion. Formula (1) gives pointwise
`r xi_theta=(m/k)xi_z`, so the complete observed moments obey

    G_z=(m/k) delta Z_mass, S_z=(m/k) delta P_z.         (5)

Here delta Z_mass is the axial material first-position moment and
delta P_z its physical momentum. This linked translation is retained;
the mode is not advertised as a pure internal spin with zero momentum.
Its initial G is the actual particular Lin displacement (1). No
integration constant or additional rotor has been appended.

For any fixed tag with `B Ttag>0`, equations (2),(4) give

    G_z=J_tag(k)theta, S_z=J_tag(k)theta_T,
    J_tag(k)=rho pi epsilon^2 m Z0 B Ttag(k)/Gobs(k)>0. (6)

The axial overlap Z cancels in this literal current/angle ratio but
stays in the mass (3). The initially simplest positive b_tag therefore
already gives positive observed current on the entire mode history.
The next construction supplies one fixed tag with constant current
coefficient through the actual second carrier jet.

## 3. Four fixed radial controls, with no moving tag derivatives

Fix a positive reference B_* and a positive finite target ratio tau.
Solve these four linear equations for one b_tag:

    B=B_*,
    partial_k^j Ttag(k0)=tau partial_k^j Gobs(k0),
                                             j=0,1,2. (7)

These are physical tag moments, not frequency or action fits. The
actual pressure normalization A(k) and its derivatives occur in each
side of (7). Their common nonzero factor and triangular derivative
matrix can be canceled equivalently before choosing the radial controls;
none is silently held constant.

The four radial row functions are linearly independent on any open
cavity annulus. To prove this without a numerical small determinant,
remove that triangular A(k) factor. The rows are

    r^(m+1), r I_m(k0 r), r^2 I_m'(k0 r), r^3 I_m''(k0 r).

In the convergent Bessel series, the coefficient of
`r^(m+1+2j)` in the last three rows is a nonzero common coefficient
times a quadratic polynomial in `n=m+2j`, with derivative weights
`1,n/k0,n(n-1)/k0^2`. A linear relation on an open interval extends
analytically; for j>=1 it would make this quadratic vanish at every
such n. Already the three values j=1,2,3 have Vandermonde determinant
16 before the nonzero powers of k0 and series coefficients. Thus all
three coefficients vanish, then so does that of r^(m+1). This is an
exact independence proof, not a generic-rank assumption.

Select four interior evaluation radii with nonzero row determinant and
replace them by narrow smooth radial bumps. Continuity preserves the
determinant. The exact finite matrix inverse solves (7), with bounded
signed coefficients. A sufficiently small epsilon then ensures the
positive-density condition without changing these moment equations.
The reference Q0 is nonzero; the actual mode and its frequency have not
been altered to achieve the observation.

For this SAME stationary tag, with all derivatives taken at fixed
chi_r, chi_z and b_tag,

    J_tag(k)=Delta+O((k-k0)^3),
    Delta=rho pi epsilon^2 m Z0 B_* tau>0.             (8)

Consequently the entire fixed-time carrier-two histories obey

    G_z=Delta theta+O((k-k0)^3),
    S_z=Delta theta_T+O((k-k0)^3),
    G_z(0)+integral_0^T S_z=Delta theta(T)+O((k-k0)^3). (9)

Time derivatives of these remainders are controlled on every fixed
window by the analytic mode branch. At k0 the identities are exact
for all linear-mode times. The measured positive overlap is
`eta=Delta/M(k0)>0`; no unnecessary unit-overlap assumption is imposed.
Reducing a common physical tag fraction lowers Delta while its
normalized angle and full phase mass remain unchanged.

## 4. Conditioning, continuation and the exact licensed object

This construction orders its choices: first the admissible profile and
a sufficiently large but finite noncritical carrier; then the actual
radial mode, a fixed quiet tag region and axial window; then its finite
four-row matrix; finally epsilon for positive density. The pressure
amplitude in the cavity can be very small at high k. Its actual nonzero
amplitude, c(k0), the matrix inverse norm and the reference quadrupole
therefore remain observation-conditioning constants. They cannot be
omitted from a subsequent thin-ring or finite-time relative-error bound.
No uniform lower bound as k tends to infinity is claimed or needed for
this finite-parameter existence theorem.

The positive result is a smooth pressure-localizable straight annular
Euler field, an actual positive-action noncritical m=2 optical mode,
and a fixed positive stationary ordinary material quadrupole with full
displacement/spin/current identities through carrier two. Its actual
linked axial tag centroid/momentum is given by (5). The mode is finite
action per axial period/fiber; a finite-bandwidth R3 packet has separate
bandwidth errors. The tag lies in the stationary cavity, not in a claimed
nonzero EPS central core. Global closed-ring perturbation theory and
array/continuum interaction must retain all poloidal pressure channels.
The construction supplies those next routes with actual mode and
observation data; it does not claim to have already closed them.
