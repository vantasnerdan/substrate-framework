# Constructive actual passive-current controls on the same Euler field

This is a fixed-window physical-output construction. Its actual energy
and remaining action/optical obligations are stated at the end. The
field remains0188's A=1/100, B=1 finite periodic stationary Euler field.

## 1. An exact Euler invariant sector with nonzero measured current

For any smooth periodic g(Y,Z), the perturbation (g,0,0) is solenoidal.
The full Euler pressure for this perturbation vanishes: X independence
gives div[Tg,0,0]=0 and (g,0,0).grad u=0. Thus

    exp(tL)(g,0,0)=(exp(-tT)g,0,0)                  (1)

EXACTLY, not through a finite-mode approximation. Let
f=kappaY VX+kappaX VY and prepare the correlated return

    r_V=f(g,0,0).

Whole-law contraction gives E[f²]=1/5 and cancels the other components
of c_X against f. Its actual observed current is therefore

    P_g(t)=(A/5)<sinZ exp(-tT)g>.                    (2)

The microscopic field producing(2) is explicitly known and carries its
own energy. No harmonic is declared to be a physical angle or rotor.

## 2. Actual low frequencies from the existing wrapped streamlines

On the component Y in(0,pi), fix c=psi with |c|<B-A. Then

    cosY=(c-AcosZ)/B,
    Z_t=-sqrt[B²-(c-AcosZ)²],
    period(c)=integral_0^(2pi) dZ/sqrt[B²-(c-AcosZ)²]. (3)

These are actual streamline equations of T. For c approaching
c_s=B-A from below, put delta=c_s-c and z=Z-pi. The speed squared is

    2B delta+AB z²+O(delta²+delta z²+z^4).

The local period integral diverges logarithmically. Hence
omega(c)=2pi/period(c)>0 tends continuously to0. In the neighborhood
c>A, period'(c)>0 directly from differentiating(3); here B>2A, so that
neighborhood exists. Thus every sufficiently small positive frequency
is attained on a regular wrapped streamline, away from the separatrix.

Choose the angle theta with theta_t=omega(c), theta=0 at Z=0. The
symmetry of the speed in(3) gives Z=-pi at theta=pi. Thus sinZ<0 for
0<theta<pi, and its first sine coefficient

    s1(c)=2<sinZ(c,theta) sin theta>_theta

is strictly negative. The second half has the same sign in the product.
This proves a genuine physical observation overlap. It is not an assumed
nonzero matrix row. Both omega and s1 are smooth on any closed regular
streamline band chosen before the long-wave limit.

## 3. Finite-width smooth packets and exact energy

The normalized coarea measure on this component is

    dmu=mu(c) dc dtheta/(2pi), mu(c)=period(c)/(4pi²).

Take a smooth compactly supported nonnegative band density eta(c) with
integral eta dc=1, and choose

    g(c,theta)=G(c)cos theta,
    G(c)=10 w eta(c)/[A mu(c)s1(c)].                  (4)

It extends by zero to a smooth periodic finite-energy field since the
support avoids the band boundary and the separatrix. Exact angular
orthogonality in(2), not a mode-frequency fit, gives

    P_g(t)=w integral eta(c) sin[omega(c)t] dc.        (5)

Its physical mean is zero; its initial current is zero. The actual
whole-law squared velocity norm is

    N_g=E||r_V||²
       =(10 w²/A²) integral eta(c)²/[mu(c)s1(c)²] dc. (6)

This is finite for every selected finite band, but need not be small.
Disjoint bands have additive energies. The initial V gradient kinetic
energy is rho k² sum N_g/2 and is retained. Whole-field TR makes the
return even while the D material rate is odd, so their averaged initial
cross energy vanishes for this coherent selection. The even gradient
mass in(6) does not vanish under pairing.

The exact initial material V configuration remains zero, and each packet
has zero mean. Therefore the full initial common phase is still rho J;
this does not identify the physical measured velocity momentum with the
complete canonical momentum after evolution.

## 4. Constructive fixed-window output approximation

Fix T, a finite derivative order r, and a desired odd smooth output H(t)
on [-T,T]. The intended H is the DIFFERENCE between the target common-V
current and the actual uncorrected Euler output(8) of
`actual-velocity-output.md`, not an empirical comparator. Whole-field TR
makes this difference odd. Approximate it in C^r by an odd polynomial
p(t)=sum_(j=0)^m p_j t^(2j+1). This is constructive, for example by
polynomial approximation to its highest derivative followed by integration
with the actual initial derivatives and odd symmetrization.

Choose omega_l=l epsilon, l=1,...,m+1. The exact Vandermonde equations

    sum_l w_l (-1)^j omega_l^(2j+1)/(2j+1)!=p_j,
    j=0,...,m,

have a unique solution. For fixed p and m their weights are
O(epsilon^(-(2m+1))). Taylor's remainder gives

    ||sum_l w_l sin(omega_l t)-p(t)||_(C^r[-T,T])
       <=C(p,m,r,T) epsilon²,                       (7)

after increasing m to cover r if needed. Every selected frequency lies
in the actual low-frequency streamline interval for sufficiently small
positive epsilon. This does not create new oscillator input dynamics.

Choose disjoint regular bands around their streamline values. For
|omega(c)-omega_l|<=delta_l, differentiation of sin(omega t) gives
the explicit band error at derivative order j

    |w_l| delta_l [j Omega_max^(j-1)+T Omega_max^j],  (8)

with the j=0 bound |w_l|delta_l T. The band widths can be selected after
the finite weights so that the sum of(8) is as small as prescribed.
Each resulting eta,G is smooth and the energy(6) is finite. Equations
(1)–(8) construct actual Euler initial data giving any prescribed
fixed-window C^r output accuracy. They neither assert an exact linear
output for all time nor conceal the norm cost of this approximation.

Applied to H=-a t-R_V,bare(t), this supplies a genuine common-V physical
mean output accompanying0188's positive D coefficient, to the stated
controlled finite-window accuracy. The inherited phase and energy are
those of the actual packet preparation, not those of a guessed wave.

## 5. Actual tag separation and the remaining action join

Choose all packet bands below c_s=B-A. A material tag wholly in an
elliptic core psi>c_s+delta0 stays there because psi is conserved along
the stationary flow. The added velocity(1) vanishes identically on its
trajectories. The zero-initial-displacement Lin response to this added
velocity therefore vanishes there as well. At first spatial order, its
literal material angle, spin, displacement dipole G, and local shape
rows on that tag are all ZERO. This is exact support/transport separation,
not an identification of a canonical scalar with a tag angle.

The full finite-K pressure can affect higher spatial rows; no claim that
the complete finite-K perturbation vanishes on the tag is made. Those
rows are controlled only after the finite packet field and cell are fixed.
No same-EPS spectral or homogenization transfer follows from(1) alone.

Most importantly, the output approximation does not make the added
gradient mass(6) or the integrated-current connection disappear. A full
joint action must either retain that actual state or show that its energy
and phase corrections are controlled at the chosen physical branch and
packet-own scales. A small macro k can make its absolute gradient energy
small, but is not by itself an equality of second-jet action coefficients.
That is the exact remaining join, rather than an isolated-mean autonomy
condition added to the parent objective.
