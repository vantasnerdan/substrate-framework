# Executed packing repair: positive density at fixed finite-time accuracy

The single-ring padding limit in the previous file dilutes the optical
population. This is the registered third candidate, motivated by that
specific failure. It preserves closed compact rings at every finite
radius, instead of multiplying a vanishing density by an unnamed factor.

## 1. Exact geometry and a count in physical volume

Fix the profile, finite carrier p0, packet length L, normal preparation
radius B and time window from the first file. Choose normal spacing
d>4B and arc spacing D sufficiently larger than L. For a large R,
take coaxial copies of the actual CLV ring with major radii
R_j=R+j d in [R,2R] and center heights z_l=l d in an interval of
length R. Each copy has its own epsilon=1/R_j and the SAME physical
core units and limiting profile. Their compact velocity/pressure
regions and their normal preparation tubes are disjoint. Their sum
is EXACT stationary Euler and is still globally axisymmetric.

On each ring place equally spaced rotated copies of the finite arc
preparation and its tag, with count M_j=floor(2pi R_j/D). The tags
are actually disjoint and nonnegative. The preparation is band-limited
around n/R_j=p0, not sharply truncated in arc length; its overlapping
Schwartz tails are retained. Its KKS overlap matrix is a small
perturbation of the diagonal one when D/L is sufficiently large.

Place the entire cluster in a periodic cube of side P=cR, where c is
a fixed sufficiently large geometric constant, for example c=16.
The compact clusters in different cubes are disjoint. The number of
marked arc packets satisfies

    N_R=3pi R^3/(d^2 D)+O_(d,D)(R^2),
    N_R/P^3 -> 3pi/(c^3 d^2 D)>0.                          (1)

This follows by summing 2pi R_j/D over the R/d radial positions and
R/d heights. Every packet has its actual finite tag mass and actual
full phase action; no amplitude proportional to the box volume enters.

## 2. Why normal pressure coupling is uniformly small, not absent

All harmonics used by these preparations satisfy |n|>=c0 p0 R.
In the annular region R/2<r<3R containing the cluster, the scalar
pressure/vector-potential operator has the positive angular term

    -Delta_n=-partial_r^2-r^-1 partial_r-partial_z^2+n^2/r^2,
    |n|/r >= c1 p0.                                       (2)

For a source supported in one fixed normal tube, test the full
equation with an exponential normal-distance weight. Choose its
logarithmic slope below c1 p0/4 inside that annular region and make
the weight constant beyond it. Cauchy--Schwarz and absorption in
the positive term (2) control both the potential and its gradient.
The escape across the annular boundary costs exp(-c2 p0 R).
Equivalently, the radial angular barrier has integral
|n| |log(r/r_source)|. It prevents a purported cheap pressure path
through the remote axis or infinity. Fixed derivative commutators
and the shifts n to n+/-1 for vector components preserve these bounds.

Consequently the FULL Biot--Savart/pressure map from one normal tube
to a disjoint tube at normal distance b has operator norm bounded by

    C_(s,p0,B) exp[-c3 p0(b-2B)]
                              +C_s R^a exp(-c2 p0 R),     (3)

in the corresponding finite Sobolev inventories. Polynomial factors
from finitely many derivatives can be absorbed by reducing c3.
The cylindrical volume measure is retained; r/R stays between fixed
positive bounds on every source and target tube. This is an estimate
for the whole-space inverse, not an artificial outer cylinder wall.

The linear vorticity equation is particularly useful here. Perturbed
vorticity stays in the disjoint invariant background support regions.
Its diagonal block is the actual isolated-ring Euler evolution; its
off-diagonal blocks are exactly the other tubes' Biot--Savart velocities
and gradients in `omega dot grad v-v dot grad omega`. They are bounded
operators in H^s by (3). Local transport/stretch has the uniform
coefficient energy bound, since disjointness keeps the background's
global C^s norm uniform in the number of rings. The row and column
sums of the first term of (3) over the normal two-dimensional lattice
are bounded by C exp[-c3 p0(d-2B)], independently of R. The total
escape contribution is still a polynomial times exp(-c2 p0 R).
The Schur test and Duhamel therefore bound actual cluster evolution
against the direct sum of actual isolated-ring histories uniformly
on the fixed time window. Nonlocal coupling is small but retained.

## 3. Arc localization and the weighted spatial-jet bound

Normal separation alone would not control common-lab derivatives:
two packets can be far apart along the same circle. The packet band
supplies the second estimate. In scaled angular frequency q=n/R,
the actual pressure-sector inverse on the fixed nonzero band has
uniform finite q derivatives by differentiating its resolvent and
using (2). The velocity/vorticity generator's q derivatives are
bounded on the required Sobolev inventory: angular advection is
multiplication by i q times a bounded coefficient, and differentiated
pressure inverses retain the same gap. Duhamel gives the corresponding
finite q derivatives of the actual propagator for finite time.
This continuous-q interpolation is an operator estimate; only its
integer n samples are used as physical ring fields.

Summation by parts in those smooth sampled spectral coefficients
therefore gives, for every selected finite A, a bound proportional to
(1+distance_arc/L)^(-A) for separated packet histories. The same
bound applies to the finite p and time derivatives. Poisson images
around the circle are included, with the shortest periodic arc distance.
No finite-speed propagation is asserted for Euler pressure.

Combining this with (3), the matrix of packet interaction/observation
errors E_ab has a uniform weighted Schur bound

    sup_a sum_b (1+|x_a-x_b|^3)|E_ab|
      <= C_T [exp(-c3 p0(d-2B))
                +C_(A,L) D^(-A+4)+R^a exp(-c2 p0 R)].     (4)

The same-ring overlap terms use the Schwartz arc bound. Normal
distances and their third powers are summable against (3). Taking A
large enough gives the displayed arc sum. The exact constants include
packet/tag conditioning and the finite p derivative inventory. The
third moment is stated so two spatial derivatives have a controlled
Taylor remainder, not just two formal coefficient sums.

Prepare a common laboratory input by multiplying each actual packet
at x_a by exp(i K dot x_a), using the local axis t_a in
p_a(K)=p0+t_a dot K. For an axial vector input its amplitude is
n_a dot Phi, with n_a the actual axial angle axis of that tag. These
are actual initial displacement columns; all discrete coefficients
remain in the same fixed high angular band. In cross forms the
phase derivative costs x_a-x_b, not the absolute cluster radius.
Thus (4) is exactly the estimate needed for the common-K derivative
action and observations. It is not an inference from integer n alone.

## 4. Periodic clusters, complete mean and a positive density bound

The cluster remains axisymmetric before its periodic repetition, so
the actual isolated-cluster histories retain their large angular band.
Their exterior spherical expansion starts at degree >=c0 p0 R.
Choose the fixed cube ratio c large enough that each other cluster
lies beyond a fixed larger concentric ball. The convergent harmonic
addition expansion gives an exterior bound with factor b^(-c0 p0 R),
b>1, times polynomial powers of R and the finite Sobolev norm. There
is no factorial constant substituted for this geometric-series bound.

The phased periodic-copy comparison of the preceding file applies.
Its projector/gauge derivative factors are only polynomial in P=cR,
so the remote-cluster error, through two ray derivatives, is still
bounded by R^a exp(-c4 p0 R). The Bloch zero harmonic is kept as
P_kappa exactly as there. Whole O(3)/parity averaging gives the
actual isotropic second jet, with the complete stress and shape rows.

Choose d and D once, after the prescribed finite-time physical error
margin and all packet conditioning constants, so (4)'s first two
terms are below that margin. Then send R large, without increasing
d,D,L or changing tag masses. Both the geometry error and the
remote-cluster error tend to zero. Initial cross forms, their actual
evolution and their derivative-level observations remain uniformly
close to the sum of the original positive packet contributions.

For example, if the relative error is below one quarter, the actual
whole-law phase and current densities have positive lower bounds
proportional to

    j >= const J_packet/(c^3 d^2 D)>0,
    Delta_density >= const Delta_tag/(c^3 d^2 D)>0.          (5)

The full Haar factor is 1/3 for both un-reconstructed action and
current; the physical angle estimator uses the inverse orientation
second moment, as in 0202. No energy of the linked polar companion
is subtracted. The entire ambient fluid still has density rho, not
the occupied tag fraction.

This is a density-preserving family at FIXED prescribed finite-time
accuracy. Tightening that accuracy can require larger d,D,L and can
change the bound (5); no uniform positive bound in a simultaneous
zero-error limit is asserted. At every selected finite R all cores
are actual compact closed rings. A local limit as R tends to infinity
can have straight, rather than closed, cores and is not named an EPS
limit. Likewise this theorem supplies positive finite-time action and
actual common-K jets, not an exact autonomous optical band or the
parent's full coupled continuum. It repairs the dilution mechanism
without hiding the remaining dynamical and topological distinctions.
