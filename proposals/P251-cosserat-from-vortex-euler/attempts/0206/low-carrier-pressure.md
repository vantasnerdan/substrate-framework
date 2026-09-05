# Low-carrier pressure representation and the exact Euclidean reduction

## 1. The bounded axisymmetric slow operator

At k=n/R tending to zero, use the actual velocity components (a,b)
=(v_theta,v_z), not a material displacement whose axial component is
of order 1/k. On the active disk the limiting radial velocity is
v_r=k h, where

    h(r)=-i/r integral_0^r s b(s) ds,
    P(r)=-integral_r^a_c 2Omega(s)a(s) ds.          (1)

The pressure condition P(a_c)=0 is the LIMIT of the actual exterior
problem: the nonzero exterior return is of order k² log k in this
pressure row. It is not an impermeable radial wall. With tau=k t,
the complete leading equations are

    a_tau=-iW a-Z h,
    b_tau=-iW b-W' h-iP.                           (2)

Every W' and pressure term is present. The maps in (1) are compact on
the physical weighted L²(r dr) space: their integral kernels are
Hilbert--Schmidt on the bounded disk, or obtained as limits of those
kernels with the regular axis condition. Thus (2) is multiplication
by -iW plus a compact matrix integral operator. Its essential transport
range is -i[0,W0].

Substitution a=-Z f, b=P/(c-W)-W'f and h=-i(c-W)f
reduces its eigenvalue -ic to precisely the k=0 Sturm problem in 0201.
In particular its principal c0>W0 is simple, separated from the ENTIRE
slow transport interval, and has positive inherited scaled phase energy.
This is the actual unbounded-column long-wave operator, not an assigned
local sound speed or a continuation at fixed k.

For a physical solenoidal m=0 vorticity, the cross-sectional integral
of axial vorticity is zero at every k!=0. At k=0 retain that limiting
constraint. The transverse vector integral is zero by angular symmetry.
These two moment cancellations are what permit (1); allowing a new
circulation-changing exterior harmonic would define a different limit.

## 2. What exact symmetry removes on the n=1 ring

On the actual global ring, let E_R be the real four-dimensional span
of transverse translations and transverse rotations. The differentiated
solution family in `low-harmonic-action.md` proves that E_R is invariant
under the full Euler generator, with its exact translation/tilt Jordan
block. Its KKS matrix is nondegenerate: the translation/rotation block
is I0 times a planar rotation, and the rotation/rotation block is the
actual J0. No assertion that the latter vanishes is made.

Therefore the KKS-orthogonal complement E_R^perp is also invariant.
Indeed for x in that complement and e in E_R, Hamiltonian skewness
gives Omega(Lx,e)=-Omega(x,Le)=0. The exact projection onto E_R is
obtained by inverting its finite KKS matrix against the four rows
Omega(e_j,x). This preserves the entire exterior response: its entries
are the actual Euclidean impulse/angular-impulse moments. It does not
mean-centre a truncated velocity or freeze an external torque.

This quotient removes the m=+/-1 translation zero modes in the straight
limit BEFORE attempting a fast inverse. All other nonaxisymmetric
straight channels are separated from zero. The m=1 radial zero is the
unique translation proved in 0201; m>=2 and m<=-2 are excluded at zero
by the full radial boundary identity. The compact transport-resolvent
argument then supplies a bounded fast inverse on that quotient.

For n=2 there is no corresponding exact invariant Euclidean subspace.
The transverse bending coordinates remain in the slow block. Reusing
the n=1 quotient for n=2 would remove physical degrees of freedom.

## 3. The actual Schur target, and the terms which need bounding

Let S be the axisymmetric profile space (1)-(2), and F the remaining
space after the exact n=1 Euclidean reduction. The full physical
ring operator in these coordinates has blocks L_SS,L_SF,L_FS,L_FF.
For Re z initially in the genuine Euler resolvent half-plane, its
exact effective equation is

    [z-L_SS-L_SF(z-L_FF)^(-1)L_FS]s=0.             (3)

This is a FULL pressure expression. An independent assertion that the
bare L_SS is the entire operator would discard its last term. In the
time domain the same object is the Volterra kernel

    K_R(t)=L_SF exp(t L_FF)L_FS,                    (4)

plus the actual prepared complement term L_SF exp(t L_FF)f0. A finite
window proof can estimate these response rows directly; it need not
prove whole-spectrum stability or an all-time isolated pole.

The concrete sufficient scaled limit for the n=1 pole is

    R L_SS -> A0,
    R L_SF(z/R-L_FF)^(-1)L_FS -> 0                 (5)

uniformly on the FIXED contour around -ic0 of (2), together with the
actual phase/current pullbacks. Here A0 is the right side of (2) with
n=1. These are not consequences of the old fixed-k estimate: the
old pressure constant depends on a positive carrier lower bound.

There is a constructive way to expose the remaining estimates. Expand
the EXACT 0186 kernel at fixed integer n, not at a fractional angular
harmonic. In the intermediate range, the above zero vorticity moments
remove the monopole; the second-order axial/curvature remainder then
has the logarithmic integral integral_a^R ds/s. Terms linear in geometric
curvature are odd under simultaneous inversion of both meridional
coordinates. The leading background/border correction of 0195 has that
same odd parity. Its direct S-to-S average therefore vanishes, while
its S-to-F and F-to-S blocks remain and belong in (3).

This identifies the hoped-for bounds

    L_SS=R^-1 A0+O(R^-2 log²R),
    ||L_SF||+||L_FS||=O(R^-1 log R)                (6)

after the correct bounded identification of the Euclidean quotient.
Equation (6) is an ESTIMATE TO COMPLETE, not a measured fact: the
singular divergence reconstruction, moving source Jacobian and the
finite KKS projection must be controlled in the same operator norm.
The zero-moment/parity argument by itself does not prove all those
norms. This is the route's exact current pressure remainder, rather
than a general statement that the low-carrier mode cannot exist.

## 4. Why the observation needs its own same-order correction

The 0201 displacement of the limiting principal mode has xi_z=O(R).
Consequently multiplying an O(log R/R) background correction by the
displacement does not automatically give a vanishing global spin error.
The global position carries another R lever arm. The actual material
centroid, covariance and spin rows must therefore use the corrected
finite-ring mode and the exact Euclidean projection together.

For the full invariant tube fraction, the exact current remains

    delta S_D=partial_t G_D+2rho integral_D xi cross u,
    G_D=rho integral_D (x-X_D) cross(xi-delta X_D).                 (7)

The mean u over this actual stationary material domain is zero in its
physical translating frame; the actual global impulse is not its mass
centroid momentum. Using only the limiting local m=0 torsion in (7)
would miss the curvature and orbital returns responsible for any global
transverse response. Conversely the nonzero section flux in 0201 is
not already that material current. Equations (3)-(7) name the full
operator/observation calculation activated by the positive slow pole.

The main positive results at this boundary are the exact same-Euler
Euclidean phase/current reduction, the explicit compact slow radial
operator with a separated positive pole, and the actual positive
nonrigid finite-core bending FORM. The full low-carrier ring pole and
its positive measured material spin remain the specific construction
in (5)-(7); no all-time or finite-window solution is inferred from
the restricted bending frequency alone.
