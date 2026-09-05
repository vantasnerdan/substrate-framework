# Conditional compact assembly, actual chart, and fixed positive density

The input is exactly the finite-R compact stationary solution frozen in
README. It is still being constructed in 0263. This file proves the
downstream geometric and chart implications without assuming the unresolved
joint acoustic Schur coefficient from 0250.

## Exact disjoint assembly

Let the compact velocity and pressure, with exterior pressure set to zero,
be u_R,p_R. Choose a cubic cell of side L>2(R+a+delta), where a+delta
bounds the meridional support radius. Periodize translated copies:

    u_per(x)=sum_{m in Z^3}u_R(x-Lm),
    p_per(x)=sum_{m in Z^3}p_R(x-Lm).                    (1)

The supports are disjoint, all physical derivatives vanish on their
boundaries, and the sum is locally finite. Thus div u_per=0 and
(u_per dot grad)u_per=-grad p_per exactly: every nonlinear cross term
vanishes pointwise. This is compact-velocity gluing, which would have
failed for the noncompact velocity tails of the isolated 0211 ring.
Constant positive mass density fills the entire cell, including its
quiescent part. Equation (1) is one smooth stationary Euler field.

Choose the positive core tube inside one regular level curve of phi_R.
Its meridional center revolves into a closed, nonzero, Euclidean circular
vortex line. The tube is bounded in one Euclidean lift of the cell and
contractible in the periodic ambient space. It is therefore a different
topological supplier from C-CST-017's noncontractible axial tube.

## Actual frequencies and a volume action-angle chart

On a fixed regular annulus write psi_R=R phi_R and F=lambda psi_R,
with g=0. The exact velocities are

    |u_pol|=R|grad phi_R|/r,
    dot varphi=lambda R phi_R/r^2.                     (2)

For C_R(h)={phi_R=h}, define

    T_R(h)=integral_C r dl/(R|grad phi_R|),
    Omega_1=2pi/T_R,
    Delta varphi=lambda h integral_C dl/(r|grad phi_R|),
    Omega_2=Delta varphi/T_R.                           (3)

Choose theta1=Omega_1 times elapsed poloidal travel time, and subtract
from varphi the periodic primitive of
(dot varphi-Omega_2)/Omega_1. The resulting theta2 is periodic and gives
exactly u=Omega_1(h)partial_theta1+Omega_2(h)partial_theta2. The primitive
exists because its angular mean is zero by (3). Axisymmetry removes
theta2 dependence from the coordinate coefficients.

The Euclidean volume form is r dr dz dvarphi. Coarea and (2) give

    dVol=(R/Omega_1(h)) |dh| dtheta1 dtheta2.            (4)

The shear of theta2 has determinant one and does not change this form.
Any monotone action I(h), including the physical flux action, therefore
gives J(I)dI dtheta1 dtheta2 with J positive and independent of angles.
This is the actual volume chart required by the vector Kelvin lift.

Frequency variation is checked separately from return-map twist. In the
straight core limit Phi(s)=c J0(lambda s),

    Omega_1,0(s)=c lambda J1(lambda s)/s
       =c lambda^2/2-c lambda^4 s^2/16+O(s^4),
    d_s Omega_1,0=-c lambda^4 s/8+O(s^3).               (5)

Choose a small fixed annulus 0<s_-<s_+ before R. Its frequency derivative
has a strict nonzero margin. Differentiation of the regular contour
quadratures (3) preserves this margin for sufficiently large finite R.
Also Omega_2 is positive and bounded above and below by positive
constants times 1/R on this annulus. Thus the acoustic band Omega_1(I)
and the optical band Omega_1(I)-Omega_2(I) overlap on a nonempty compact
subinterval for one sufficiently large R. This proves actual frequency
overlap, independently of a Floquet-stability assertion.

For fixed toroidal index n=-1, every poloidal Fourier denominator is
l Omega_1-Omega_2. With R large enough, l!=0 gives a lower bound
|l|inf Omega_1/2, and l=0 gives inf Omega_2>0. For n=0 the actual
zero-mean sources in 0260 have no l=0 term and the same l!=0 bound.
Consequently the exact cohomological inverse, including its angular
shear term, exists on this same fixed annulus with finite smooth costs.
These constants can grow with the chosen finite R. No radius grows
with preparation accuracy.

## Exact core identity and the pressure-transfer boundary

Because g=0 and F=lambda psi_R, the vorticity is exactly lambda u on
the selected core, even though phi_R is not the literal circular profile.
The actual I in (4) is invariant under that field. Hence the two neutral
tangent polarizations and their isovortical lift use an exact eikonal and
an exact local constant-curl identity, not a nearby periodic Beltrami
approximation or an assumed persistence of an invariant action.

The periodic Euler projector differs from the whole-space projector by a
smooth image kernel on two source neighborhoods contained in this cell.
Indeed its local diagonal singularity is the same fundamental-solution
singularity; after that term is subtracted, the remaining kernel solves
the homogeneous elliptic equation locally. Equivalently, periodic
Fourier multipliers give the full parametrix directly. Every image-kernel
pairing with a compact high-action carrier can be integrated by parts in
I to arbitrary fixed order, with the actual smooth kernel derivatives
in the constants. The periodic zero mode is retained as its explicit
ray-wise projector; its smooth source mean has the same oscillatory
estimate. The action tag's resonant coefficient is kept separately.

This supplies the correct full-pressure setting for the finite-window
recursion. It does not turn a local leading gain into an invertible
joint physical gain: that determinant is still evaluated with the same
Bloch phase, material moments and all odd couplings. The outstanding
0250 Schur computation is a real supplier obligation here.

## Fixed positive law and twist

Take uniform translation phase in the cell and Haar rotations/reflections
of the whole field, cell, tags and source law. Add the stipulated
time-reversal pairing. Translation averaging and steady advection of
the action tags make the law stationary; the whole rotation average
is isotropic. Positive fractions in the tag law stay positive. A fixed
nontrivial smooth tag chi supported in the invariant tube has

    tagged volume fraction = integral_cell chi dVol/L^3 >0.     (6)

The geometric tube volume and its core length are finite and positive,
with length density 2pi r_core/L^3. Both (6) and this length density
are constants after R,L and the tag are fixed; neither tends to zero
along a later preparation-accuracy sequence. A literal positive
per-tube measured coefficient j_tube, once derived by the physical
supplier, gives isotropic density E[f] j_tube/(3 L^3)>0 in the vector
Haar convention. This last formula preserves the measured-coefficient
hypothesis rather than identifying volume with inertia.

The exact flux-action twist is the 0258 contour formula on this same
field. Its straight-core limit dq/dJ_u=lambda^2/(8c) is nonzero; the
fixed-annulus derivative persists by the assumed core convergence.
One may select a Diophantine boundary rotation in its interval image.
The inner g=0 equation is analytic and uniformly elliptic away from
the axis, so the local analytic persistence license uses this actual
core and its positive flux form. No analytic extension across the
flat compact-support boundary is asserted or needed.

`route_verdict: established as stated for the conditional exact compact
assembly, actual chart/frequency overlap, cohomological source domain,
and fixed positive geometric density`

`evidence_scope: conditional on the finite-R stationary supplier; full
joint material gain, measured-coefficient normalization and remaining
source/action/current joins retain their separate constructions`
