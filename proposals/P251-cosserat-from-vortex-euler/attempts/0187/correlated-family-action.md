# Coherent cancellation with positive microscopic action and actual energy

This is the algebraic interface for actual families of stationary Euler
optical modes. It does not supply their same-field existence. It changes
the microscopic preparation, not the user's macro objective: each family
receives the SAME macro angle/rate inputs, but its initial measured angle
and rate are multiplied by a specified signed response amplitude.

## Complete common-input preparation

Let two families r=1,2 have positive probability w_r, sum w_r=1, positive
physical phase masses j_r(p), and a common nonzero physical frequency
nu>0 at their reference carriers p_r. Write v_r=nu_r'(p_r),
b_r=nu_r''(p_r). Use a smooth real amplitude A_r(p), whose first three
jets are a_r,d_r,e_r. These are prepared state correlations, not signed
probabilities. Actual real signed carrier bands and whole-field O(3)
averaging are those of0174. In a transverse/longitudinal scalar channel
their common-K factor is c_s=1/5 or3/5.

For the SAME initial macro pair(q0,r0), a local mode is prepared as

    theta_r(t,p)=A_r(p)[cos(nu_r(p)t)q0
                         +sin(nu_r(p)t)r0/nu_r(p)].       (1)

The macro measured angle is the probability average with the geometric
axis reconstruction, not the action-weighted average. Choose

    sum w_r a_r=1, sum w_r e_r=0.                         (2)

The latter avoids an unneeded second-gradient initial observation filter;
it is realized by locally linear A_r. The physical angle/rate is then
initially exactly(q0,r0) through K². Unit rigid-rotation response of the
tags is unchanged; opposite signed OPTICAL perturbations are not a rigid
rotation of all material. Thus this is not the older special preparation
requiring every local optical angle to equal the macro angle.

## Observed rows and the complete initial phase

For h(nu)=cos(nu t) or sin(nu t)/nu, actual carrier differentiation gives

    sum w_r (A_r h(nu_r))'' = B h_nu+V h_nunu,
    V=sum w_r a_r v_r²,
    B=sum w_r(a_r b_r+2d_r v_r).                         (3)

If V=0, both observed rows are precisely the Taylor rows of a stationary
oscillator with nu_eff=nu+c_s K² B/2 through K². Consequently W=1 and

    nu_eff²=nu²+c_s K² nu B+O(K³).                      (4)

This removes the second-jet temporal variance, not just its initial
Taylor coefficients in time. The expansion is uniform on fixed time
windows with its needed derivatives. It is not an all-K closure.

The inherited phase is averaged BEFORE eliminating these rows. Put
h_r=j_r A_r². The exact restricted initial form has mass

    J(K)=J0+c_s K² J2/2,
    J0=sum w_r j_r a_r²>0, J2=sum w_r h_r''.             (5)

Signed amplitudes are squared here. Equations(4),(5) give the stationary
physical action M=J(K), stiffness J(K)nu_eff² through K². All gradient
inertia stays present and is transformed together with the potential.

## The full physical Euler energy is a separate exposing check

A moving prepared-history action alone need not identify its Hamiltonian
with Euler energy. At time zero the ACTUAL sum of mode energies is

    E=1/2 sum w_r h_r(p)[r0²+nu_r(p)² q0²].             (6)

After the common-K expansion, it equals the stationary observed action's
energy if and only if the following additional scalar vanishes:

    D=sum w_r {h_r(v_r²+nu b_r)+2nu h_r' v_r}
                                            -J0 nu B=0. (7)

This retains derivatives of both mode mass and prepared amplitude.
When(7) holds, equality continues on the actual histories through K²:
Euler energy is conserved and the measured rows obey(4). No energy is
discarded by identifying an arbitrary observed frequency with a modulus.
This licenses an inherited-energy second-jet response on this prepared
class, not a globally invariant microscopic two-plane at finite K.

## Constructive signed correlation and positive curvature

For v1²!=v2², (2) and V=0 uniquely give

    w1 a1=-v2²/(v1²-v2²),
    w2 a2= v1²/(v1²-v2²).                              (8)

One amplitude is negative, while J0 remains strictly positive. Choose
both v_r nonzero. Conditions B=B_* and D=0 are two LINEAR equations for
d1,d2 after j_r,j_r',a_r,v_r,b_r are fixed. Their coefficient determinant
is a nonzero multiple of

    w1 w2 nu v1 v2 (j2 a2-j1 a1).                      (9)

The opposite signs in(8), j_r>0, make this nonzero. Hence any specified
finite B_*>0 is realized by finite initial preparation slopes, without
changing probability, adding a force, or fitting an empirical comparator.
The curvature c_s nu B_* is a PREPARATION-DEPENDENT derived coefficient;
this argument does not predict a unique intrinsic material modulus.
Carrier dependence of A encodes a real spatial-gradient preparation and
is part of the stated initial data, not an invisible fitting constant.

An exact algebraic example uses w1=w2=1/2, j1=j2=1, j1'=j2'=0,
nu=1, (v1,v2)=(2,1), b1=b2=0. Then

    (a1,a2)=(-2/3,8/3), (d1,d2)=(9/20,1/10),
    B_*=1, J0=34/9, J2=17/80.

These numbers exercise the interface; they are NOT asserted to be the
jets of an actual Euler branch. Actual families must supply those jets.

## Actual spin/current and remaining geometry

Suppose an actual stationary tag in each family has the0181 identities
G_r=Delta_r theta_r, S_r=Delta_r theta_r,t through carrier two, with
constant Delta_r>0. Equalizing Delta_r=Delta>0 by reducing positive tag
fractions preserves these identities and their phase masses are recomputed
in the correspondingly normalized angle coordinates. Then the whole-law
current is exactly G=Delta q, S=Delta q_t at this order, despite signed
local optical amplitudes. Its physical overlap is eta=Delta/J0>0;0182
keeps that overlap in the physical translation reconstruction. Equal tag
moments do not assert equal phase masses or erase the axial shape terms.

This suggests actual0181 branches with equal gaps after Euler amplitude
scaling and unequal slopes after geometric scaling. Two geometrically
scaled copies suffice here because the amplitude derivatives in(7),
unlike a frozen-amplitude variance subtraction, remain independent
preparation variables. The initial registered comparison of different
dimensionless branches remains an alternative. Actual mode normalization,
current/return maps and admissible common-K potentials must be supplied
before applying this interface to a same-field continuum.

The positive response/energy algebra is established as stated. Exact
Euler family realization, finite-action packet transfer, coupled acoustic
memory and stationary EPS same-field geometry remain active obligations;
this formula is not a replacement for any of them.
