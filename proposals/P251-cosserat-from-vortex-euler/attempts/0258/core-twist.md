# Exact retained-core twist, with physical flux normalization

Choose c,lambda>0 and restrict to a small regular annulus in the exact0253
core Phi(rho)=c J0(lambda rho). Orient the poloidal angle along its flow so
its frequency is positive. The straight velocity has axial speed lambda Phi
and poloidal speed -Phi'=c lambda J1(lambda rho). Its auxiliary axial circle
has length2pi R. Its rotation number, poloidal turns per axial turn, is

    rot_R(rho)=R q(rho),
    q(rho)=J1(lambda rho)/(rho J0(lambda rho)).             (1)

This circle is a straight periodic surrogate, not a compact Euclidean ring.
The finite-ring transfer below uses the actual cylindrical field equations.

Define the physical velocity flux through the enclosed cross-section and
its normalized action by

    F_u(rho)=integral_disk u_axial dx dz
             =2pi c rho J1(lambda rho),
    J_u=F_u/(2pi)=c rho J1(lambda rho).                     (2)

Because omega=lambda u on the core, F_omega=lambda F_u and
J_omega=lambda J_u. Neither flux is unweighted geometric area.

The Bessel series give

    q(rho)=lambda/2+lambda^3 rho^2/16+O(rho^4),
    J_u(rho)=c lambda rho^2/2-c lambda^3 rho^4/16+O(rho^6).
                                                                  (3)

Consequently the core limits of the normalized twist are

    dq/dJ_u=lambda^2/(8c),       dq/dF_u=lambda^2/(16pi c),
    dq/dJ_omega=lambda/(8c),     dq/dF_omega=lambda/(16pi c).
                                                                  (4)

All are nonzero. In the action-angle convention with angle of period2pi,
the return angle is2pi rot_R, so its derivative with respect to J_u is
2pi R dq/dJ_u, with central limit pi R lambda^2/(4c).
Changing orientation changes its sign and does not change nonvanishing.
Equation(4) states each normalization explicitly rather than importing an
unmatched numerical torsion factor from another section convention.

## Actual finite-ring transfer

Suppose0255 supplies an axisymmetric finite-R solution with psi_R=R phi_R,
F(psi)=lambda psi, and core convergence phi_R->Phi in C4 on a neighborhood
of a fixed small regular annulus. Here R tends through sufficiently large
values only to prove existence of at least one suitable finite template;
the selected template is then fixed before all preparation accuracy limits.
This is a conditional implication, not a claim that0255 already supplied
such a family. In fact an explicit sufficiently small C4 error at one
large R suffices.

Let C_R(s)={phi_R=s} be a regular closed meridional contour. The exact
poloidal velocity and angular toroidal velocity are

    |u_pol|=R |grad phi_R|/r,
    dot(varphi)=lambda R s/r^2,       r=R+x.                (5)

Over one poloidal traversal, the toroidal angle increases by

    Delta varphi_R(s)
      =lambda s integral_C_R(s) dl/[r |grad phi_R|].       (6)

The exact rotation number is rot_R=2pi/Delta varphi_R. Therefore

    rot_R/R
      =2pi/[lambda s integral_C_R(s)
                         (R/r) dl/|grad phi_R|].          (7)

The velocity-flux action of the disk D_R(s) enclosed by that contour is

    J_u,R(s)=(lambda/(2pi)) integral_D_R(s) (R/r) phi_R dx dz.
                                                                  (8)

For the limiting circle, (7) gives exactly (1) and (8) exactly (2).
The regular-level implicit-function theorem and differentiation of these
contour/area integrals show C1 convergence in any fixed regular action
interval of the normalized rotation number as a function of J_u,R. C4
convergence is more than enough for those derivatives away from the center.
The required denominators, gradient and action derivative stay separated
from zero on that fixed annulus. No estimate is taken through a shrinking
critical-level interval.

Choose this annulus small enough that dq/dJ_u is bounded away from zero
by (4). Equations(7)-(8) then preserve its sign for one sufficiently close
finite-R core. The axisymmetry gives exact invariant stream surfaces on
that actual ring, and the velocity-flux action is the invariant area of its
toroidal Poincare section. Its return map has the corresponding nonzero
action derivative. Multiplication by constant lambda supplies the vorticity
flux convention without changing the field lines or the conclusion.

Using the exact period quadrature avoids an erroneous exponential-in-R
trajectory-comparison estimate: the normalized quantity in (7) converges
directly by local contour geometry. The finite-ring return need not be
uniformly close in angle to the straight surrogate after a long orbit.

## Scope and consequence

The exact Bessel-core flux twist and the conditional finite-R persistence
are established as stated. No independently tuned twist parameter is needed
when the actual compact supplier retains this core in the specified norm.
On the open g=0 core, the equation is a uniformly elliptic linear
constant-curl stream equation away from the axis. Once the construction
gives convergence on a slightly larger core neighborhood, interior elliptic
regularity can upgrade lower-order convergence to the required derivative
control; the neighborhood and uniform coefficient bounds remain explicit.

This removes an unnecessary extra scalar tuning obligation. It does not
prove the degenerate outer inverse, the compact finite solution, global
stability, density, or the actual full-pressure response/action/current
join.0255 still earns those by constructing the finite stationary template;
0250/0257 continue the dynamical and physical-observation construction.
