# Rotational-coherence removal and complete-fluid removal are different

This applies to0102's full canonical conditional action, including its
remaining finite derivative kinetic terms and actual observation maps.
Its microscopic base is one stationary Euler field law. Three kinds of
geometry serve different purposes: core-angle pairs, angle-gradient
attachments, and the independent STF strain attachments. Their thinning
parameters are declared as part of the ensemble, not silently identified.

## 1. Remove rotational coherence while retaining positive translation shear

Fix the background covariance, finite good-patch geometry and STF population.
Independently thin the core-angle AND angle-gradient populations with one
parameter a in[0,1]. A retained point keeps all its precomputed moments and
reactions; no coefficient is refitted after thinning. Campbell averaging
then gives exactly

    j(a)=a j0, kappa(a)=a kappa0,
    C_spin(a)=a C_spin,0, M_spin,grad(a)=a M_spin,grad,0.

Every term involving Phi or its derivatives carries the same factor a,
including higher spatial terms of the exact neighbor-difference action.
The associated mean/current filters vanish with a as well. The negative
macro gradient inertia -j(a)|beta_dot|^2/2 is retained before this limit.
The original displacement density rho is independent of a because the
background and ambient fluid have not been deleted.

The fixed and linear shear corrections of all thinned reactions are
uniformly bounded for a in[0,1]. Their positivity/inverse bounds are the
per-cell ones already established before averaging; thinning does not
invert a vanishing averaged reaction matrix. Choose the fixed STF attachment
amplitude large enough that its positive quadratic Schur contribution
dominates this bounded family, including the negative bare covariance
shear. Thus the actual shear mu(a) stays positive on the full interval,
and mu(0)>0 is derived from the surviving STF action, not borrowed from
the old affine vorticity-pushforward formula.

At a=0 remove the absent Phi coordinate in the UNREDUCED action. To leading
second spatial order the surviving incompressible equation is

    rho U_tt=mu(0) Delta U-grad p, div U=0.

Its explicitly retained strain-gradient inertia supplies higher dispersive
corrections, if that fuller physical-field action is being used. In the
same derivative-normalized convention as0102 the displayed equation is
the Navier--Cauchy sector. No longitudinal compressional Euler mode is
created. The P242 correspondence is this conditional translational sector,
not equality to an unrelated five-scalar filament coefficient.

Taking the a>0 spin equation, dividing by j(a), and then setting a=0 would
leave a spurious oscillator whose whole action and physical observation
weight vanished. The actual limit is coordinate deletion, as the already
reviewed0063 unreduced-matrix argument establishes.

## 2. Remove the whole vortical background and all attached structures

Take u_b=b u_1 with b>0, preserving the same nonzero curl eigenvalue.
Every sample remains stationary Euler and p_b=b^2 p_1. Its covariance
and the bare material shear scale as b^2. On fixed compact generator
geometry, Vop and KKS scale as b, while H_orbit and material K scale as
b^2. Therefore B^2/P, where that internal chart remains valid, need not
vanish as b decreases: a ratio of two small coefficients is not a
license to retain a rotor when its vorticity-angle chart disappears.

For the complete-fluid removal family also thin EVERY attached population,
including STF attachments, with a parameter c->0. Keep finite geometry
and bounded structural amplitudes, or state an explicit bound on their
population-weighted forms. Then all internal action weights and their
observation rows vanish with c; all elastic terms vanish with b^2 c or
with the bare b^2 covariance. At b=c=0 the defining action contains only
the original incompressible fluid kinetic term. Delete all absent angular
and attachment coordinates before any Schur inverse involving them.

The resulting linear continuum equation is exactly the linearized Euler
fluid sector rho U_tt=-grad p, div U=0. The microscopic equation throughout
the construction was Euler; this statement is its unstructured limit at
the declared linear continuum scope, not a derivation of nonlinear Euler
from a linearized constitutive formula.

These two limits preserve the positive claims without equating removal of
coherent rotation, removal of strain-bearing structures, and removal of
the energetic Euler background. They do not set a surviving quadratic
fluctuation to zero merely because its coherent first moment is zero.
