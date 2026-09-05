# Actual periodic image forcing and the steady continuation interface

Candidate A concerns the compact-vorticity closed Euclidean ring, not
candidate C's noncontractible periodic tube. This source computes the
actual first image fields and separates the remaining steady inverse.
No nonlinear stationary superposition or unproved implicit theorem is
asserted.

## Exact local periodic Green expansion

Let G_P solve -Delta G_P=delta-P^-3 on the cubic torus of side P with
zero mean. In a ball of radius less than P/2,

    G_P(x)=1/(4pi|x|)+P^-1 g(x/P),

where g is smooth, even and cubic-invariant and Delta g=1. Its exact
Taylor structure is

    g(z)=g0+|z|^2/6+c4 H4(z)+O(|z|^6),
    H4=z1^4+z2^4+z3^4-(3/5)|z|^4.

The coefficient c4 is the actual lattice Green coefficient, not an
adjustable modulus. Every fixed derivative of the remainder is bounded
on a smaller ball. The quadratic coefficient follows from Delta g=1
and cubic symmetry; H4 is harmonic. This uses the full periodic inverse,
including its zero-mean condition.

For compact smooth divergence-free vorticity omega supported in that
ball, int omega=0 by int div(x_i omega)=0. Write the actual impulse
I=(rho/2)int y cross omega(y)dy. The complete image velocity is

    u_im(x)=int grad[P^-1 g((x-y)/P)] cross omega(y)dy.

The quadratic term therefore gives the exact leading UNIFORM velocity

    u_im,2=-2 I/(3rho P^3).

It can be included in the actual translating frame; it is not a shape
force. A point-vorticity monopole estimate would miss this cancellation.
After this term, for a thin ring of radius R and circulation Gamma,

    |u_im-u_im,2| <= C |Gamma| R^4/P^5,
    |grad u_im| <= C |Gamma| R^3/P^5

on a fixed core neighborhood when P=dR and d is sufficiently large.
These follow by subtracting the constant-in-y kernel before using
int omega=0, so the bound uses int |y||omega|=O(|Gamma|R^2), not a
fictitious nonzero monopole. Finite-core moment corrections are retained
by their actual vorticity integral; no filament is needed for the bounds.

## A sharper exposed circular-ring term

For the circular filament only as an explicit leading geometric test,
omega(dy)=Gamma R(-sin theta,cos theta,0)dtheta and
y=R(cos theta,sin theta,0). The quartic image term on z=0 has

    (u_im,4)_z=-(12pi c4 Gamma/5P^5)
                       [R^2(x^2+y^2)+R^4/2],
    (u_im,4)_x=(u_im,4)_y=0.

On the actual reference circle this is a constant axial drift. One
derivation integrates the polynomial exactly; equivalently planar
Green's formula reduces it to -Gamma int_disk Delta_perp H4(x-y).
The angular-four polynomial piece is planar harmonic and cancels in
that integral. Hence the first quartic lattice anisotropy does not by
itself force a fourfold deformation of a circular filament. Higher
terms and finite-core/profile corrections remain genuine sources.

## The precise next stationary construction

The true periodized ring velocity solves div u=0 and curl u=omega, but
its original vorticity generally has a nonzero Euler transport residual
[u_im,omega]. A constant drift can be absorbed in the frame. The
remaining source is smooth in the core and has the scale above; it is
not zero just because supports of distinct images are disjoint.

The reviewed0220 operator is a full positive-frequency response inverse
on a specified curved source graph. Its contour estimate is not, by
itself, an inverse of the steady operator at frequency zero. A nonlinear
steady continuation requires deriving its zero-frequency source-range
estimate, retained translation/rotation and circulation-profile kernels,
and sufficient tame regularity for iteration. In particular the flat
vorticity tail makes a naive unweighted elliptic inverse inappropriate.

This is the method-repair target: use the actual vanishing factors in
[u_im,omega], quotient the geometric neutral rows, and test the weighted
zero-frequency operator rather than import the positive optical contour.
The image computation gives an actual compatible forcing and its leading
frame cancellation. It does not claim the stationary density array yet.
