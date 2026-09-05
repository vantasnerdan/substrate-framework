# Exact curved bending and its full Euler residual

This is new work in 0213. The smooth global ring is the actual 0211/0195
Green family, in its translating frame. Fix the core profile, taper and
one finite toroidal harmonic first; then increase the ring radius R.
All constants below may depend on that fixed profile and harmonic.
The vorticity is supported in a tubular neighborhood of fixed radius a,
has uniformly bounded fixed spatial derivatives, and the decaying lab
velocity is v. The stationary-frame velocity is u=v−U e_z, where
U=O(log R/R). No exterior velocity is cut off.

## 1. Repairing a first-order shear in the bending lift

Write x=r−R, and let q(phi), Z(phi) be smooth mean-zero periodic radial
and vertical displacements. Take their mean-zero antiderivatives Q'=q,
H'=Z and set A=q'+Q. The following physical cylindrical components are
exactly divergence-free:

    xi_r = q−z Z/R + [x²(q''+q)/2+zx(Z''+Z)]/(Rr),
    xi_phi = −Q−x A/R−z Z'/R,
    xi_z = r Z/R.                                             (1)

An explicit vector potential is

    A_r=0,
    r A_phi=−R q z+(r³−R³)Z/(3R),
    A_z=x Q+x²(q'+Q)/(2R)+zx Z'/R−z H.                        (2)

Direct cylindrical curl gives (1). Choose a smooth cutoff equal to one
on an open neighborhood of the entire vorticity support, with compact
support away from the cylindrical axis. Then curl(cutoff A) is a global
smooth compact divergence-free generator. Its induced Euler velocity
is still the complete-space V xi=P(xi cross omega); only the generator,
not that velocity, is compactly localized.

The less complete lift xi=(q−zZ/R,−Q+zH/R,rZ/R) has shear
E_rphi=(q'+Q)/(2R) at the centerline. Its displacement and divergence
are correct, but this shear gives the wrong residual order. Formula
(1) adds the physical cross-section rotation and its exact divergence
return. It is a method repair, not a change of the measured centerline.

Let E=sym Dxi in an orthonormal cylindrical frame. E vanishes exactly
on the centerline, and on the fixed core

    ||E||_{C^j} <= C_j/R².                                   (3)

For example its R²-leading coefficient, with x,z fixed, is

    [[x(q+q'')+z(Z+Z''), 0, x(Z+Z'')/2],
     [0, −x(q+q'')−z(Z+Z''), 0],
     [x(Z+Z'')/2, 0, 0]].

For q=a cos phi,Z=0 it is exactly a e_x on the whole source. For
q=0,Z=b cos phi it is exactly (−b/R)e_y cross position. These exact
translation and rotation germs give zero strain, not just a small
approximation. In particular (3) does not manufacture a force on a
rigid mode.

## 2. Exact Kelvin and full-pressure identities

Use [a,b]=(a dot grad)b−(b dot grad)a. Define

    V xi=P(xi cross omega),
    A_E xi=V xi−[u,xi].                                      (4)

V uses the full R³ Leray projection with decaying induced velocity.
For a compact generator the independent Cartesian identity is

    A_E xi=−2P[(sym Dxi)u].                                  (5)

Indeed xi cross curl u−[u,xi]
=grad(xi dot u)−(Dxi+Dxi^T)u. This is also the actual fixed-Kelvin
reconstruction sign; replacing A_E by its negative changes evolution.

Steady Euler gives [u,omega]=0. Curl Vxi=[omega,xi], so Jacobi and the
linearized vorticity equation give L V=V A_E. The zero harmonic velocity
condition at infinity fixes the curl ambiguity. Thus for every smooth
prepared generator history,

    (partial_t−L)Vxi=V(partial_t xi−A_E xi).                  (6)

This is the actual coadjoint forcing class, not an independently
postulated body force. Its physical dimensions are consistent:
xi is a displacement, A_E xi a velocity, and V A_E xi an acceleration.

At a rigid translation germ, A_E xi=0 on the active source. At a
rotation germ xi=J position, Vxi=Jv−(J position dot grad)v, while
[u,xi]=Ju−(J position dot grad)u. Consequently

    A_E xi=U J e_z                                           (7)

on the source. The actual translating-frame rotation has the neutral
translation Jordan response already derived in 0206; setting (7) to
zero would lose it. Compact completion does not change Vxi or its
source germ, because omega vanishes where that completion differs.

## 3. Actual Biot–Savart commutator estimate

Let K be the complete Biot–Savart matrix kernel, v(x)=integral K(x−y)
omega(y)dy. Integration by parts in the compact vorticity source gives
the exact source-point identity

    A_E xi(x) = −Dxi(x)u_infinity
      + integral [K(x−y)Dxi(y)−Dxi(x)K(x−y)]omega(y)dy
      + integral DK(x−y)[xi(x)−xi(y)]omega(y)dy.              (8)

Here u_infinity=−U e_z. This formula uses only the generator germ on
the source and hence is independent of its arbitrary compact return.
It retains the induced pressure and the complete distant ring.

For fixed x in a slightly enlarged core, subtract the rigid affine
field xi(x)+O_x(y−x), where O_x=skew Dxi(x). The resulting rigid
integrand is exactly zero by rotational covariance of K. On source
points at distance d<R/4, the explicit formula (1) and its first two
derivatives give

    |Dxi(y)−O_x| <= C(a+d)/R²,
    |xi(y)−xi(x)−O_x(y−x)| <= C(ad+d²)/R².                  (9)

The bounds follow in overlapping tubular charts from the first-order
Taylor formula, (3), and ||D²xi||<=C/R² there. They do not use a
small-strain bound on an artificial cutoff outside the source.

Since K=O(d^-2) and DK=O(d^-3), the near integrand in (8) is bounded by
C |omega(y)|(a+d)/(R²d²). The volume of a source shell of radius d is
O(d² dd) for d<a and O(a² dd) for a<d<R/4. Integrating gives
C(1+log(R/a))/R². For d>=R/4, |xi|<=C, |Dxi|<=C/R and the whole source
volume is O(R); the contribution is O(R^-2). Finally
|Dxi u_infinity|<=C |U|/R=O(log R/R²). Therefore

    ||A_E xi||_{L-infinity(core)} <= C(1+log R)/R².          (10)

The same estimate holds in any fixed core C^j norm if the fixed-profile
Green family and q,Z have the corresponding extra derivatives. To see
the near-kernel derivative bound without differentiating an absolute
d^-3 majorant, split at a fixed core radius, integrate a derivative
onto the smooth source in that near integral, and subtract the same
local rigid Taylor polynomial. All near terms contain E or D²xi,
each O(R^-2); the distant terms can be differentiated directly.
The fixed profile supplies uniform local source derivatives. There is
no new logarithmic loss from this operation.

In particular an actual bending phase with frequency
O(log R/R²) has eta=partial_t xi−A_E xi of that same C^1 order. The
measured core displacement is O(1). This is a full Euler residual
estimate, not yet a claim that its accumulated error is small on the
reciprocal-frequency time interval.

## 4. The new seed really licenses the weighted physical response

The 0211 cutoff obeys h(t)=exp(1−1/t)(1+O(t)) at t downarrow0.
Its defining integral therefore gives, by one endpoint integration by
parts (or l'Hopital with positive functions),

    H_delta(phi) = (e/delta) phi² exp(−delta/phi)(1+O(phi)),
    H_delta H_delta''/(H_delta')² -> 1.                     (11)

At the actual simple radial outer flux zero, phi' tends to a nonzero
constant and phi'' stays bounded. Since H_delta/H_delta' is O(phi²),

    w w_rr/w_r²
      = H_delta H_delta''/(H_delta')²
        +(H_delta/H_delta') phi''/(phi')² -> 1,             (12)

where w=lambda H_delta and the reflected axial velocity is W=−w.
Thus the physical-velocity Hardy step in 0206 applies to THIS new
seed, not merely to its old exponential example. This conclusion uses
the explicit tail, not smooth compact support alone.

For the straight axisymmetric component of eta, the exact 0206 forced
Sturm transformation consequently yields the source estimate

    |ell_eta(Y)| <= C_eta sqrt(Q_c(Y)),
    C_eta=||sqrt(s) w' eta_theta/Omega||_2
       +sqrt(2)||s^(3/2)(eta_z'−ik eta_s)||_2.              (13)

The actual source in (10), projected onto this poloidal component,
has C_eta=O(log R/R²) at k=n/R. Its reflected straight pressure/velocity
response has the corresponding uniform weighted bound even when
c=omega/k tends to zero from the positive sector. The c=0 exterior
zero-velocity quotient from 0206 is retained. Formula (13) is a direct
licensed use of the computed residual; it is not yet an estimate for
the entire curved fast-channel inverse.

## 5. Next exact reaction object

The remaining reaction is the full curved complement Schur form.
The n=1 Euclidean subspace is removed by its exact inherited KKS
projection, including (7). For n>=2 the cross-section translations
are the bending coordinates themselves; they cannot be discarded as
global rigid symmetries. The fast complement includes all other
poloidal harmonics and the axisymmetric weighted threshold channel.

The response-specific norm (13) is stronger than a crude 1/omega
bound, but importing it into the curved inverse requires controlling
curvature-induced feedback in that SAME weighted source space. In
particular, a small unweighted velocity operator error alone does not
control the flat-tail generator norm. The next calculation is the
flux-coordinate quadratic-form pullback and its off-diagonal pressure
rows. No positive pole or acoustic-window conclusion is inferred from
(10) alone. The actual n=2 global quadrupole and its zero linear global
spin selection from 0206 also remain distinct physical observations.
