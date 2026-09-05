# A positive full-Euler force-free mode above the axial transport range

## 1. Actual endpoint and full radial equations

Take the smooth force-free column that is the straight limit of 0195's
global endpoint. In physical cylindrical coordinates,

    u=r Omega(r)e_theta+W(r)e_z,
    Z=2Omega+r Omega', Kappa²=2Omega Z,
    W W'=-r Omega Z=-r Kappa²/2.                   (1)

Z is positive in its fixed active disk, flat to zero at its outer radius
a_c, with Omega bounded below by Omega_min>0 there. Outside, Z=W=0 and
Omega=Gamma/(2pi r²). The source repair of 0195 preserves Omega'<=0:
it is confined to the outer flat tail where Z<2Omega has a strict margin.
In the fixed inner core Omega=Omega_c>0 and W is positive with a strict
maximum W0 at the axis. All derivatives are those of the actual endpoint;
W is not a fitted constant axial offset.

For pressure P divided by density and exp(i m theta+i k z-i omega t),
put s=omega-m Omega-k W and f=xi_r. Full momentum and incompressibility
give, without dropping axial shear,

    v_r=-i s f, v_theta=m P/(r s)-Z f, v_z=k P/s-W'f,
    f'=-(1/r+2m Omega/(r s))f+(m²/r²+k²)P/s²,
    P'=(s²-Kappa²)f+2m Omega P/(r s).               (2)

The explicit W' terms cancel from the two first-order equations, not
from the reconstructed velocity or KKS. Equations (2) agree with the
fully derived 0140 system, now without its small-shear assumption.

Choose m=0, k>0 and write omega=k c. With y=r f the complete radial
problem becomes

    [(c-W)² y'/r]'+[Kappa²-k²(c-W)²]y/r=0.         (3)

Regularity is y=O(r²) at zero. In the actual irrotational exterior it
is y=C r K1(k r), with the full exponentially decaying pressure/velocity,
not an imposed wall. This exterior is used throughout the proof.

## 2. A simple principal branch c>W0

Let H_k be the completion of compact regular radial functions under

    ||y||_k²=integral_0^infinity (|y'|²+k²|y|²)/r dr.

The weighted Hardy inequality integral |y|²/r³ <= integral |y'|²/r
follows by completing the square in integral |y'-y/r|²/r. In particular
the source quadratic form

    b[y]=integral Kappa² |y|²/r

is compact on H_k: its near-axis contribution is at most a constant
times delta²||y||_k², and ordinary compactness applies on each remaining
bounded annulus. There is no source at infinity.

For c>W0 define

    a_c,k[y]=integral (c-W)² (|y'|²+k²|y|²)/r.

This is coercive. Its positive scalar Sturm Green function and b define
a compact positivity-improving Birman--Schwinger operator. Its largest
eigenvalue Lambda_1(c,k) is simple; its eigenfunction is strictly positive.
The min-max ratio b/a_c,k shows strict decrease in c>W0, and Lambda_1
tends to zero as c tends to infinity.

At c=W0 use a nonzero regular trial y supported in a sufficiently small
inner disk of radius d. Since W0-W=O(r²), its a-form is O(d²) times its
b-form for a fixed scaled trial; Kappa²=4Omega_c² there. The k² term
is smaller still. Therefore a_W0,k[y]<b[y], for all sufficiently small
fixed k. By continuity Lambda_1(c,k)>1 for some c>W0. There is one and
only one principal c_*(k)>W0 with Lambda_1=1. Its eigenfunction has no
radial node, its full Evans root is simple and its velocity is an actual
smooth Euler normal mode. The vorticity perturbation is compactly
supported, while the actual exterior velocity is retained.

The branch stays finite as k decreases. Indeed

    b[y] <= ||Kappa²||_infinity a_c² ||y||_k²

by Hardy and compact active support. Thus
c_*(k) <= W0+a_c sqrt(||Kappa²||_infinity), a uniform finite bound.
The analogous k=0 weighted problem allows a constant y at infinity;
the same compact form and principal-root argument gives c_*(0)>W0.

The finite-core Robin matching to r K1(k r) has its usual k² log k
small-k expansion. Interior coefficients are smooth in (c,k²), and the
principal Evans derivative is nonzero. Hence

    c_*(k)=c_*(0)+O(k²|log k|),
    omega_*(k)²=c_*(0)² k²+O(k^4|log k|),           (4)

with the corresponding two differentiated bounds. The squared-frequency
curvature is positive for fixed sufficiently small k>0. The optical gap
at that FIXED carrier is nonzero; no carrier-independent gap at k=0 is
claimed. The axial phase/action measure is fixed when taking carrier
derivatives, not reset to a k-dependent period.

There is also an exposing exact failed alternative. At c=0, (1) and
integration by parts give

    a_0,k[y]-b[y]
      =integral W²[(y'-y/r)²/r+y²/r³+k²y²/r] >0   (5)

for nonzero admissible y. The negative-phase-speed attempt therefore
cannot use the same trial argument. The positive branch above W0 is a
different actual spectral construction, not a sign change by convention.

## 3. Full-poloidal isolation, including the translation modes

For small k, choose the fixed principal frequency inside

    k W0 < omega_*(k) < Omega_min/2,
    k W0 < Omega_min/2.                            (6)

The full transport bands are m Omega(r)+k W(r). All m=0 bands lie
below the first bound, positive m bands above Omega_min, and negative
m bands below -Omega_min/2. Thus the mode is outside EVERY transport
band, not only its axisymmetric sector.

Off these bands the full active-vorticity Euler operator has the
transport-plus-relatively-compact-pressure representation of 0195.
Only genuine eigenvalue coincidences need examination. The general
real-frequency identity derived directly from (2) is

    [r Re(f conjugate(P))]'
      =r (m²/r²+k²)|P|²/s²+r(s²-Kappa²)|f|².       (7)

The complete exterior potential gives Re(f conjugate(P))<0 at a_c.
For m>=3 or m<=-2, (6), Omega'<=0 and Kappa²<=4Omega² make the right
side strictly positive, excluding a mode at omega_*. Only m=-1,1,2
remain beyond the principal m=0 mode.

At k=omega=0 the m=2 determinant is nonzero: s²=4Omega²>=Kappa²,
so (7) and the exterior matching exclude a nonzero regular solution.
Its finite-core coefficients are regular at that point. It therefore
has no root in the small omega=O(k) neighborhood.

The m=1 zero mode is the exact Euclidean translation

    f0=1, P0=-r Omega².                            (8)

It is unique at zero by the regular first-order radial problem. Set
omega=k c+o(k). Differentiating (2) at zero gives delta s=c-W. The
interior derivative of the Wronskian r(f0 delta P-P0 delta f) is
IDENTICALLY zero, including the W contribution. At the origin that
Wronskian vanishes. The full exterior boundary condition is
P/f=-r s² at k=0; since W(a_c)=0, its first derivative gives

    r(delta P+r Omega² delta f)=2r² Omega c.

Thus c=0. Its derivative with respect to omega is nonzero, so the
translation root is simple. The same calculation for m=-1, or real
conjugation with k reversed, gives the same zero first derivative.
Their small-k frequencies are O(k²|log k|), separated from (6).
This uses the actual exterior W=0, not a locally chosen Galilean offset.

Within m=0 the principal root is simple and cannot coincide with another
radial eigenvalue. The resulting pole is therefore simple in the FULL
poloidal space. It has an honest isolated contour at each selected fixed
k>0, available for a later ring-kernel continuation.

## 4. Full KKS sign, including axial vorticity shear

The actual material displacement columns are

    xi=(f,i B,i C), B=-2Omega f/s, C=k P/s².

With the fixed axial phase measure L_z, the full Euler KKS coefficient is

    beta=2pi rho L_z integral f(Z B+W' C) r dr.      (9)

The W'C term is essential. Substituting y=r f, then using (3) and the
full exterior decay to integrate by parts, gives EXACTLY

    beta=-(2pi rho L_z/k)
           integral (c-W)(y'^2+k²y²)/r dr <0.      (10)

The cross term from differentiating 1/(c-W) cancels the full W' term;
discarding the latter would destroy this identity. The laboratory phase
Hamiltonian is -beta omega times the positive two-coordinate norm/2,
and is positive because omega=k c>0. The scalar mode action after
elimination has its actual positive inertia. This action is not an
assertion that the material tag's axial spin equals canonical momentum.

In fact m=0 enforces v_theta=-Z f and

    delta_material(r u_theta)=r v_theta+r Z f=0.    (11)

Each material particle retains its axial angular momentum on this
axisymmetric Kelvin perturbation. A literal material torsion/current
interpretation must respect (11). The next companion derives the
different, nonzero Eulerian section angular momentum and its actual
advective/pressure flux before making any observable identification.
