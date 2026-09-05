# The complete mean-metric threshold response

This is an actual-action preconditioner for the full curved Schur
problem. Restrict contravariant velocities to poloidal order zero in
the volume/action-angle coordinates of `actual-operator.md`. Their
kinetic metric is the poloidal mean of the ACTUAL physical metric.
This restriction does not yet eliminate the other velocity harmonics
or the actual exterior harmonic return. Those two corrections are
the remaining full Schur calculation, not assumed absent here.

## 1. Exact metric and covariant background identities

Reflection across the meridional symmetry section makes the mean
I-theta and I-y metric rows zero. Write the remaining positive metric

    g_bar=[[A,0,0],[0,B,C],[0,C,D]],
    E=D−C²/B>0.

The background is u=(0,Omega,W), omega=f u, with actual mean metric
covectors satisfying

    (B Omega+C W)'=f W,
    (C Omega+D W)'=−f Omega,
    C Omega+D W=G(H),  H'=−Omega.                         (1)

Primes in this document mean I derivatives. Define

    T=B' Omega+C' W,
    J=C' Omega+D' W−C T/B,
    d=c−W,      temporal frequency omega_mode=k c.

Differentiating (1), rather than supposing a Euclidean radial
profile, gives the exact cancellation

    J+f(Omega+C W/B)=−E W'.                               (2)

## 2. All forced Kelvin, pressure and velocity rows

Let eta=(eta_I,eta_theta,eta_y) be the solenoidal coadjoint source at
poloidal order zero. The ACTUAL unit-volume constraint is
eta_I'+ik eta_y=0. Define eta_theta_cov=B eta_theta+C eta_y.
A full coadjoint generator and pressure are

    Xi_I=Y,
    Xi_y=iY'/k,
    Xi_theta=−i[(T Y−eta_theta_cov)/(B d)+C Y'/B]/k,
    p=i[J Y+d E Y'−E eta_y]/k.                           (3)

They solve the complete forced Kelvin/pressure system precisely when

    (d² E Y')'+[d J'+f W T/B−k² A d²]Y
      =f W eta_theta+(f W C/B)eta_y
         +d(E eta_y)'−ik A d eta_I.                     (4)

The corresponding contravariant PHYSICAL velocities are

    v_I=−ik d Y−eta_I,
    v_y=(dY)'−eta_y,
    v_theta=−f W Y/B−C v_y/B.                            (5)

Equations (1)-(2) imply, exactly,

    g_bar v=Xi cross omega−grad_k p,
    v=−ikd Xi−(0,Omega'Y,W'Y)−eta.

Thus (4) retains angular-axial metric coupling, pressure, the forced
particular terms and all generator components. The checked flat
limit is A=1/(2I), B=2I, C=0,D=1; it returns the original reflected
Sturm equation after I=s²/2. The pressure check uses both identities
in (1); omitting them was the recorded first-verifier repair.

## 3. Why the metric cross row cannot destroy the flat-tail norm

Put W=G(H) a(I), where the actual period average in (10) of
`actual-operator.md` gives smooth positive a=1+O(delta_R). Here
delta_R=R^-2 times a fixed finite logarithmic factor records the
first nonzero averaged curvature remainder; no pointwise relative
comparison of shifted flat functions is used. The exact covariant
row in (1) gives

    C=G(H) c_bar(I),
    c_bar=(1−D a)/Omega=O(delta_R).                       (6)

This factor holds exactly wherever the active flux coordinate is
defined. In particular C vanishes with the SAME flat G as W.
The remaining metric estimates on a fixed outer annulus are

    A=(1+O(delta_R))/(2I), B=2I(1+O(delta_R)),
    D=1+O(delta_R), E=1+O(delta_R),
    J=O(delta_R)(|w'|+w),         w=−W>=0,               (7)

with fixed derivative versions where needed. Formula (7) for J
follows directly from its definition, (6), and G=−w/a. In the
central disk these are regular polar metric estimates; the angular
degeneracy of B is the ordinary axis regularity, not a new boundary.

## 4. Relative coercivity, including the apparently dangerous J'

Use the reference quadratic form with the ACTUAL nonnegative w:

    Q_ref(Y)=integral [(c+w)² |Y'|²
       +k²(c+w)² |Y|²/(2I)+w w' |Y|²/I] dI.            (8)

Integration by parts gives the exact positive decomposition

    Q_ref=integral w²[|Y'−Y/(2I)|²+|Y|²/(4I²)
                       +k²|Y|²/(2I)] dI
       +integral (c²+2cw)[|Y'|²+k²|Y|²/(2I)] dI.         (9)

The regular-axis and true decaying-column return make the boundary
term zero. At c=0 the exterior zero-velocity quotient is retained.
For c>0 this is the full positive Hardy norm. The actual tail w is
G times the smooth positive a, so w w''/(w')² tends to one as in
0213. The tail Hardy estimate gives

    integral (w')² |Y|² <= C Q_ref(Y).                   (10)

The form controls integral (c+w)²|Y'|² and the corresponding Hardy
norm of (c+w)Y. Constants depend on the fixed profile, not on small c.

The form from (4) is

    Q_bar=integral [d²E |Y'|²
             +(k² A d²−d J'−f W T/B)|Y|²] dI.

Its J' row is controlled after its actual integration by parts:

    −integral d J' |Y|²
       =integral J[w'|Y|²+2d Re(Y_bar Y')].             (11)

By (7), (9)-(10) and Cauchy–Schwarz, (11) is bounded in absolute
value by C delta_R Q_ref. Bounding J' pointwise before this step
would invent a false loss at the flat edge.

Finally f=−G'/Omega and W=−w give

    f W T/B = w G' B'/B−w² G' C'/(Omega B).

Since G=−w/a and C=G c_bar, its difference from −w w'/I is a sum
of terms with O(delta_R) coefficients multiplying w w', w²,
w²(w')² and w³w'. All are bounded as quadratic-form multipliers by
(9)-(10), with the regular central weights supplied by B~2I.
The kinetic coefficient differences in (7) have the same bound.
Consequently

    |Q_bar−Q_ref| <= C delta_R Q_ref.                    (12)

This is an actual relative form estimate for the mean-metric
preconditioner, not an unweighted profile perturbation. Taking R
large makes its full positive-sector coercivity strict. A fixed
narrow complex sector follows by the same sesquilinear estimate and
the real-part version of (9).

The source in (4) is bounded in the dual Q_ref norm by fixed
source Sobolev norms: fW has the factor wG', C has (6), and dY has
the Hardy control. Thus the complete physical response (5) and the
coadjoint force from (3) retain the uniform-velocity and
k^-1 polylog(1/c) graph bounds of 0213. The particular pressure and
angular-axial terms have not been discarded to obtain them.

## 5. Remaining full physical return

The mean metric uses the inherited positive kinetic action of the
zero-poloidal velocities. It is a comparison block, not a claim that
the physical pressure stays in that harmonic. To turn (12) into the
full ring inverse, retain the nonzero-poloidal metric/pressure blocks
and the difference between the true exterior harmonic return and
the straight decaying return used in (9). Their coefficients and
moment cancellations are the next calculation. Positivity of this
comparison form does not itself certify an actual Euler pole.
