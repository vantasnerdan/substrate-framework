# Actual flux-action twist and persistence of the closed tube

This child of0145 treats the explicit constant-lambda field, not an
arbitrary local force-free ansatz. It constructs one unknotted invariant
solid torus with an elliptic core and robust boundary. No arbitrary-knot
or optical-mode theorem is inferred. The parent's registered0145
boundary was validated before this calculation.

## 1. Exact finite-radius section and normalization

Take lambda>0, kappa,k>0 with kappa²+k²=lambda², and

```
psi=A r J1(kappa r) cos(kz),
u=(-psi_z/r, lambda psi/r, psi_r/r).
```

Let j be any sufficiently large positive zero of J0, and put R=j/kappa.
J1(j)!=0 by uniqueness of the Bessel ODE. Choose
`A=U/[lambda J1(j)]`, U>0, so u_theta(R,0)=U. The circle r=R,z=0
is an exact periodic core. In the local coordinates x=r-R,z its
normalized flux function is

```
h_R(x,z)=psi/(UR/lambda)=f_R(x)cos(kz),
f_R(x)=(R+x)J1(kappa(R+x))/[R J1(j)].
```

Exactly, f_R(0)=1, f_R'(0)=0 and

```
f_R''-f_R'/(R+x)+kappa² f_R=0.                (1)
```

On every fixed compact x-interval away from r=0, ordinary parameter
ODE estimates give
`f_R=cos(kappa x)+O(1/(kappa R))` in every fixed scaled C^m norm.
This does not require a numerical Bessel-zero asymptotic. For example
the first correction is

```
f_R=cos(kappa x)
 +[x cos(kappa x)-sin(kappa x)/kappa]/(2R)+O(R^-2).
```

At kappa=k the limiting normalized stream function is
`h=cos(kappa x)cos(kappa z)`. Fix once a closed regular level interval
`0<h_-<h_+<1`, contained in the central positive cell. Its level curves
are smooth simple closed curves with gradient uniformly separated
from zero. For sufficiently large R the same is true for h_R and a
slightly larger level annulus. Rotating these disks about the axis
gives actual embedded unknotted solid tori; u_theta>0 throughout them.

## 2. Actual poloidal transit, toroidal advance and flux action

The finite-radius poloidal velocity is

```
(u_r,u_z)=UR/[lambda(R+x)] J grad h_R.
```

Here J(x,z)=(-z,x). Since h_R has a maximum at the core, the local
flow is clockwise in the (x,z) chart. On the level h_R=h the actual
toroidal speed is `theta_dot=UR h/(R+x)²`. Therefore the toroidal
angle advanced during ONE poloidal circuit is exactly

```
Delta_theta_R(h)=lambda h integral_(h_R=h)
                    dl/[(R+x)|grad h_R|].          (2)
```

The toroidal-section return map, theta increasing by2pi, has positive
poloidal cycle count
`rho_R(h)=2pi/Delta_theta_R(h)`. We orient its poloidal angle along
the flow. Reversing that angular convention changes the signed rotation
number and the twist sign, not twist nonvanishing.

The ACTUAL invariant section measure is u_theta dx dz, not unweighted
area and not r times that measure. Define its outward-increasing action

```
I_R(h)=1/(2pi) integral_(h_R>=h) u_theta dx dz.
```

Coarea gives the exact finite-R identity

```
I_R'(h)=-UR h/(2pi) integral dl/[(R+x)|grad h_R|]
       =-UR Delta_theta_R(h)/(2pi lambda),
I_R'(h) rho_R(h)=-UR/lambda.                  (3)
```

This identity keeps the nonuniform transit and the physical flux in
the same action. Analytic action-angle coordinates on each regular
annulus conjugate the return map to this rotation. In flow-oriented
angle coordinates the positive section form can be written
`d(angle) wedge dI`; either ordering convention is legitimate provided
the rotation sign is kept consistent.

## 3. Exact Taylor--Green limit and strictly nonzero twist

For kappa=k, set X=kappa x,Z=kappa z and
`Omega=U kappa²/lambda`. The limiting material equations are

```
X_dot=Omega cos X sin Z,
Z_dot=-Omega sin X cos Z,
h=cos X cos Z.
```

In one quadrant, `X_dot=Omega sqrt(cos²X-h²)`. Substituting
`sin X=sqrt(1-h²)sin t` in the quadrant transit integral gives

```
T(h)=4 K(1-h²)/Omega.                         (4)
```

K denotes the complete elliptic integral with PARAMETER m, not modulus.
The limiting toroidal advance is `U h T(h)/R`. Consequently

```
F(h):=lim rho_R(h)/R=pi kappa²/[2lambda h K(1-h²)],
I_0'(h)=-2U h K(1-h²)/(pi kappa²).            (5)
```

Let m=1-h² and let E(m) be the complete elliptic integral of the second
kind. Differentiating their defining integrals, or using the exact
elliptic derivative identities, gives

```
d/dh[h K(1-h²)]=(K(m)-E(m))/m,
K(m)-E(m)=m integral_0^(pi/2)
                 sin²t/sqrt(1-m sin²t)dt >0.       (6)
```

Thus F'(h)<0 and I_0'(h)<0, and the actual flux-action twist is

```
dF/dI = pi² kappa^4/(4lambda U)
        (K-E)/[m h³ K³] >0.                  (7)
```

Its center limit is `kappa^4/(2lambda U)`; at kappa=k this is
`lambda³/(8U)`. Equation(7), rather than an unweighted radial
frequency, is the nonzero torsion supplier.

All contour integrals and their first two h derivatives in (2)--(3)
depend smoothly, indeed analytically, on the nearby analytic defining
function on the fixed regular annulus. One explicit way to see this
is to identify nearby contours using the bounded normal vector field
`grad h/|grad h|²`, then differentiate the integrands on the same
parameter circle. The gradient lower bound and a fixed larger annulus
bound every denominator. Equation(1) therefore gives

```
||rho_R/R-F||_C²(h_-,h_+) + ||I_R'-I_0'||_C¹
    <= C/(kappa R).                          (8)
```

The exact limiting twist in (7) has a strictly positive minimum on
this compact interval. Choose finite R large enough for (8) to preserve
half that margin. Then

```
(1/R) d rho_R/d I_R = dF/dI +O(1/(kappa R)) >0.  (9)
```

The return ANGLE is2pi rho_R, so its action derivative is also nonzero.
No unspecified small error is compared to a vanishing endpoint twist;
the annulus and its positive margin were fixed first.

## 4. Arithmetic-free selection of the core and boundary

The exact core transverse linearization has frequency U kappa k/lambda,
and its toroidal speed is U/R. Hence its return rotation number is

```
rho_core=R kappa k/lambda=j k/lambda.          (10)
```

At the single choice kappa=k, the number j/sqrt(2) is NOT declared
irrational or nonresonant. No such Bessel-zero arithmetic is needed.

For the original single-mode family, keep lambda and the large integer
index fixed, and vary k/lambda in an arbitrarily small open interval
about1/sqrt(2). The core and regular annulus persist, as does the strict
twist margin. Equation(10) is strictly monotone in this continuous
parameter. Choose a value avoiding half-integers, or choose a Diophantine
core rotation if desired. The core is then strictly elliptic, with
multipliers different from +1 and -1. Its transverse jet is as close
to circular as desired. With this field now fixed, (9) makes the boundary
rotation a local diffeomorphism of the level h. Diophantine numbers have
full measure in every open interval, so some h in(h_-,h_+) supplies a
Diophantine invariant boundary. This chooses an actual flux level, not
a Floquet winding.

### Optional exact-circularity repair within global constant-lambda fields

If exactly circular core linearization is required simultaneously,
there is an explicit small two-mode extension, still globally constant
lambda. Fix a small b>0 with b<lambda and consider

```
psi=r J1(kappa r)cos(sqrt(lambda²-kappa²)z)
   + t(2/b)r J1(br)cos(sqrt(lambda²-b²)z).         (11)
```

At z=0 impose a critical circle r=R and equal transverse Hessian
entries. Using the exact Bessel equation, these conditions are

```
C1=kappa J0(kappa R)+2t J0(bR)=0,
C2=(lambda²-2kappa²)J1(kappa R)
   +t(2/b)(lambda²-2b²)J1(bR)=0.             (12)
```

At the original circular point t=0,kappa=lambda/sqrt(2),R=j/kappa,
the determinant of partial(C1,C2)/partial(kappa,t) tends as b->0 to

```
2 kappa J1(j)(4-j²) !=0.                     (13)
```

All positive zeros j of J0 are larger than2: for 0<=x<=2 the alternating
series, with z=x²/4<=1, gives
`J0(x)>=1-z+z²/4-z³/36>0`. Also J1(j)!=0 as noted above. Thus a finite
small b has nonzero determinant, with no assumption about rationality
or fractional parts of any Bessel zero. The ordinary analytic IFT
solves (12) for kappa,t as R varies on an open interval about the
chosen large radius. Both modes have nonzero radial wavenumber; their
global velocity and derivative bounds are finite once b is fixed.

Normalize the overall amplitude so u_theta at the core equals U>0.
At a critical point the constant-lambda equation and isotropic Hessian
give `Hess_(r,z) psi=-(lambda² psi/2)I`. Therefore the EXACT circular
core rotation number is now `lambda R/2`, strictly varying with R.
Choose a nearby finite R with a nonresonant, even Diophantine, core
number. The strict annular twist survives this small field change,
and select a Diophantine boundary level as before. This optional
extension preserves exact circularity without unproved arithmetic.
It is declared as a two-mode extension, not attributed to the original
single-mode formula at kappa=k.

## 5. Robust invariant tube and actual nearby flux measure

For the selected field, choose a slightly larger positive-flow solid
torus around the boundary. Its analytic toroidal-section return map
preserves the positive physical measure u_theta dr dz. The boundary
is analytically conjugate to its selected Diophantine rotation, and
(9) gives nonzero normal torsion in action-angle coordinates.

The applicable primary result is EPS, arXiv:1210.6271v2, Theorem7.6
and the measure identification in equations(7.28)--(7.32), already
archived under `../../sources/1210.6271.pdf`. The required theorem and
measure argument were read directly for this child. The theorem applies
to analytic area-preserving return maps close in a sufficiently high
C^k norm, with nonzero boundary torsion and a Diophantine boundary.

A nearby analytic divergence-free field has a nearby positive section
measure, generally NOT the same u_theta dr dz. On a slightly larger
disk, the difference of the two analytic area forms is exact. Moser's
near-identity measure identification conjugates the new return map to
one preserving the old measure. The section and return map remain
defined because the positive toroidal component and the finite return
time persist at this fixed geometry. Apply Theorem7.6 after that
identification. It supplies a nearby invariant torus with the same
Diophantine boundary rotation. Suspension gives an embedded invariant
solid-torus boundary in the actual nearby field.

The core return fixed point persists separately by IFT because its
multiplier is not +1; the strict ellipticity margin keeps the nearby
area-preserving fixed point elliptic and different from ±1. A chosen
Diophantine core is not claimed to remain Diophantine under EVERY
perturbation; that infinite nonresonance property is not open. The
robust conclusion is the nondegenerate elliptic periodic core and KAM
boundary just stated. The solid torus remains unknotted under the
small ambient deformation.

All approximation/KAM thresholds are chosen AFTER the finite geometry,
twist margin and nonresonance choices. The parent's global/periodic
constant-lambda approximation may then be taken below these finite
thresholds; this child does not assume an approximation-uniform radius
or a hidden Bessel-zero Diophantine bound.

Verdict: established as stated for the actual flux-action twist,
arithmetic-free parameter selection and robust closed unknotted tube.
The exact-circular two-mode extension is an additional constructive
option, not needed if arbitrarily near-circular single-mode geometry
suffices. No numerical spectrum, small-ratio numerical design or
arbitrary-knot conclusion is involved.
