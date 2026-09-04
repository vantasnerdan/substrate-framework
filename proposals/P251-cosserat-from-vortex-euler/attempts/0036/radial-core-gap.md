# Exact smooth radial-core separation

This is an independent analytic lemma, not a consequence of patch convergence.
It concerns planar Euler on the plane and the tangent to a fixed vorticity
rearrangement class. It establishes the fast-core estimate needed by a future
coupled-core reduction. It does not assert that the nonlinear six-core steady
state or its slow invariant subspace has already been constructed.

## Object, signs, and space

Let `w(r)` be a nonnegative, nonincreasing, nonzero, smooth compactly supported
radial vorticity, with support in `r <= R`. Set

```
a(r) = -w'(r)/r >= 0,
M(r) = integral_0^r s w(s) ds,
Omega(r) = M(r)/r^2,
Gamma = 2 pi M(R).
```

The proof below also applies to compactly supported, radially regular `C^2`
profiles: its integrations and bounded kernels need no infinite-order
smoothness. In particular the `U_+^4` seed used in
smooth-polygon-construction.md is `C^3` at its transverse edge and satisfies
these hypotheses before its subsequent C-infinity smoothing. This extension
does not use a discontinuous patch limit.

The velocity is `u = J grad psi`, `J(x,y)=(-y,x)`, and `Delta psi=w`.
Thus the base angular velocity is positive `Omega`. On the core,
`Omega_min=Gamma/(2 pi R^2) <= Omega <= Omega_max=w(0)/2`:
`Omega'=(w-2 Omega)/r <= 0`, since `w` is nonincreasing.

For each angular Fourier number `m>=1`, use the Hilbert space
`L^2({a>0}, r dr)`, with perturbation `eta=sqrt(a) f exp(i m theta)`.
The closure in this norm of smooth isovortical perturbations is the space at
issue. Axisymmetric changes of the radial profile change Casimirs and are not
tangent perturbations of this radial rearrangement orbit. They are not declared
fast or discarded from an unrestricted Euler initial-value problem.

The positive inverse of `-Delta_m`, regular at zero and decaying at infinity,
has kernel, relative to measure `s ds`,

```
G_m(r,s) = (1/(2m)) (min(r,s)/max(r,s))^m.
```

Linearizing `partial_t w + u.grad w=0` gives the self-adjoint frequency
representation

```
partial_t f = -i m H_m f,
H_m = Omega - sqrt(a) G_m sqrt(a).
```

In particular `H_m` has units of inverse time, not mechanical energy. It is the
symmetrized linear generator on the stated vorticity tangent space. Converting
it to a physical action requires the Euler Poisson form; no mass is inferred
from this formula alone.

## Translation identity and Schur estimate

The radial Poisson equation gives the exact identity

```
G_1(a r) = M/r = Omega r.
```

Indeed `-Delta_1(M/r)=-w'=a r`; the solution has both required endpoint
conditions. Equivalently, splitting the Green integral at `s=r` and integrating
by parts gives the same expression. Therefore
`f0=sqrt(a) r` is the translation kernel of `H_1`.

Define

```
B_m = Omega^(-1/2) sqrt(a) G_m sqrt(a) Omega^(-1/2),
h = sqrt(Omega) sqrt(a) r.
```

Then `B_1 h=h`. The kernels are nonnegative and symmetric. Since
`G_m <= G_1/m` pointwise, the weighted Schur test with positive weight `h`
gives `||B_m|| <= 1/m`. Consequently, for every `m>=2`,

```
<f,H_m f> >= (1-1/m) integral Omega |f|^2 r dr,
m H_m >= (m-1) Omega_min.
```

This includes continuous as well as point spectrum. No finite-basis spectrum
approximation is used.

## Explicit gap after the m=1 translation kernel

The Doob transform `P=h^(-1) B_1 h` is a reversible Markov operator on the
probability measure

```
d pi(r) = Omega(r) a(r) r^3 dr / Z,
Z = integral_0^R Omega a r^3 dr = integral_0^R r w(r)^2 dr > 0.
```

The last identity follows by integrating `-M w' = r w^2 - (M w)'` and
using `M(0)=0`, `w(R)=0`. Its transition density with respect to `pi` obeys

```
dP(r,s)/d pi(s)
  = Z G_1(r,s)/(Omega(r) Omega(s) r s)
  = Z/(2 Omega(r) Omega(s) max(r,s)^2)
  >= delta,
delta = Z/(2 Omega_max^2 R^2) > 0.
```

The identity holds almost everywhere on `{a>0}`; zeros of `a` are outside the
weighted state space. The lower bound extends by limits where needed.
Normalization gives `delta<=1`. Hence
`P=delta Pi+(1-delta)Q`, with `Pi` projection on constants and `Q` another
reversible Markov contraction. This proves

```
<g,(I-B_1)g> >= delta dist(g,span(h))^2.
```

The gauge matters. The stronger weighted estimate

```
<f,H_1 f> >= delta integral Omega |f|^2 r dr
```

holds under `integral Omega f sqrt(a) r (r dr)=0`. This is not ordinary
orthogonality to `f0`. For the ordinary invariant complement `f perpendicular
f0` in `L^2(r dr)`, one instead uses the projection comparison:

```
<f,H_1 f>
 >= delta inf_c ||sqrt(Omega)(f-c f0)||^2
 >= delta Omega_min inf_c ||f-c f0||^2
 = delta Omega_min ||f||^2.
```

Thus `ker H_1=span(f0)` and the ordinary orthogonal complement has the
strict spectral gap `delta Omega_min`. The comparison does not require
equating the two orthogonality conditions.

## A C-infinity compact core with an explicit, nonnumerical bound

For

```
w(r) = c exp(-1/(1-r^2/R^2))  (r<R),  0 (r>=R),
c>0,
```

the preceding assumptions hold. `Omega_max=c/(2e)`. On
`0<=r<=R/sqrt(2)`, `w>=c exp(-2)`, so

```
Z >= c^2 exp(-4) R^2/4,
delta >= exp(-2)/2,
gap_fast >= exp(-2) Gamma/(4 pi R^2).
```

All nontranslation, nonaxisymmetric isovortical core frequencies are therefore
separated from zero by a constant of order `Gamma/R^2`, with an explicit
positive lower bound. The bound is independent of a thin transition-layer
limit, quadrature, roundoff, and any point-vortex comparator.

## What the lemma activates next

For six separated cores of radius `epsilon` and fixed molecular separation
`d`, the candidate fast scale is `Gamma/epsilon^2`, whereas the centroid
scale is `Gamma/d^2`. The estimate above controls the radial reference
operator. To transfer it to a coupled exact equilibrium one must construct
that equilibrium in coordinates that retain its isovortical space, control
the full operator perturbation (including transport), and perform the slow
Schur reduction with the physical symplectic form. Mere uniform velocity or
Holder convergence does not supply those operator estimates.

Two executable representations remain: (1) nested level sets preserving all
superlevel areas, whose nonlinear composition estimates still require proof;
(2) a local smooth vorticity-versus-relative-streamfunction nonlinearity,
using fixed disjoint cutoff disks to avoid differentiating a moving flat edge.
Neither representation has yet earned a coupled existence theorem here.
