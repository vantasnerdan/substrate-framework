# Smooth polygon continuation by a bordered Green equation

This is a constructive analytic route, independent of a patch-to-smooth
spectral inference. It produces some smooth monotone cores, not continuation
of every prescribed radial profile. No empirical coefficient is chosen.
The argument uses ordinary ODE continuation, compact Green operators, and the
implicit function theorem. Its critical kernel computation is given explicitly.

## 1. A radial core with a full bordered inverse

Begin with the radial solution of

```
Delta U + U_+^4 = 0,  U(0)=A>0, U'(0)=0.
```

While `U>0`, its radial mass `m(r)=integral_0^r s U(s)^4 ds` is positive
and increasing, and `U'=-m/r`. It crosses zero at a finite radius `R`: if it
remained positive, any positive lower bound on `m` after a fixed positive
radius would force `U` below zero logarithmically. At the crossing `U'(R)<0`.
Beyond `R`, continue by `U=-m(R) log(r/R)`. The circulation is
`Gamma=2 pi m(R)>0`.

The scaling family is

```
U_A(r)=A U_1(A^(3/2) r),  Gamma_A=A Gamma_1.
```

For the linearized operator `L=Delta+4 U_+^3`, its only regular radial
homogeneous solution, up to scale, is

```
Z0=(2/3) U+r U'.
```

It has logarithmic coefficient `-2 m(R)/3`, which is nonzero. Thus no
nonzero radial solution is bounded at infinity. The `m=1` bounded kernel
consists precisely of the two translations, and there is no bounded kernel
for `m>=2`: apply radial-core-gap.md with
`w=U_+^4`, `Omega=-U'/r`, and `a=4 U_+^3 Omega`. In that notation the
elliptic Birman--Schwinger operator is exactly `B_m`. Its `m>=2` norm is at
most `1/m`; its `m=1` unit eigenvalue is simple. Fourier decomposition
therefore proves that the full bounded kernel consists of translations.

Now choose once and for all a nondecreasing `C-infinity` step `gamma`, zero
on `t<=0`, one on `t>=1`, and strictly between zero and one in between. Set

```
F_eta(s)=s^4 gamma(s/eta),  eta>0.
```

This is a smooth nonnegative monotone function, flat at `s=0`, and
`F_eta -> s_+^4` in `C^1` on bounded intervals: the function and derivative
errors are bounded respectively by constants times `eta^4` and `eta^3`.
The radial initial-value problem `Delta U+F_eta(U)=0` converges with its
central-value derivative on bounded intervals to the preceding one. This can
be seen first at the regular singular origin from the Volterra equation

```
U(r)=A-integral_0^r t log(r/t) F_eta(U(t)) dt,
```

and then by ordinary ODE continuation away from zero. The transverse first
zero and its slope persist. Since `d Gamma_A/dA=Gamma_1>0` at `eta=0`,
the central value can be adjusted for small `eta` so that `Gamma` stays
exactly fixed. Hence the new compact core has the same circulation and a
`C-infinity` vorticity `F_eta(U_eta)`.

Here is the precise invertibility statement that survives this smoothing.
Take a fixed disk `B_L` containing the core with an annular collar where
`U_eta<0`. On the subspace `X=C^{1,alpha}(B_L)` even under `y2 -> -y2`,
use `N(x)=log|x|/(2 pi)` and the map

```
(v,mu,lambda) -> (
 v+N*(F_eta'(U_eta) v)-mu+lambda y1,
 integral F_eta'(U_eta) v,
 integral y1 F_eta'(U_eta) v).
```

It is a Fredholm operator of index zero: the Green part is compact on `X`,
and the borders are finite rank. At `eta=0` its kernel is zero. To see this,
pair its first component with `F_0'(U_0) partial_1 U_0`, the left
translation kernel. The constant pairs to zero and

```
integral y1 F_0'(U_0) partial_1 U_0 = -Gamma,
```

so `lambda=0`. The mass constraint cancels the logarithmic part of the
global extension `v=mu-N*(F_0'v)`. It is bounded, hence is a translation
by the kernel calculation above. Reflection leaves only `partial_1 U_0`;
the centroid constraint removes it by the same nonzero integral. Then
`v=mu=0`. Thus the border is invertible.

The potentials `F_eta'(U_eta)` converge uniformly on `B_L`, so the
bordered operators converge in norm. Their inverses persist for sufficiently
small positive `eta`. Fix one such `eta` from now on. It is an existence
parameter, not a fitted molecular coefficient. This closes the radial
circulation and translation difficulties for an actual `C-infinity` core.

For clarity, uniform potential convergence is enough for the asserted
operator norm: on a fixed disk, the logarithmic potential maps `L-infinity`
boundedly into `C^{1,alpha}` for every fixed `0<alpha<1`. Splitting its
gradient difference integral into distances below and above a displacement
`h` gives respectively bounds of order `h` and `h log(L/h)`, both bounded
by a constant times `h^alpha`. Thus multiplication by the potential difference
followed by the Green operator is small in `X -> X` norm. Compactness of
the unperturbed operator follows from its additional regularity for the
`C^{1,alpha}` input and the compact Holder embedding on the fixed disk.

## 2. Exact nonlinear polygon equation on fixed disks

Fix any `N>=2`, radius `d>0`, and the preceding core. Let `Q_j` rotate by
`2 pi j/N`, and let the `j`th core be

```
x=Q_j(d e1+epsilon y),
omega_epsilon(x)=epsilon^(-2) F_eta(U(y)),   y in B_L.
```

Outside the disjoint disks set the vorticity to zero. The negative collar
ensures this extension is smooth. The unknown profile `U` is even in `y2`;
impose its exact mass `Gamma` and centroid zero. Put
`D_j=d(e1-Q_j e1)`. On `B_L`, solve

```
U(y)-mu + integral N(y-z) F_eta(U(z)) dz
 + sum_{j=1}^{N-1} integral
   [N(D_j+epsilon(y-Q_j z))-N(D_j)] F_eta(U(z)) dz
 - Omega [epsilon d y1+epsilon^2 |y|^2/2] + lambda y1 = 0,

integral F_eta(U)=Gamma,
integral y1 F_eta(U)=0.
```

All integrals are over the fixed disk. External kernels are analytic in
`epsilon` for small `|epsilon|`; the self Green operator gains derivatives.
Smooth Nemytskii composition with `F_eta` makes the displayed map smooth
on a small neighborhood in `X x R^2`. The annular negative collar is an open
condition in this space. At `epsilon=0`, its derivative in `(U,mu,lambda)`
is precisely the invertible border just proved. The implicit function
theorem gives a unique nearby `(U,mu,lambda)(epsilon,Omega)`. Reflection
preserves the equation since it permutes `j` with `N-j`.

The auxiliary `lambda y1` is removed, not retained as an external force.
The exact polygon identity is

```
sum_{j=1}^{N-1} Gamma grad N(D_j)
  = Gamma (N-1)/(4 pi d) e1.
```

At first order in `epsilon` the external residual is therefore a multiple
of `y1`. The unique bordered solution has no first-order profile change,
and

```
lambda(epsilon,Omega)
  = epsilon [Omega d-Gamma(N-1)/(4 pi d)] + O(epsilon^2).
```

Because `lambda(0,Omega)=0` for every nearby `Omega`, the quotient
`lambda/epsilon` extends smoothly to zero (Hadamard's integral formula).
Its derivative with respect to `Omega` is `d`, nonzero. A second ordinary
implicit function theorem therefore gives an exact rate `Omega_epsilon`
near `Gamma(N-1)/(4 pi d^2)` for which `lambda=0`. On this branch
`U_epsilon=U_eta+O(epsilon^2)` in `C^{1,alpha}`; elliptic bootstrapping
gives smooth profiles and the same order in every fixed finite regularity
norm. No common-angle force or frozen background is supplied to this step.

## 3. Why this is an exact Euler equilibrium

Let `psi=N*omega_epsilon`, with the logarithmic normalization appropriate
to the nonzero total circulation. The displayed equation with `lambda=0`
is exactly

```
U(y)=constant-[psi(x)-Omega_epsilon |x|^2/2]
```

on each core, after absorbing the subtracted constants and the self
`Gamma log(epsilon)/(2 pi)` into `mu`. Therefore
`omega_epsilon` is a function of the relative streamfunction on every core.
It is zero in the collar and exterior. Consequently

```
(J grad psi-Omega_epsilon Jx).grad omega_epsilon = 0
```

everywhere. Rotation of this profile is a genuine smooth Euler solution on
the plane, consisting of `N` disjoint, compactly supported, equal-circulation
cores. The total kinetic energy has the usual whole-plane logarithmic
infrared divergence for nonzero total circulation; finite changes of its
renormalized interaction energy are well-defined. No finite total energy
claim is made.

For `N=6` this supplies the exact smooth finite-core base state needed for
the two-triangle molecular route. Its finite-dimensional affine action is
derived separately in finite-core-angle-action.md. Existence alone is not
asserted to prove an invariant optical eigenspace of the full Euler PDE.

## Analytic scope and next review

The positive object is an exact smooth polygon branch for sufficiently small
but nonzero smoothing and concentration parameters, with circulation fixed.
The proof is local and nonquantitative: it gives existence neighborhoods,
not numerical maximum core radii. Its key inputs are explicitly proved kernel
and border statements rather than a presumed spectral-persistence theorem.
Independent review should check the function-space map, bordered index and
left-kernel pairing, circulation normalization, and the division of the
centroid solvability condition by `epsilon`. The symbolic receipt checks
normalizations and scaling; it is not presented as a mechanized proof of
the Banach-space implicit function theorem.
