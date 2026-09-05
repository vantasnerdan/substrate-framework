# A genuine Euler Floquet packet on the EPS periodic core

## Source reconciliation and the constructive ordering

The archived primary [EPS, arXiv:1210.6271v2](https://arxiv.org/abs/1210.6271)
does not put a Lundquist oscillator inside its knotted tube. Its §7
chooses `λ=ε³` for convenience, and the leading thin-tube field is
harmonic. Proposition 7.12 obtains elliptic core monodromy from
`exp(±i∫τ dα)` with total torsion outside `πZ`. That geometric statement
alone is not a theorem about linearized Euler perturbations.

There is, however, a source-supported constructive repair. Theorem 6.8
permits **every** nonzero sufficiently small `|λ|`, at a fixed sufficiently
thin tube, with the explicit estimate

```
||vα-hα||Ck + ε||vy-hy||Ck ≤ Ck,Λ ε|λ|.
```

Section 7 expressly permits any nonzero `λ=O(ε³)`, and Remark 7.11
states that the persistence argument works for small arbitrary `λ`,
including the harmonic limit, with error terms `O(ε²+|λ|)` and
`O(ε³+|λ|)`. Therefore use this order:

1. Choose a source-admissible analytic knot with total torsion outside
   `πZ`, and perform the source's small generic adjustment for Diophantine
   boundary rotation and nonzero normal torsion.
2. Fix one sufficiently thin **positive** `ε` with a strict elliptic core
   trace margin and robust boundary torus for its harmonic field.
3. Only then choose nonzero `|λ|` smaller than both the source persistence
   tolerance and the finite gyro bound derived below.
4. Apply the source's Theorem 8.3 global approximation with that SAME
   eigenvalue, small enough in the required local `Ck` norm to retain both
   the invariant tube and all strict estimates below.

No unknown `exp(C/ε)` is assumed to be dominated by `ε³`. At the frozen
tube every needed constant is finite, and `λ` remains a free smaller
parameter. This constructs a decaying global EPS Beltrami field, not a
remote added Lundquist field. Its core has period `T` and particle
monodromy with eigenvalues `1,e^{±iθp}`, `θp∉πZ`.

## 1. Actual Euler amplitude and the periodic Kelvin-compatible covector

Along the actual periodic core `x(t)`, put `A(t)=Du(x(t))` and let
`Mdot=A M`, `M(0)=I`. Choose the unique real left eigenvector `k0` of
`M(T)` at eigenvalue one, normalized by `k0·u(x0)=1`. This pairing is
nonzero because the eigenvalue one is simple and its right eigenvector
is the nonvanishing flow direction. Define

```
k(t)=M(t)^(-T) k0,       kdot=-A^T k,       k(T)=k0.
```

The actual short-wave Euler amplitude equations are

```
bdot=-A b+2k(k·A b)/|k|²,           k·b=0.             (1)
```

They are derived again by the pressure recursion below, rather than
imported as a spectral conclusion. These are the geometric-optics
equations used in the primary Euler literature, e.g.
[Friedlander–Strauss–Vishik, Theorem 4.1](https://content.ems.press/assets/public/full-texts/serials/aihpc/14/2/4398291/online/10.1016-s0294-1449-97-80144-8.pdf);
no instability theorem from that paper is used.

Since `ω=λu`, the exact Cauchy equation on the stationary trajectory is
`ωdot=Aω`, and

```
s=k·ω=λ k·u=λ ≠0
```

is conserved. For the material amplitude

```
a=k×b/s,       k·a=0,
b=-(s/|k|²)k×a,
adot=A a-(s/|k|²)k×a.                                (2)
```

The signs and the coefficient follow by differentiating `k×b`, using
`tr A=0` and `(A-A^T)v=ω×v`. In particular (2) is the actual leading
Lin equation `adot-Aa=b`, not a selected positive-Hamiltonian pair.
The map between the two transverse amplitudes is invertible through
every finite interval because `s` never vanishes. Its constants may grow
as `λ→0`; `λ` is fixed before choosing the packet scale.

## 2. Why the Euler return is elliptic on this selected EPS field

Let `C=A-(s/|k|²)[k×]` and `Fdot=C F`, `F(0)=I`. Then

```
tr C=0,       kdot=-C^T k,
det F(t)=1,   k(t)^T F(t)=k0^T.
```

Thus `F(T)` preserves the actual plane `k0⊥`, and its transverse return
`FE` has determinant exactly one. The quotient multiplier is one;
`tr FE=tr F(T)-1`. This exact determinant, not a numerical unit-modulus
test, is what turns a strict trace estimate into ellipticity.

The gyro perturbation has operator norm
`|s|/|k|≤|ω|≤|λ|U*`. Fix finite bounds `|u|≤U*`, `||A||≤Γ*`, `T≤T*`
in a sufficiently small neighborhood of the fixed harmonic template and
its continued core. The particle return has a uniform strict margin
`m*=2-|tr ME|>0` after reducing that neighborhood if necessary. Duhamel
gives, in physical Euclidean coordinates,

```
||F(T)-M(T)|| ≤ |λ|U*T* exp[(Γ*+|λ|U*)T*].            (3)
```

Choose the nonzero eigenvalue so that `|λ|U*≤1` and

```
3|λ|U*T* exp[(Γ*+1)T*] < m*/2.                       (4)
```

All choices are compatible with Theorem 6.8 and Remark 7.11 at the fixed
ε. The global approximation is then chosen sufficiently accurate to
preserve (4). Hence `|tr FE|<2`, `det FE=1`, and its eigenvalues are
`e^{±iθE}`, with `0<θE<π`. This proves a genuinely Euler-compatible
elliptic amplitude sector **on the same EPS periodic core**.

In a periodic orthonormal frame for `k(t)⊥`, real Floquet theory gives
`a(t)=P(t)R(νt)aF`, with `P(t+T)=P(t)`, invertible, and
`ν=θE/T>0` after a fixed real conjugacy. The physical signal has the
periodic geometric modulation `P`; it is not a pure rigid lab-frame
oscillator. An optical time `2π/ν` is finite. No full Euler spectral
stability or invariant finite-dimensional PDE subspace follows.

## 3. Compact high-order Euler packets, with exact initial Kelvin data

Take a small material neighborhood of `x0` strictly inside the invariant
tube. Give it a smooth compact envelope, and choose an initial phase with
gradient `k0` at its center. Transport phase and envelope by the actual
flow: `Dφ=0`, `k=∇φ`, `D=∂t+u·∇`. Shrink the neighborhood so that `|k|`
and `|k·ω|` stay bounded away from zero throughout the prescribed finite
optical time. Its support stays inside the material invariant tube.

The following explicit recursion constructs the full pressure response
to arbitrary finite order. With `z=iN`, write formal compact amplitudes

```
v=e^{iNφ} Σ z^-j bj,
π=e^{iNφ} Σ z^(-j-1) πj.
```

Set `b_-1=0`, `π_-1=0`, `cj=-div b_(j-1)`. Enforce

```
k·bj=cj,
πj=[-2k·A bj-k·∇π_(j-1)-D cj]/|k|²,
D bj=-A bj-kπj-∇π_(j-1).                            (5)
```

These are linear transport equations with known smooth lower-order
forcing. Differentiating the constraint gives (5) exactly. For `j=0`,
it is (1). Each step has finite smooth coefficients on the compact swept
support. It is a constructive finite recursion, not a numerical closure.

To make the truncated velocity **exactly** solenoidal, set `d_-1=0` and

```
dj=-k×[bj-curl d_(j-1)]/|k|²,
VN=curl[e^{iNφ} Σ_(j=0)^m z^(-j-1) dj].              (6)
```

The required transverse compatibility follows from `curl k=0` and the
constraint in (5). Formula (6) agrees with the formal velocity through
order `m` and has one explicit extra compact term of order `N^(-m-1)`.
Applying the full Euler operator and the pressure sum through `πm`
leaves a compact residual bounded in L² by `Cm N^(-m-1)`. Leray
projection retains its whole-space exterior pressure reaction and cannot
increase that norm. Eikonal transport cancels the potentially large
`N Dφ` terms; every remaining coefficient is a fixed smooth profile.

For exact Kelvin preparation start with the compact solenoidal initial
displacement

```
ΞN=curl[e^{iNφ0} a_potential/(iN)],
a_potential=-k×a0/|k|².
```

Its leading amplitude is `a0`. Define the actual initial velocity
`vN(0)=P(ΞN×ω)`. The initial coefficients in (5) are not freely reset:
write `ΞN×ω=e^{iNφ0} Σ z^-j fj`, and recursively set

```
q_-1=0, b_-1=0,
qj=[k·fj-k·∇q_(j-1)+div b_(j-1)]/|k|²,
bj=fj-kqj-∇q_(j-1).                                 (7)
```

The difference between the initial field and (6) is an explicit gradient
plus an L² remainder `O(N^(-m-1))`. The whole-space Leray contraction
therefore proves the same initial error without discarding its exterior
tail or invoking a local inverse curl. Take real or imaginary parts to
obtain real perturbations. This is exactly prepared Kelvin data for the
actual EPS field, not just a leading-order isovortical symbol.

The full Euler energy estimate proves

```
sup_(0≤t≤2π/ν) ||vN-VN||2 ≤ Cfinite N^(-m-1).         (8)
```

Reconstruct the exact Lin displacement from initial value `ΞN`. Its
amplitudes obey `D aj=A aj+bj`, with their actual prepared initial
coefficients. The same transport estimate gives (8) for displacement
and material velocity after a fixed finite change of constant. Thus the
entire complement is controlled, not frozen. All constants can depend
on the fixed tube, λ, envelope, finite order and interval; `N` is then
chosen finite and large enough.

## 4. Actual core angle and full material moments

Choose the phase at the parcel center to be `π/4` and take imaginary
parts. A small initial ellipsoid of size `ell/N`, transported with the
actual base and perturbed flows, lies inside the invariant EPS tube.
Let its covariance be `Q(t)/N²`, where to leading small-parcel order
`Qdot=AQ+QA^T`, and its mass is `MN=O(N^-3)`.

The leading physical material-vorticity tilt is

```
qphysical=N cos(π/4) (s/|ω|) P_(ω⊥) a,
Φphysical=(ω/|ω|)×qphysical.                          (9)
```

The map from `a∈k⊥` to this transverse angle is invertible because
`k·ω≠0`. It has the same nonzero Floquet quasifrequency with periodic
geometry modulation. The exact observation uses a positive transported
weight and includes `∫wt[curl v+ξ·∇ω]`; it is not a Fourier coefficient
or a fixed-point angle imposed on a remote cage.

The actual centered material definitions give the leading moment rows

```
U=sin(π/4) a,
δP=MN sin(π/4) adot,
δS=(MN/N) cos(π/4)
       [a×A(Qk)+(Qk)×adot].                          (10)
```

They include the deformation of the tag and the ambient background
velocity gradient. In particular the spin is not inferred from
`∫r×v` alone. Unlike the special straight rotational comparison in 0109,
the general knot supplies this time-dependent tensor map, not the same
scalar `MR²/5` inertia.

Nonzero spin is constructive, not presumed. At time zero rotate coordinates
so `k` points along `e3` and a chosen nonzero amplitude along `e1`, with
their harmless positive scales restored afterward. Then
`b=-(ω·e3)e2` and `A21-A12=ω·e3≠0`. The linear map acting on `q=Qk` is

```
Lq=e1×Aq-q×[-(Ae1+b)],
L=[e1×]A-[Ae1+b ×].
```

Its entries satisfy `L31+2L13=A21-A12≠0`; hence `L` is not zero. At
least one of `q=e3,e3+e1,e3+e2` has `Lq≠0`. Each is realized by the
explicit positive covariance

```
Q = [[1+q1², q1 q2, q1],
     [q1 q2, 1+q2², q2],
     [q1, q2, 1]],     Qe3=(q1,q2,1).
```

Its Schur complement is `I2`, so it is strictly positive. Scale by a
sufficiently small positive `ell²` and realize it as the covariance of
a solid ellipsoid. This gives a nonzero initial full physical spin row.
Continuity and a finite bound on its derivative give a positive explicit
lower bound for its squared time integral over the optical interval.
Momentum also has positive squared time integral: a nonzero elliptic
Floquet amplitude cannot have `adot` identically zero, because its return
has no eigenvalue one. Individual components may vanish at particular
times; a constant scalar spin/rate law is not claimed.

The small ellipsoid gives a relative Taylor error `O(ell²)` in (9)–(10),
plus `O(1/N)` from spatially varying coefficients. Unknown Euler errors
are integrated over the **actual** parcel. Its volume sensitivity is
`N^(3/2)`; its spin lever arm is `O(1/N)`; the smooth positive angle
weight has gradient norm `O(N^(5/2))`. Consequently (8) gives relative
errors `O(N^(1/2-m))` for both moment and angle rows against their
nonzero leading scales. Choose `m≥3`, then fix `ell>0` small enough
for the strict spin/angle margins, and finally choose finite `N` large
enough. The actual material momentum and spin have positive finite-time
squared norms, with quantitative error, inside this SAME EPS tube.

The exact moving-domain balances remain
`Pdot=-∮p n` and `Sdot=-∮p r×n`. The full pressure in (5)–(8) and its
external reaction supply those exchanges. The surrounding tube is not
replaced by a slipping wall. Kelvin preparation and actual Euler
transport preserve the initially tangent vorticity of the material tube;
the perturbed tube need not remain stationary.

## Temporal observation and the action boundary

The frequency `ν=θE/T` is the **comoving amplitude** Floquet
quasifrequency. The actual Euler perturbation retains its entire phase
`exp(iNφ(t,x))`, with `Dφ=0`, as well as periodic geometry modulation.
A fixed spatial probe also sees carrier advection. If a separate global
longitudinal mode is sewn around the closed tube with integer winding,
that carrier contributes the circulation harmonics `n 2π/T` and their
Floquet sidebands. The local compact packet constructed here does not
itself impose such a global winding or discard its advected phase.

In particular `ν²` is not identified with the optical gap of conditional
C-CST-009. An autonomous relative rotor would require its actual Euler
KKS/Hessian restriction, Krein sign, and material observation map, with
any moving-frame generator or Routh subtraction performed on that same
action. Elliptic Floquet multipliers and nonzero material spin alone do
not establish a positive canonical rotor inertia. The physical moment
map in (10) is retained precisely to expose this distinction.

## Verdict and remaining parent construction

Established at finite-time geometric-optics scope: a source-constructible
EPS invariant tube with an actual Euler elliptic Floquet core amplitude,
arbitrarily accurate compact Kelvin-prepared Euler packets realizing its
physical material angle, and nonzero actual parcel momentum and spin.
The original near-Lundquist insertion route was replaced by a same-core
spectral construction with a quantitative gyro margin.

This is not yet a nonlinear invariant finite-dimensional Euler manifold,
a scalar isotropic inertia on an entire knot, or the full coupled
finite-wave-number continuum. The exact Galilean extension in the other
artifact supplies independent centroid initial data at macro `k=0`.
The actual slow-envelope reaction and the moment/tensor joining are
the next constructive consumers; a positive cage Hessian alone does not
replace them.
