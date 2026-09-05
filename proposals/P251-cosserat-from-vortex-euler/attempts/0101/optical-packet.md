# An actual Kelvin-prepared finite-time optical tilt packet

## 1. Positive statement and observable

For every prescribed relative error `ε>0` there is a smooth stationary
Beltrami Euler field on R³, containing an EPS knotted invariant tube, and a
two-dimensional family of exactly isovortical perturbations whose **actual
linearized Euler evolution** has a physical, positive-weight finite-core
vorticity-axis tilt within `ε` of a circular oscillator throughout one
period. The packet lies in a distant, finite straight flow-box of this same
field, not in the knotted torus. The carrier and localization parameters
are finite and chosen from the estimates below. No invariant finite-mode
Euler ansatz, favorable Schur truncation, MHD dynamics, spectral stability,
or uncomputed small-memory premise is used.

The tilt can be observed with a materially transported positive weight,
including its first-order moving-weight correction (§7). In an axial
Galilean frame its frequency is `Ω=Uλ/2`; relative to the local corotating
material axes it is `2Ω`. It is not the large axial Doppler frequency `Uk`.
The estimate is finite-time, not an exact global eigenfunction or a
long-time invariant rotor. In particular it does not yet identify the
angle or inertia of the EPS knotted tube or complete the parent continuum.

## 2. Exact base and finite parameters

Fix `λ,U>0`. The globally smooth Lundquist field

```
uL = U[J1(λr)eθ + J0(λr)ez],       curl uL = λuL,
pL = -ρ|uL|²/2
```

has a finite smooth axial core. In the frame moving with `Uez`, compare it
with `uR=Ω ez×x`, where `Ω=Uλ/2`. Write `δ=λa≤1/2`. On `r≤2a`, the
convergent alternating Bessel series and their differentiated recurrence
give, with `s=λr≤1`,

```
|J0(s)-1| ≤ s²/4,   |J1(s)-s/2| ≤ s³/16,
|(J1(s)-s/2)'| ≤ 5s²/16.
```

The derivative estimate follows from `J1'=J0-J1/s`. Thus for
`d=uL-Uez-uR`, conservative bounds are

```
||d||∞ ≤ 2Uδ²,          ||Dd||∞ ≤ 3Uλδ.                 (1)
```

Choose a fixed smooth radial `F:R²→[0,1]`, equal to one on the unit disc
and zero outside radius two. Choose the explicit odd compact function

```
H(ζ) = ζ exp[-1/(1-ζ²)]  (|ζ|<1),    H=0 otherwise;
h=H'.
```

All profile integrals below are fixed finite constants, independent of
`a,k,δ,U`, and all denominators are strictly positive. Write
`IF=∫F²`, `IG=∫|∇F|²`, `IF2=∫|D²F|²`, and
`IH,Ih,Ih1,Ih2=∫H²,∫h²,∫h'²,∫h''²` on the appropriate dimensions.

## 3. Exact packet, return, and pressure residual

Let `Jb=ez×b`, `b(t)=R(-Ωt)b0`, with horizontal `b0`, and put `Fa=F(x⊥/a)`.
Define the compact potential and its actual divergence-free velocity:

```
Aapp = k^-1 Jb Fa h(kz),
vapp = curl Aapp
     = -b Fa h'(kz) + ez k^-1 (b·∇Fa)h(kz).             (2)
```

The vertical return is essential. Direct differentiation gives the exact
uniform-rotation linearized Euler residual *before* pressure projection:

```
∂t vapp + (uR·∇)vapp + (vapp·∇)uR
 = -2Ω k^-1 (Jb·∇Fa)h(kz) ez.                         (3)
```

Thus horizontal columns solve the desired rotation exactly; the compact
solenoidal return supplies a computed, nonzero residual. Apply the full
R³ Leray projector `P` to (3); its exterior pressure response is retained,
and its L² operator norm is one. It is not replaced by a local projector.

The time-independent norm `N=||vapp||2` satisfies exactly

```
N² = |b0|²[a² IF Ih1/k + IG Ih/(2k³)],
||Rrotation||2² = 2Ω² |b0|² IG Ih/k³.
```

Consequently, for `ka≥1`,

```
||Rrotation||2/N ≤ Ω Crot/(ka),
Crot = sqrt(2 IG Ih/(IF Ih1)),
||∇vapp||2 ≤ CD k N,                                  (4)
CD = sqrt(IG/IF)+sqrt(Ih2/Ih1)
     +sqrt(IF2 Ih/(IF Ih1))+sqrt(IG/IF).
```

The four terms in `CD` bound horizontal-profile derivatives, axial
derivatives, transverse return derivatives, and axial return derivatives
separately. These deliberately conservative triangle bounds avoid relying
on unproved directional Hessian cancellations. Equations (1)–(4) imply

```
||RL||2/N ≤ Ω[Crot/(ka) + 4CD δ(ka) + 6δ].             (5)
```

Set `ka=δ^-1/2`, so `a=δ/λ` and `k=λδ^-3/2`. Every quantity is finite for
fixed positive `δ`, and the right side is at most
`Ω(Crot+4CD+6)sqrt(δ)`. This is an analytic parameter selection, not a
numerical wavelength or box-size convergence claim.

## 4. Exactly Kelvin-prepared initial data

At time zero let

```
Wh = -b0 Fa h(kz)/k,
Wz = (b0·∇Fa)H(kz)/k²,
ξ0 = curl W/(2Ω).                                     (6)
```

Both `W` and `ξ0` are smooth and compact, `div ξ0=0`, and `∂z W=vapp(0)`.
For the exact uniform-rotation vorticity `ωR=2Ωez`,

```
curl(ξ0×ωR) = 2Ω ∂z ξ0 = curl vapp(0),
P(ξ0×ωR)=vapp(0).                                     (7)
```

The final equality uses the global decaying divergence/curl solution, not
an equality before Leray projection. Explicitly

```
ξ0,h = [-Jb0 Fa h' - J(D²Fa)b0 H/k²]/(2Ω),
ξ0,z = (Jb0·∇Fa)h/(2Ωk).
```

For `ka≥1`, the same separated norms give

```
||ξ0||2 ≤ Cξ N/Ω,
Cξ = [1+sqrt(IF2 IH/(IF Ih1))+sqrt(IG Ih/(2IF Ih1))]/2.
```

For the *actual* background vorticity `ω`, use initial velocity
`v(0)=P(ξ0×ω)`. This is exactly on its Kelvin leaf. Equivalently, push
`ω` by the volume-preserving flow generated by `ξ0` and invert curl; the
derivative is this initial velocity. No arbitrary initial Kelvin momentum
is reset to the oscillator value. Actual linearized Euler then preserves
this preparation. Its Lin displacement solves
`ξt+[u,ξ]=v` with initial value (6).

## 5. The same stationary field contains EPS topology

Let `uE` be the campaign's already sourced EPS field of curl eigenvalue
`λ`, with robust knotted invariant torus and the stated derivative decay.
Choose `U>0` small enough in the required finite `Cm` norm to preserve that
torus after adding any translate of `uL`. The sum

```
u(x)=uE(x)+uL(x-c),    curl u=λu,
p=-ρ|u|²/2                                             (8)
```

is one exact stationary Euler field. The nondecaying Lundquist component
is intentional: (8) is smooth and bounded with bounded derivatives, not
a claim of finite total background energy on R³. Perturbation energies
and all packet norms here are finite. EPS persistence is a local robust
topology import; no approximate superposition of nonlinear solutions is
being used, since the common Beltrami eigenvalue makes (8) exact.

In coordinates `y=x-c-Ut ez`, the EPS correction is
`uE(c+Ut ez+y)`. During `T=2π/Ω` its center travels the finite distance
`UT=4π/λ`. Choose `c` sufficiently far from the EPS core that, throughout
this finite swept flow-box,

```
||uE||∞ ≤ Uδ²,      ||DuE||∞ ≤ Uλδ,
||ωE||∞ ≤ Ωδ.                                         (9)
```

Smaller bounds used for the transported observation below can be imposed
at the same finite choice of `c`. All relevant EPS derivatives decay, so
these are achieved conditions, not hypotheses imposed on an unknown
Euler solution. Equations (5),(9) yield

```
||Ractual(t)||2/N ≤ Ω C sqrtδ,
C = Crot+6CD+8.
```

Also `||ω-ωR||∞≤4Ωδ` on the support of (6), using the Bessel bounds and
(9), whence
`||v(0)-vapp(0)||2/N ≤ 4Cξ δ`.
The initial Leray tail outside this support is included in that inequality.

## 6. Full Euler evolution and quantitative optical-period control

Let `v` be the unrestricted linearized Euler solution about (8). Bounded
smooth coefficients give the usual finite-time energy evolution in L²;
pressure is obtained with the same whole-space Leray projection. In the
moving coordinates choose a finite global bound
`Γ≥||DuE||∞+||DuL||∞`, independent of `δ,c`. Subtract the actual equations
for `v` and (2), integrate against their difference, and use incompressible
advection and pressure orthogonality. This proves, not assumes,

```
sup(0≤t≤T) ||v-vapp||2/N
 ≤ exp(ΓT)[4Cξ δ + 2π C sqrtδ] =: Eδ.                 (10)
```

The residual is the actual omitted infinite-dimensional Euler dynamics;
the estimate controls all of it through one optical period. For any
target `e>0`, for example choose

```
sqrtδ ≤ min(1/sqrt2, e exp(-ΓT)/(4Cξ+2π C)).           (11)
```

The possibly large `exp(ΓT)` is retained. This is an existence/control
result, not a claim that the resulting scales are numerically efficient.
For a prescribed finite number of periods replace `T` by that finite
time and repeat the same choice. No uniform infinite-time conclusion
follows. For each fixed set of finite scales, this linear solution is the
derivative of actual smooth Euler solutions for sufficiently small
perturbation amplitude over the fixed time interval; no amplitude-uniform
high-frequency nonlinear lifespan is claimed.

## 7. A positive-weight physical tilt, including material transport

Choose smooth nonnegative normalized `w⊥` supported in radius `1/2`,
radial, and `wz` supported in `|ζ|<1/10`. Set
`w0=a^-2 k w⊥(x⊥/a)wz(kz)`. These are positive sampling weights, not a
Fourier demodulator. On their support `Fa=1`, and

```
curl(vapp)⊥ = -k Jb h''(kz),
ch = -∫wz h'' > 0,
qapp(t) = k ch Jb(t)/(2Ω),
qapp,t = -Ω Jqapp,         qapp,tt + Ω²qapp=0.         (12)
```

Positivity is explicit: in the interior

```
h''(ζ) = -2 exp[-1/(1-ζ²)]
          (3ζ⁸+24ζ⁶-26ζ⁴+3)/(1-ζ²)^6;
3ζ⁸+24ζ⁶-26ζ⁴+3 ≥ 3-26/10000 > 0  (|ζ|≤1/10).
```

For a fixed spatial weight the observation is the derivative of the
normalized vector `∫w0 ω`, projected onto the plane perpendicular to its
actual background direction. Its denominator is `2Ω[1+O(δ)]`, hence
nonzero. Integration by parts, including the full exterior velocity,
gives

```
|∫w0 curl(v-vapp)| ≤ ||∇w0||2 ||v-vapp||2,
||∇w0||2 ≤ Cw k^(3/2)/a,
N ≤ CN |b0|a/sqrt(k),
CN = sqrt(IF Ih1 + IG Ih/2),
Cw² = ||∇w⊥||2²||wz||2² + ||w⊥||2²||wz'||2².
```

Therefore its relative tilt error is at most a finite profile constant
times `Eδ+δ`, uniformly as the packet scales shrink.

For an actual tagged observation let `wt` be the scalar weight transported
by the unperturbed flow, and transport it also in the perturbed family.
The derivative of the weighted vorticity vector is exactly

```
δV(t) = ∫wt [curl v + ξ·∇ω] dx.                       (13)
```

This follows from `δw=-ξ·∇wt` and integration by parts with `div ξ=0`.
The second term is the moving-tag contribution and is retained. A radial
weight is invariant under `uR`. In scaled coordinates `(s,ζ)=(x⊥/a,kz)`,
the actual base flow on its swept support is `Ω(Js,0)` plus a smooth
perturbation of `C²` size at most `Cb Ω sqrtδ`: the horizontal Lundquist
correction is `O(Ωδ²)` and its vertical correction is
`O(kUδ²)=O(Ωsqrtδ)`. Choose the EPS translate so its **scaled** `C²`
correction has this bound too. For each fixed finite `a,k` this uses only
finitely many decaying EPS derivatives and is achievable by a finite `c`.

The integral equation for this scaled flow and its derivative, followed
by Gronwall on `0≤Ωt≤2π`, gives
`||wtilde_t-wtilde_0||H1 ≤ Cflow sqrtδ` and
`||wtilde_t||H1≤Cflow`, where the tildes remove the normalization
`a^-2 k`. Here `Cflow` depends only on the fixed profiles and a uniform
scaled `C²` bound, not `a,k`. The support remains in `r<a`, so `Fa=1`,
by shrinking `δ` once if needed. One explicit admissible constant is
obtained by bounding the flow derivative by `exp(2π(1+Cb))` and applying
the chain rule to the fixed weight's first two derivatives. Thus weight
transport adds at most `Ctag sqrtδ` to the relative estimate from (12).

Finally the exact Lin equation and the L² energy bound imply

```
||ξ(t)||2 ≤ exp(ΓT)[Cξ N/Ω + T N(1+Eδ)].
```

On the observation support `||∇ω||∞≤Cω Ωλ`, with `Cω` independent of
`δ` after the EPS placement. Since `||wt||2≤Cw0 sqrt(k)/a`, the ratio of
the second term of (13) to the leading vector `k ch |b0|` is at most

```
Cmat (λ/k) exp(ΓT)[Cξ+2π(1+Eδ)] = O(δ^(3/2)).
```

Combining these explicit bounded factors, actual **materially observed**
linear tilt differs from (12) by at most
`Cobs[Eδ+sqrtδ+δ^(3/2)exp(ΓT)(1+Eδ)]` in relative amplitude throughout
the period. Every constant is determined by the displayed profiles,
bounded background derivatives, and finite integration interval; taking
`δ` smaller achieves any prescribed error. No missing boundary transport
is being hidden in a fixed spatial observation. Perturbation amplitude
can independently be chosen as `|b0|=2Ω q*/(k ch)` for small physical
tilt `q*`, without changing any relative bound.

## 8. Physical and parent boundaries

Because `h=H'` has zero integral, `∫Aapp=0` and the leading complete packet
angular impulse is `∫x×vapp=2∫Aapp=0`. Its finite-core tilt is balanced by
returns. Neither a nonzero total-spin inertia nor a tube-spin constitutive
law is inferred from (12). The packet is a genuine smooth-Euler optical
orientation signal with controlled material reconstruction; it is not
the already solved parent Euler-to-Cosserat mechanism.

The geometry is a shrinking-radius, axially shorter packet in a smooth
Lundquist-like core of the same field that carries an EPS tube elsewhere.
Identifying this sector with an angle of that EPS torus, obtaining a
nonzero matched spin, and obtaining the parent's spatial constitutive
coupling remain distinct constructions. The next candidate is an EPS
interior flow-box with a verified compatible rotational jet; the present
proof does not assume that jet exists.

Route verdict: **established as stated**, finite-time Kelvin-prepared
smooth-Euler material optical tilt, with explicit parameter-dependent
error. Parent completion does not follow.
