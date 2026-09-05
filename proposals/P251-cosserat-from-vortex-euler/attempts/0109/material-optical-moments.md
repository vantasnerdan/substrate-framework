# Genuine Euler optical motion with nonzero material-parcel spin

## 1. What is constructed

For every prescribed positive relative error, a finite transported parcel
in one smooth stationary EPS-compatible Beltrami Euler field has an actual
Kelvin-prepared linearized Euler motion whose physical core orientation,
centroid, momentum, and intrinsic angular momentum follow the explicit
nonzero optical family below throughout one period. The inertia is derived
from the full material moment, including moving-boundary terms. The full
Euler and material-transport residuals are estimated; no projected ODE is
substituted for Euler evolution.

This is a tagged near-axis parcel in a straight Lundquist-like flow-box,
not the distant invariant EPS torus and not a claim that an arbitrary tag
defines a persistent vortex structure. Its translation and angle are
linked in a two-dimensional optical family, not independent continuum
fields. Ambient compensation and pressure exchange remain present.

## 2. The compact net-spin obstruction and the executed continuation

There is a real obstruction to simply changing the envelope while keeping
both the Kelvin generator and velocity compact about uniform vorticity.
Suppose `v,ξ∈C∞c(R³)`, both solenoidal, and
`curl v=2Ω ∂z ξ`. Put `G(x,y)=∫v(x,y,z)dz`. Integrating this curl equation
gives `∂x Gz=∂y Gz=0` and `curl2 Gh=0`. Compactness gives `Gz=0`, while
integrated incompressibility gives `div2 Gh=0`. The compact harmonic vector
field `Gh` is zero. Next `Lh=∫z vh dz` is compact and
`div2 Lh=∫vz dz=0`, so `∫Lh dxdy=0`. Consequently

```
∫v dx=0,                    ∫x×v dx=0.                 (1)
```

This closes candidate A only in that compact, constant-vorticity,
fixed-Kelvin class. It says nothing comparable for every variable-vorticity
field or noncompact Leray response. It is not an Euler no-go for physical
subparcel spin. Candidate B below explicitly transports the parcel and
computes its compensating ambient exchange. No finite-domain slip wall or
extra circulation is inserted to evade (1).

## 3. A new inner envelope and the exact uniform-rotation comparison

Keep the profile notation, `Jb=ez×b`, and Bessel estimates of 0101, but
choose a **new** odd compact primitive

```
H(ζ)=χ(ζ)(ζ-ζ³/6),
χ even and smooth, χ=1 on |ζ|≤1/2, χ=0 on |ζ|≥1;
h=H',       h'=-ζ, h''=-1 on |ζ|≤1/2.
```

Here `F(x⊥/a)` is radial, one on radius one, zero outside radius two;
`b(t)=R(-Ωt)b0` is horizontal. All new profile norms are evaluated with
this `H`, not borrowed numerically from the different 0101 envelope. Set

```
vapp = curl[k^-1 Jb Fa h(kz)],
W = -b Fa h(kz)/k + ez(b·∇Fa)H(kz)/k²,
ξapp = curl W/(2Ω).                                   (2)
```

They are smooth compact solenoidal fields. `∂zW=vapp`, hence
`P(ξapp×2Ωez)=vapp` exactly. On the inner slab where `Fa=1`,

```
vapp=b k z,        ξapp=q(t)z,
q=k Jb/(2Ω),       Φ=Jq=-k b/(2Ω),
qdot=-ΩJq,         Φddot+Ω²Φ=0.                        (3)
```

`q` is the actual axial-vorticity **tilt**, while `Φ` is its infinitesimal
physical rotation vector: `Φ×ez=q`. This distinction is essential for
the sign and phase of the physical spin. The inner displacement shears
the fluid; it is not declared to be a rigid rotation of the whole parcel.

Select constants `c=1/4`, `ell=1/8` and the initial material ball

```
D0 = {x²+y²+(z-zc)² < R²},
zc=c/k,     R=ell/k,
M=4πρR³/3,  σz²=R²/5,    j=Mσz²>0.                 (4)
```

For `ka≥1`, this lies strictly inside the polynomial core of (2). It is
invariant under the comparison flow `uR=ΩJx`, with centroid `zc ez`.
The true selected material set is instead transported by the actual Euler
flow; its control is proved below. There is no imposed physical boundary.

## 4. Full material moments and the positive physical inertia

For any material parcel with baseline centroid `X`, mean velocity `V`,
and displacement `ξ`, differentiation of the exact definitions gives

```
U=δX=<ξ>,            δP=M< Dt ξ >=M Udot,
δS=∫D ρ[(ξ-U)×(u-V)+(x-X)×(Dt ξ-Udot)] dx.           (5)
```

For the invariant comparison parcel, integration by parts also gives
`δS=ρ∫D[(x-X)×ξt+2ξ×uR]`, as in 0084. Applying it to (3) yields

```
U=zc q=-zc JΦ,
δP=M zc qdot=-M zc JΦdot,
δS=Mσz² Jqdot=j Φdot.                                (6)
```

All three responses are nonzero for `b0≠0`. Specifically,
`|U|=zc k|b0|/(2Ω)`, `|δP|=M zc k|b0|/2`, and
`|δS|=Mσz² k|b0|/2`; none relies on a symmetry-breaking numerical error.
The `2ξ×uR` integral vanishes by the transverse first moments. The other
term uses `∫(z-zc)z=Mσz²/ρ`.

The fixed-domain Eulerian velocity moment alone is **twice** the answer:
`ρ∫D(x-X)×vapp=2j Φdot`. The missing moving-boundary contribution is
`-j Φdot`, and their sum is (6). Thus `j` is not a borrowed rigid-body
inertia: the ball's transverse rigid inertia is `2MR²/5`, whereas the
actual prepared fluid mode has `j=MR²/5`. Nor is an old coadjoint Hessian
used to infer its physical spin.

## 5. Ambient pressure exchange is part of this same motion

For the actual material parcel, the exact Euler momentum balances are

```
Pdot = -∮∂D p n dA,
Sdot = -∮∂D p (x-X)×n dA.                            (7)
```

In the comparison flow `pR=ρΩ²(x²+y²)/2`. The pressure correction of the
explicit comparison packet constructed below vanishes on the parcel.
The first variation of its force is `-MΩ²U`. The first variation of its
torque is

```
-ρΩ²∫D [(ξ-U)×x⊥ +(x-X)×ξ⊥] dx
 = -j Ω²Φ = d(j Φdot)/dt.                             (8)
```

The first integral term is zero by transverse symmetry; the second uses
the same positive axial variance as (6). Thus the restoring force and
torque are actual background pressure traction on the displaced material
boundary. They are neither an isolated free rotor nor an inserted spring.
For the true Euler solution (7) holds exactly, with its complete pressure
reaction. The finite-time moment approximation below retains that solution,
rather than replacing its pressure by `pR` everywhere. Opposite tractions
act on adjoining fluid; the compact comparison packet still obeys (1).

## 6. Stronger pressure and Lin residuals: the analytic repair

The unprojected rotation residual from (2) is
`-2Ω ez (Jb·∇Fa)h/k`. Define the explicit pressure per unit density

```
πapp = 2Ω (Jb·∇Fa)H(kz)/k².
```

An exact cancellation gives

```
vapp,t+(uR·∇)vapp+(vapp·∇)uR+∇πapp
 = 2Ω H(kz) ∇⊥(Jb·∇Fa)/k².                          (9)
```

Leray projection still acts on the right side, including its exterior
reaction; its L² contraction gives a genuinely stronger estimate. Write
`N=||vapp||2` and use the fixed profile integrals of 0101, now with this
new `H`. For `ka≥1`,

```
N²=|b0|²[a² IF Ih1/k+IG Ih/(2k³)],
||Rvelocity,rotation||2/N ≤ Ω Cp/(ka)²,
Cp=2 sqrt(IF2 IH/(IF Ih1)).                            (10)
```

Direct differentiation of the **full** compact material field also gives

```
ξapp,t+[uR,ξapp]-vapp
 = H(kz) J(D²Fa)Jb/k².                               (11)
```

Hence its L² norm is at most `N CL/(ka)²`, `CL=Cp/2`. The residual
vanishes in the inner parcel, but its full support is included in the
global error bound. This is not a zero-complement assumption.

Use finite profile constants from separated derivatives to bound
`||ξapp||2≤Cξ N/Ω`, `||Dξapp||2≤Cξ1 kN/Ω`, and
`||Dvapp||2≤CD kN`. For example the first two are obtained from (2) by
triangle inequalities using derivatives of `F` through order three and
`H` through order three; all powers of `1/(ka)` are bounded by one.

## 7. The actual same-field Euler solution and finite-time bounds

As in 0101 use `u=uE+uL(·-cE)` with the same curl eigenvalue `λ`, choosing
the Lundquist amplitude small enough for EPS persistence and the translate
far enough to make the EPS correction small on the finite swept flow-box.
The field is exactly stationary, smooth and bounded, though its background
energy on R³ is not finite; every perturbation norm used here is finite.
Its knotted torus is elsewhere. In the axial Galilean frame write
`u=uR+d` on the packet support. Put `δ=λa≤1/2`. The achieved bounds are

```
||d||∞≤3Uδ², ||Dd||∞≤4Uλδ, ||ω-2Ωez||∞≤4Ωδ,
Ω=Uλ/2,       T=2π/Ω.
```

Prepare the exact actual Kelvin state with
`v(0)=P[ξapp(0)×ω]`, `ξ(0)=ξapp(0)`. Let `v` be its unrestricted
linearized Euler solution and `ξt+[u,ξ]=v` its actual material
reconstruction. On a fixed finite interval, these are the derivatives of
actual smooth Euler flows with sufficiently small perturbation amplitude.
Choose a finite global derivative bound `Γ≥||Du||∞`, independent of the
packet scales and translation. The full Euler energy estimate yields

```
sup ||v-vapp||2/N ≤ Ev,
Ev=e^(ΓT){4Cξδ+2π[Cp/(ka)²+6CD δka+8δ]}.             (12)
```

Using (11), the actual Lin residual of `ξapp` has norm at most `N RL`,
where

```
RL=CL/(ka)²+6Cξ1 δka+8Cξδ.
```

Subtract the two Lin equations and use the same transport estimate:

```
sup Ω||ξ-ξapp||2/N ≤ Eξ=2πe^(ΓT)(Ev+RL).              (13)
```

On the parcel and a surrounding core ball the material-velocity difference
also has the bound

```
||Dt ξ-Dt,R ξapp||2/N ≤ Em,
Em=Ev+(Γ/Ω)Eξ+8Cξδ+CL/(ka)².                         (14)
```

This follows from `Dt ξ=v+(Du)ξ` and (11); it retains both the velocity
and reconstruction errors. Set the finite scales

```
a=δ/λ,       ka=δ^-1/3,       k=λδ^-4/3.             (15)
```

Then `Ev,Eξ,Em≤C* δ^(2/3)`, with explicitly bounded fixed-profile and
finite-time Gronwall constants. This improves the older first-order
pressure bound precisely where the physical moment needs it. Neither
the large Gronwall factor nor the inverse moment sensitivity is dropped.

## 8. Actual transported parcel, angle and moment sensitivity

The material ball has radius `ell/k`, not the packet's larger transverse
support radius `a`. On this smaller ball introduce isotropic coordinates
`y=k(x-zc ez)`. The actual base flow differs from `ΩJy` by a `C²`
perturbation of size `O(Ω λ/k)` on a fixed containing ball: the Lundquist
axial correction is `kU O((λ/k)²)=O(Ω λ/k)`, while the horizontal
correction is smaller. Choose the EPS translate so its scaled `C²` norm
has the same finite bound. All required physical derivatives decay; this
is achieved by placement, not an unknown trajectory hypothesis.

The flow integral equation and its first derivative imply that the scaled
transported parcel differs from the rotating reference ball in `C¹` by
`Cflow λ/k` through `T`. Its volume is exactly constant. Its diameter is
at most `C/k`; its relative symmetric-difference volume is at most
`Cflow λ/k`. Its centroid and baseline mean velocity are retained when
using (5). These bounds follow with constants depending on the fixed
scaled ball, `ΩT=2π`, and a uniform scaled `C²` bound, not on `a,k`.

The elementary Cauchy–Schwarz sensitivity is important:

```
sqrt(|D0|) ~ k^-3/2,
N ≤ CN |b0| a/sqrt(k),
N/(|b0|sqrt(|D0|)) ≤ C (ka).                          (16)
```

Using (12)–(14) in (5), and bounding the *known smooth comparison fields*
on the small symmetric difference of the domains, gives the actual bounds

```
|U-Uapp|/|Uapp| ≤ Cmom[(ka)(Eξ+Ev)+λ/k],
|δP-δPapp|/|δPapp| ≤ Cmom[(ka)(Em+Eξ)+λ/k],
|δS-δSapp|/|δSapp| ≤ Cmom[(ka)(Em+Eξ)+λ/k].            (17)
```

The constants contain the fixed nonzero `c,ell`, so no vanishing leading
moment is hidden in a relative estimate. The spin test has lever arm
`O(1/k)` and background centered velocity `O(Ω/k)`, giving exactly the
same `ka` sensitivity, not an unrecorded stronger scale loss. Actual
velocity and displacement errors are integrated over the actual parcel;
only the explicit comparison fields are moved across the small symmetric
difference. This avoids requiring a trace estimate on the unknown error.

For physical orientation, take a smooth positive radially symmetric
sampling density carried with the parcel and initially supported strictly
inside this ball. Normalize its weighted actual vorticity vector and let
`Φphysical=J` times its linear tilt. The exact variation includes
`∫wt[curl v+ξ·∇ω]`, not merely a fixed spatial curl. Since the comparison
tilt is constant on the ball, it is exactly (3) for any such weight.
The normalized weight has gradient norm `O(k^(5/2))`, so integration by
parts and (16) give

```
|Φphysical-Φapp|/|Φapp|
 ≤ Cobs[(ka)(Ev+Eξ)+δ+λ/k].                           (18)
```

The moving-weight term uses `||∇ω||∞≤CΩλ` and the Lin bound, and is
`O(λa)=O(δ)` in relative amplitude. For its time derivative use the exact
material vorticity equation `Dt ω=(ω·∇)u`. Its first variation is

```
δmaterial(Dt ω)
 = (curl v+ξ·∇ω)·∇u + (ω·∇)v + ω·[(ξ·∇)Du].          (19)
```

Integrate derivatives of `v` against the transported smooth weight.
`||Du||≤CΩ`, `||D²u||≤CΩλ` on this core, the weight derivative scale is
`k`, and the actual normalized-vector derivative has a denominator bounded
below by `Ω`. This gives the same relative error bound as (18) for
`Φdotphysical`, measured against `Ω|Φapp|`; no time derivative of the
rough global error or extra carrier power is required. Combining (17)
and (18)–(19) proves

```
|δS-j Φdotphysical|/(j Ω|Φapp|)
 ≤ Cjoint[(ka)(Ev+Eξ+Em)+δ+λ/k] = O(δ^(1/3)).          (20)
```

Likewise `U+zc JΦphysical` and
`δP+Mzc JΦdotphysical` have vanishing relative errors. Select `δ>0`
small enough using the displayed finite constants that their common bound
is below any target `ε<1/2`. This preserves nonzero actual spin and
momentum throughout the period and gives a positive physical inertia with
controlled error. All choices are finite; amplitude can independently be
scaled to keep the physical tilt small. No claim of efficient computation
or of infinite-time invariant-mode closure follows.

## 9. Verdict and next consumer

Established: an exact-Kelvin-prepared actual Euler optical family with
nonzero **actual material-parcel** momentum and spin, full moving-boundary
and ambient reaction, and arbitrarily small finite-time error in its
positive spin/rate relation. The pressure parametrix repairs the moment
accuracy problem exposed by the first coarse estimate.

The entire packet still has compensating returns, and (1) explains why
they cannot simply be erased within the compact comparison class. This
does not identify the parcel with the persistent EPS knot, create an
independent six-dimensional translation/rotation sector, or prove a
homogenized continuum. Those are separate parent constructions. It does
supply genuine dynamical and physical-moment content that an arbitrary
positive projected oscillator would not provide.
