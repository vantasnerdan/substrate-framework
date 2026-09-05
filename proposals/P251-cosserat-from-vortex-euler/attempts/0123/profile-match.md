# Actual rotating-profile moment match with compact pressure completion

This is a construction on exact uniform rotation, plus an explicit
continuation criterion. The global fixed-curl torus/EPS transfer is the
separate attempt 0124. Nothing here changes the accepted conditional
claims. All formulas use the actual Euler material displacement and
complete Leray operator, not an added rigid-body inertia. The action is
quadratic and the dynamical/error statements concern actual linearized
Euler solutions, not a claimed exact finite-amplitude oscillator.

## 1. Exact core dynamics, material spin, and the factor two

Write `r=(x,y)`, `Jr=(-y,x)`, and take the stationary Euler field

```
u*=Ω(Jr,0),   ω*=2Ω ez,   p*=ρΩ²|r|²/2,   Ω>0.            (1)
```

A constant axial translation can be included by following the core;
this is a coordinate change, not a claim that a Galilean-shifted
Beltrami field still has its old curl eigenvalue. Here a transverse
amplitude q is dimensionless, h has length dimension, and h'(0)=1.
The ideal axial displacement is `ξ=(q h(z),0)`. The fixed-Kelvin
reconstruction operator is

```
Aξ=-P[(u*·∇)ξ+(Dξ)^T u*].                                (2)
```

On the transversely uniform field (2) gives exactly

```
q_dot=-ΩJq,       Φ=Jq,       Φ_dot=-ΩJΦ.                  (3)
```

The angle is the physical local rotation of the vorticity direction:
`δω=2Ω∂zξ`, so `Φ=ez×δω_perp/(2Ω)` at the core. This fixes its
normalization before any action or moment comparison.

Let the material tag in the unperturbed flow be a centered circular
cylinder, transverse area D and axial interval containing the full h.
It is invariant under (1). The full first material spin variation is

```
δS=ρ∫tag [r3×ξ_dot+2ξ×u*],  r3=(x,y,z).                  (4)
```

This includes the moving-boundary contribution. Odd h and transverse
symmetry give zero first centroid displacement and

```
B=∫z h(z) dz,
δS_perp=ρ D B J q_dot=ρ D B Φ_dot,
δI=ρ D B [q⊗ez+ez⊗q].                                  (5)
```

The last expression is the symmetric second-moment/shape variation.
It is generally nonzero and is retained in a physical mean-momentum
filter. This is shear of a tag, not its rigid rotation.

For a radial envelope g, put `C=∫g² d²r` and `A_h=∫h² dz`. The complete
compact construction below has the exact KKS scalar

```
β=2ρΩ C A_h.                                             (6)
```

In the ideal physical rotation (3), the two-angle first-order action is
`β q2 q1_dot-βΩ(q1²+q2²)/2`. Eliminating the conjugate amplitude for
one physical scalar component of Φ gives

```
I_can=β/Ω=2ρ C A_h,       j_phys=ρ D B.                    (7)
```

The factors are derived from different full integrals. Equality is the
constructive moment equation `2C A_h=D B`, not an identity obtained by
renaming momentum. Omitting the 2 in ω* would produce the old false match.

## 2. Smooth core and same-frequency axial tail

Choose an odd smooth compact core h0, equal to z near zero, with support
of scale ell. It can be chosen with `h0(z)>0` for z>0 on its support and
`A0/B0≥3/4`, by making the taper occupy a sufficiently thin fraction of
the core interval. Since g=1 on the tag disk, C≥D. Consequently

```
A0=∫h0², B0=∫z h0,
c0=2C A0-D B0 ≥ C B0/2 > 0.                              (8)
```

Fix a smooth odd φ supported in `1<|s|<2`, positive for s>0. For
sufficiently large finite R/ell, set

```
h1_R(z)=R φ(z/R),   h=h0+t h1_R,
A1=R³ a1, a1=∫φ²>0,
B1=R³ b1, b1=∫sφ(s) ds>0.                               (9)
```

The supports are disjoint and the tail leaves h'(0)=1 unchanged. Every
axial profile in the transversely uniform rotating field has the SAME
equation (3); the tail does not introduce a prescribed second frequency.
The moment equation is exactly

```
F(t)=2C A1 t²-D B1 t+c0=0,
Δ=D² B1²-8C A1 c0.                                      (10)
```

The condition `R³>8C a1 c0/(D² b1²)` makes Δ positive. Its smaller root
and its derivative are

```
t*=2c0/[D B1+sqrt(Δ)]>0,
F'(t*)=-sqrt(Δ) != 0.                                    (11)
```

These are exact profile integrals and an algebraic existence choice,
not a measured-frequency fit. At fixed core and transverse geometry,

```
R³ t* → c0/(D b1),
A1 (t*)² ~ a1 c0²/(D² b1² R³),
B1 t* → c0/D.                                           (12)
```

Thus a tail with vanishing displacement norm supplies a finite measured
spin moment because its axial lever arm grows. The physical tag contains
this tail; its spin is not borrowed from fluid outside the specified tag.
The tag's shape moment in (5) changes along with its spin.

The exact-integral oracle uses convenient C2 polynomial tapers. Smooth
compact approximations can be made while leaving h0=z near the core,
preserving oddness and positivity on each positive half-support. Their
A and B integrals converge, so the strict inequalities (8), (10) and
the simple root persist. The theorem therefore uses genuinely smooth
profiles, not the fixture's finite endpoint regularity.

The same scale separation can be obtained with a fixed finite transverse
radius and shrinking ell and R, with ell/R→0 and R/a_radial→0. Taking
R large in (12) is a dimensionless R/ell argument, not a requirement
that the tail leave a preassigned finite tube.

## 3. Actual compact divergence and pressure completion

Let `g(r)=g_*(r/a)` be smooth radial, equal to one on |r|≤a and zero
for |r|≥2a. Choose a fixed 0<σ<1; the tag disk has radius b=σa,
area D, and axial half-length L=2R enclosing the profiles. Define

```
H(z)=∫_{-∞}^z h(s) ds,
Ξ_h(q)=(g h q, -(q·∇g)H).                               (13)
```

Odd h makes H compact. Because h is positive on the positive half-axis,
`H≤0` and, by integration by parts,

```
||H||_1=-∫H=∫z h=B.                                     (14)
```

Direct differentiation gives `div Ξ_h=0`. Its axial component vanishes
throughout the tag, where it has exactly the angle, spin ansatz and
shape rows in (5). The cross product of its two transverse columns has
z component `g²h²`, so (6) is EXACT, including the entire return region.

For radial g, remove the gradient of `u*·Ξ` in (2). The exact residual
of the rotating profile is

```
AΞ_h(q)-Ξ_h(-ΩJq)
       =-2Ω P[ez (Jq·∇g)H].                             (15)
```

This is the actual pressure-completed residual, not a divergence error.
Leray contraction gives the useful global field estimate

```
||residual||_2 ≤ 2Ω |q| ||∇g||_2 ||H||_2,
||H||_2 ≤ 4R ||h||_2.                                   (16)
```

Relative to the transverse profile norm its size is O(ΩR/a). The
weighted spin estimate below is stronger than applying Cauchy--Schwarz
to (16): the radial lever arm is retained explicitly.

## 4. Complete pressure changes a linear row, so the exact initial root survives

Inside the separated tag, the unprojected source in (15) vanishes.
The residual is therefore

```
residual=2Ω∇ψ,
Δψ=(Jq·∇g)h,                                           (17)
```

with the whole-space Newton potential. No wall or discarded exterior
pressure is used. At the initial time the actual material velocity is
AΞ, not merely Ξ(-ΩJq). Radial and axial reflection symmetries imply

```
S_perp(0)=ρΩ D B_tilde(h) q,
B_tilde(h)=B(h)+[2/D]∫tag [r3×∇ψ(q=ex)]_x.              (18)
```

In particular B_tilde is a COMPUTED LINEAR functional of h. The
off-diagonal transverse row and the scalar-angle connection vanish at
this reference time by the same symmetries. One does not erase the
pressure correction by replacing B_tilde with B.

For completeness its smallness follows directly from the Newton kernel.
The collar is at transverse distance at least (1-σ)a from every tag
point. Since h is odd, the zeroth axial moment of the source vanishes.
Differentiating the kernel and pairing opposite axial points gives

```
|∂_r ψ|≤C B |z|/a³,
|∂_z ψ|≤C B |r|/a³
```

on the tag. Here `||∇g||_(L1(R²))=O(a)` and the first absolute axial
moment is B. The m=1 angular symmetry gives zero on the axis; oddness
gives zero at z=0, supplying the displayed factors. Integrating the
ACTUAL spin moment yields

```
|B_tilde-B| ≤ C B[(L/a)+(L/a)³].                         (19)
```

Both B0 and B1 obey this estimate with their positive h_i. In (10)
replace B0,B1 by B_tilde0,B_tilde1 and keep A0,A1,C unchanged:

```
F_complete(t)=2C A1 t²-D B_tilde1 t
                         +2C A0-D B_tilde0.              (20)
```

Choose R/ell with a strict discriminant margin, then a/R sufficiently
large. Equations (8), (19) preserve the positive constant term, positive
linear moment and positive discriminant. Equation (11), with these
complete coefficients, supplies a positive SIMPLE exact root. This is
the actual compact, pressure-completed initial moment match.

To see that it matches the physical scalar canonical momentum, not
merely the leading notation, take `theta=Φ_x=-q2`. Its initial angle
row is c=(0,-1), while (18) gives the physical spin row
s=(ρΩD B_tilde,0). The exact KKS form is βJ0. Thus

```
det([c;s])=β  <=>  F_complete(t)=0.                      (21)
```

The 0119 physical scalar chart gives `M=β/W` and
`I_spin=det([c;s])/W`, including the actual rate W and full moving
action. Equation (21) proves `I_spin(0)=M(0)` exactly. The pressure
correction to the initial core-angle rate is bounded by CΩB/a³, so
W has the positive sign of Ω. An equality obtained by substituting the
uncorrected Ω into both inertias would be weaker than (21).

## 5. Uniform finite-time physical error, including weighted spin and shape

The exact generator evolution is not confined to the ansatz (13).
Nevertheless its physical error is controlled for any fixed finite ΩT.
The following proof keeps the full nonlocal evolution.

Use physical-vector rotating coordinates

```
ξ(x,t)=R_(Ωt) ζ(R_(-Ωt)x,t).
```

The divergence-free Lie derivative identity converts (2) EXACTLY to

```
ζ_t=-2Ω P J ζ.                                          (22)
```

On a nonzero Fourier mode, the propagator is
`exp(-2Ωt P_k J)P_k`; equivalently restrict to k-perp and use P_k J P_k.
This is a smooth degree-zero symbol on the sphere. Its time and angular
derivatives are bounded for fixed ΩT. Its nonlocal kernel away from the
diagonal therefore satisfies

```
|∂_x^m K_t(x)| ≤ C_(m,T) |x|^(-3-m),   x!=0.             (23)
```

One proof decomposes frequency into annuli, integrates by parts on each
annulus, and sums using the symbol's degree-zero scaling; the local
delta term is separate and has no contribution across the collar gap.
This is an analytic singular-integral estimate, not a sampled propagator.

Compare the exact evolution from (13) with the rotating ansatz
`q(t)=exp(-ΩJt)q(0)`. In rotating coordinates its residual source is the
compact collar field in (15). Its L1 norm is at most
`C Ω a B |q(0)|` by (14). Rotations preserve the radial collar/tag
separation. Duhamel and (23) give, on the entire tag,

```
|∂_x^m δξ(t)| ≤ C_(m,T) ΩT B |q(0)| / a^(2+m).          (24)
```

The same estimate holds after fixed numbers of time derivatives with
the corresponding powers of Ω; differentiating the source and the
bounded angular multiplier supplies those powers. In particular the
core-angle and its first two time derivatives have dimensionless
relative errors bounded by C_T B/a³.

Apply (24) to the FULL material moment (4), including 2δξ×u*.
Since |r3|≤C(a+L), |u*|≤Ωa and the tag volume is 2LD,

```
|δS(t)| ≤ C_T ρΩ D B |q(0)| (L/a),
|δI(t)| ≤ C_T ρ D B |q(0)| (L/a),                       (25)
```

for L/a bounded above by a fixed small constant. The reference spin
scale is ρΩDB, not an absolute tolerance independent of the small tag.
Odd/azimuthal symmetry preserves the zero centroid variation. All shape
moments and their derivatives remain explicit physical observations.

For the completed-root profile, B is comparable to `(C/D)A0` and hence
to ell³ at fixed transverse shape. Thus both error scales in (24)--(25)
can be made small with ell≪R≪a. The exact evolved two-column KKS form
is conserved and nondegenerate. Apply the actual angle/rate chart of
0119 to these evolved columns; its Wronskian stays close to Ω, and its
complete scalar kinetic coefficient M(t) is positive. The physical
spin has its exact decomposition

```
S_axis(t)=I_spin(t) theta_dot+chi(t) theta,
I_spin(0)=M(0),
|I_spin(t)/M(t)-1|+|chi(t)/(Ω M(t))|
                       ≤ C_T [L/a+B/a³].                (26)
```

One may likewise bound the scalar stiffness using the next controlled
time derivative. The geometric connection chi and every residual in
(26) are retained. The exact initial root does not make them identically
zero over time.

The factor 2Ω in the COROTATING equation (22) is not the laboratory
core-angle frequency. On the ideal transverse plateau the physical
vector rotation back to the lab gives
`R_(Ωt) exp(-2ΩJt)=exp(-ΩJt)`, precisely (3). No Floquet winding,
new optical frequency, or fitted phase has been introduced.

## 6. Transfer interface, verdict, and next construction

The parent has assigned the actual fixed-lambda torus/EPS transfer to
0124. The exact reference-time condition (21) is especially transferable:
KKS, the core-angle row, and the full initial material-spin row are
continuous functionals of the smooth background on the compact support,
with the SAME full Leray operator. If their first parameter derivatives
are sufficiently close, the nonzero derivative in (11)/(20) gives a
nearby exact reference root by the implicit function theorem. A general
perturbed geometry has a chosen scalar observation/axis; additional
vector rows are retained, not inferred from the circular symmetry.

The entire finite-time estimate requires more than local field closeness:
the complete nonlocal Euler response and weighted observations must be
controlled, as they were in (22)--(25). This is the concrete transfer
task in 0124, not an extra test imposed on the proved uniform-rotation
construction. In particular an axial speed defect of size Ωλa² acts
on the short core profile through ∂z, so its relative advection error
contains `λa²/ell`, not merely λa. A compatible hierarchy is

```
δ=λa≪1, ell=a δ^(1/2), R=a δ^(1/4),
ell≪R≪a,   λa²/ell=δ^(1/2)≪1,
1/(λ R_major)≪1.
```

These are explicit candidate transfer scales, not a claim that global
EPS propagation has already been estimated here.

`route_verdict: established` for the smooth profile root, exact compact
reference-time physical/canonical match, and uniform-rotation finite-time
relative-error theorem (26). `evidence_scope:
EXACT_PROFILE_AND_REFERENCE_MOMENT_MATCH_WITH_CONTROLLED_FINITE_TIME_
UNIFORM_ROTATION_EULER_ERROR`.

The separate all-time equality route is not asserted: finite radial
localization introduces inertial-wave dispersion, and one scalar
initial root is not a functional identity in time. The next achievement
is the 0124 torus/EPS transfer with its full nonlocal error. An optional
higher-primitive return profile is available if that transfer requires
extra vanishing moments; it is not needed for the separated-collar
estimate proved here and is not silently included in this profile.
