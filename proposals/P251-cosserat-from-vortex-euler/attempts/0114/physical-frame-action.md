# Positive fixed-frame angular packet action on the same EPS tube

## 1. Freeze the physical frame and candidate before the energy sign

Orient the source knot by its positive longitudinal fluid flow. Use its
geometric normal and binormal, with nonzero curvature as in EPS §3. To
avoid a sign ambiguity, define `τE` by the **displayed source metric**:
`ds²=...+2ε²τE(y2 dy1-y1 dy2)dα+...`. Thus its reference transverse
coordinate equation is `ydot=τE J y`, with the rotation matrix
`J=[[0,-1],[1,0]]`. In the usual Frenet convention `n'=−κt+τstandard b`,
`τE=−τstandard`. The signs below use `τE`, not an unnamed torsion sign.

Let `γ0` be the selected harmonic reference core, with physical traversal
time `T0`, and define its unwrapped geometric phase

```
Θ0(t)=∫0^t τE(s)|u(s)| ds,
r0=Θ0(T0)/T0 ≠0,
F0(t)=R[Θ0(t)-r0 t].                                 (1)
```

For the thin-tube leading template, `|u|=1`; keeping the actual traversal
factor states the physical frame convention. The periodic registration
`F0(T0)=F0(0)=I` has **zero winding** because its displayed lifted angle
has zero net increment. It is fixed from the geometry before selecting
the sign of `λ`. No rotation by `2πm t/T0` is added.

The two candidates are the distinct physical fields constructed with
small `λ>0` and small `λ<0` on this same tube geometry. EPS Theorem 6.8
permits both. Choose the candidate

```
λ Θ0(T0)<0.                                          (2)
```

The reason for this choice is the derived action sign below. It is not
achieved by changing a Floquet logarithm or reversing a coordinate axis.
The opposite candidate remains a negative fixed-frame energy branch in
this comparison. Mirroring a field changes both torsion and curl sign,
so merely mirroring the same solution would not change their product.

The 0112 Floquet continuation is anchored to the **unwrapped** `r0`, not
to the principal angle in `(0,π)`: at fixed sufficiently small tube
thickness, continue the logarithm near `r0 J` as the small geometric,
Beltrami and global-approximation corrections are introduced. The physical
periodic factor is the continuation of (1). This excludes all added
integer windings before the sign is tested. The resulting signed rate
`r` obeys `|r-r0|<|r0|/2` after the construction parameters are chosen,
so `λr<0`. Its multiplier is `exp(±irT)`; the rate is not a new lab-frame
carrier frequency or the conditional continuum optical gap.

## 2. Why this registration is a physical core-angle frame

For a material amplitude `a`, the transverse physical displacement is
`P_(t⊥)a`, where `t=u/|u|`. The possible tangent component imposed by
`k·a=0` does not create an extra transverse strain: along the trajectory,
`P_(t⊥)A t=tdot`. Therefore the tangent terms cancel when differentiating
`P_(t⊥)a`. In a normal frame `E=(n,b)`, its harmonic-limit transverse
generator is exactly

```
E^T A E-E^T Edot.                                    (3)
```

The source's scaled `C1` closeness to its reference field gives (3) as
`τE|u| J+O(ε)` along the continued core. This is the coordinate equation
used in Proposition 7.12, now applied to the physical transverse quotient.
Multiplication by the physical angle normalization and by `J` changes it
only by the corresponding explicitly differentiated scalar and frame map;
these corrections are `O(ε)` near the unit-speed harmonic reference.
The Beltrami gyro term is then bounded by the fixed-tube estimates of
0112, before choosing the free small `λ`.

The physical rotation vector is independent of the sign convention for
oriented vorticity: changing `λ` reverses both the background vorticity
axis and its normalized tilt, leaving their cross product, the rotation
vector, unchanged. The action-sign comparison in (2) therefore concerns
two genuinely different physical curl fields in the same geometric frame.

## 3. Actual packet KKS form, not a guessed oscillator symplectic area

Use the compact Kelvin-prepared generator columns from 0112, scaled so
their leading physical core angles have unit size. They are real fields
of the form

```
ηi = N^-1 ai(t,x) sin[Nφ(t,x)+π/4] + lower powers of N^-1,
Dφ=0,
```

with their exact curl completions. The initial two transverse amplitudes
are oriented by the frozen physical frame, not by the sign desired for
the energy. The full whole-space Euler forms are

```
Ω(η,ζ)=ρ∫ω·(η×ζ),
H(η,ζ)=ρ∫[vη·vζ-vη·curl(vζ)/λ],
vη=P(η×ω).                                          (4)
```

Both include the actual EPS field and its exterior Leray response.
They are the stationary fixed-Kelvin energy Hessian and KKS form, not
the unreduced material metric or an appended rigid-body inertia.

For the leading oriented transverse columns, `a1×a2` is a positive
multiple of `k`, and `k·ω=λ` in the normalization of 0112. Thus

```
Ω_N=β_N J0,        J0=[[0,1],[-1,0]]=-J,
β_N=λ b*/N²+O(λ/N³),                b*>0.             (5)
```

The constant `b*` is the actual positive envelope/normalization integral
with the factor `1/2` from `sin²`; it is not a fitted inertia. The
oscillatory remainder is bounded by integration by parts using
`∇φ/|∇φ|²` on the compact support. Curl-completion terms contribute their
explicit lower-order norms. For finite sufficiently large `N`,
`sign β_N=sign λ` and `|β_N|≥|λ|b*/(2N²)`.

The leading KKS density is conserved along the packet transport.
Indeed `Dω=Aω`, `Dai=Aai+bi`, and
`bi=−(s/|k|²)k×ai`. The three `A` terms in
`D[ω·(a1×a2)]` sum to `tr A=0`; rotating both transverse amplitudes about
`k` leaves their cross product unchanged. The envelope and phase are
materially transported. Consequently the leading integral `b*` is
constant, not assumed constant because the center's monodromy is elliptic.

## 4. Cancel the carrier in the same action before estimating its sign

Write the actual reconstruction generator as `A_Euler`, so
`H=−Ω A_Euler`. Let `E_N(t)` be the completed packet columns, with the
physical Floquet registration fixed above, after an orientation-preserving
constant near-identity normalization of the two coordinates. The anchored
Floquet normal form gives the signed constant matrix `B0=rJ`. The
high-order construction supplies

```
Edot=A_Euler E-E B0+r_N.                              (6)
```

The remainder is the actual Euler/Kelvin reconstruction remainder.
It is not the pointwise center-amplitude error alone. One can take any
fixed order in 0112 and one extra derivative: its explicit pressure/Lin
recursion and the differentiated finite-time Euler estimate give
`||r_N||2≤C N^(-m-1)` after the unit-angle scaling above. Equivalently
the relative column remainder is `O(N^-m)`. All profile, `1/λ`, tube,
and finite-time constants are fixed before choosing `N`.

For the full phase action `L=−z^TΩ zdot/2−z^THz/2`, the exact moving
pullback of main's 0115 gives

```
Q=E^TΩ Edot,
Ω_E=E^TΩE,
H_eff=E^T H E+sym Q
     =-sym(Ω_E B0)+sym(E^TΩ r_N),
dotΩ_E=Q-Q^T.                                       (7)
```

This is the actual carrier cancellation: both `A_Euler` and `Edot`
contain the transported phase's potentially large advective contribution,
which cancels in the displayed identity before any error is bounded.
Bounding the original energy and the frame term separately would not
establish this small angular-energy sign.

Since (5) is a two-dimensional skew matrix,

```
-sym(Ω_E B0)=-β_N r I=:h_N I,
h_N>0                                               (8)
```

by (2). The actual remainder obeys

```
||E^TΩ r_N|| ≤ ρ||ω||∞||E||2||r_N||2
             ≤ C |λ| N^(-m-2).                     (9)
```

Choose finite `N` so this is below `h_N/4`. Then the **full** finite
packet pullback, including its exterior response and time-dependent
connection, has positive definite angular Hessian, with lower bound
`3h_N/4`. Its symplectic form is nondegenerate. This is the sign of the
actual restricted action in the frozen physical frame, not a conclusion
from a positive Floquet eigenangle.

## 5. Exact finite-interval autonomous normalization without winding

At finite `N`, (7) is initially a positive time-dependent action. Retain
that fact rather than silently replacing its coefficients by their limit.
Its exact reduced evolution is

```
B_N=-Ω_E^-1(H_eff+dotΩ_E/2)
    =B0-Ω_E^-1 E^TΩ r_N.                             (10)
```

Thus `||B_N-B0||≤C N^-m`; also `dotΩ_E` is bounded by the same actual
remainder, since `B0^TΩ_E+Ω_E B0=0`. Define the unique matrix

```
Rdot_N=B_N R_N-R_N B0,        R_N(0)=I.               (11)
```

This is an executed exact linear transport construction, not a Floquet
logarithm choice. On the fixed optical interval Duhamel gives
`||R_N-I||≤Cfinite N^-m`. For `N` large it stays in a small contractible
neighborhood of the identity; it has no extra frame winding. Direct
differentiation proves

```
R_N^T Ω_E(t) R_N=Ω_E(0)=β_N(0)J0.                   (12)
```

The constant `B0` is symplectic for this initial form. Substituting
`z=R_N y` therefore gives, up to its exact boundary derivative, the
autonomous quadratic action

```
L=β_N(0) p qdot-h0(q²+p²)/2,
h0=-β_N(0)r>0.                                      (13)
```

Eliminating the actual conjugate coordinate gives

```
p=-qdot/r,
Langle=I_N qdot²/2-K_N q²/2,
I_N=-β_N(0)/r>0,        K_N=-β_N(0)r>0.              (14)
```

No integer frequency or sign was selected in (11); its target `B0`
was frozen from the geometric continuation before (2). The positivity
was already established in the original finite-N action in §4. Equation
(11) preserves it in an autonomous finite-interval representation close
to the prescribed physical Floquet frame. It is not an assertion of an
exact nonlinear invariant finite-dimensional Euler manifold or an exact
periodic packet at finite carrier.

For the branch with time dependence `exp(-i|r|t)`, the KKS Krein pairing
is positive in this same convention: with unit circular polarization its
value `i v*Ω_E(0)v` is `-2β_N(0) sign(r)>0`. The branch and frame are
specified; adding an integer winding and then reporting its new energy
sign would be a different representation and is not used here.

## 6. Physical angle, material spin and current remain observable rows

Every physical observation is transformed along with the action:

```
Φphysical=OΦ(t) E_N(t) R_N(t) (q,p),
Sphysical=OS(t) E_N(t) R_N(t) (q,p),
Pphysical=OP(t) E_N(t) R_N(t) (q,p).                  (15)
```

The exact material rows are those constructed in 0112: the angle uses
the transported positive-weight vorticity including `ξ·∇ω`; the spin
uses the full moving-parcel derivative, not only `∫r×v`. After (14),
each row generally contains both `q` and `qdot` and periodic geometric
coefficients. The actual Euler solutions, prepared with the same initial
Kelvin data, obey these maps to the controlled finite-time error. The
near-identity correction in (11) changes all rows and their derivatives;
none is left at its old value after the action normalization.

The actual parcel angular/current multipole remains

```
Qij=Idot_ij/2-epsilon_ijm S_m/2,
induced first momentum current=-div(Idot)/2+curl S/2,
```

with its actual evolving symmetric shape tensor and ambient pressure
reaction. Galilean mean coordinates from 0112 are a separate nondecaying
mean sector and retain their finite-parcel physical mass. No orbital
carrier current or frame/Routh term is erased by declaring (13) autonomous.

In particular the canonical momentum `I_N qdot=β_N(0)p` is **not**
identified with the whole-tube mechanical spin. Nonzero parcel spin was
constructed in 0112, but matching its tensor row to (14) is a further
physical mechanism, not a consequence of positive energy. Nor is `r²`
identified with the conditional C-CST-009 optical gap.

## Result and active consumer

Established: a physically selected curl-sign candidate on the same EPS
tube gives a positive nondegenerate angular packet action in a frozen
geometric/Floquet frame; its full finite-N pullback admits the explicit
near-identity autonomous normalization above, with positive Krein sign
and a controlled actual Euler reconstruction. Carrier, moving-frame,
physical-observation and ambient-current terms are retained.

The remaining stronger consumer is a positive physical spin/current
matching and slow-envelope continuum joining of these actual packets.
The present autonomous canonical action is not by itself that mechanical
rotor/continuum construction.

### The concrete first finite-wave-number response

The candidate physical off-diagonal response is the **full centered first
moment** `Q=Idot/2-epsilon(S)/2`, not its spin term in isolation. For a
slow Fourier modulation its actual momentum-density row is `-i Q k`;
incompressibility applies `P_k` to that row. Nonzero spin proves that at
least one such response survives even with the symmetric shape term kept:

```
Σ_(j=1)^3 |P_(ej⊥) Q ej|²
 = Σ_(i≠j) |Qij|²
 ≥ (1/2) Σ_(i<j)|Qij-Qji|² = |S|²/2.                 (16)
```

Thus at least one of the three fixed coordinate wave directions has
positive time-integrated transverse current response for the parcel
constructed in 0112. This is a physical observation coupling, not yet a
proof that the mean/angle action has a surviving off-diagonal modulus
after all pressure, mass and Kelvin reactions are eliminated.

For a unit physical angle, the small parcel has mass `O(N^-3)` and its
first moment and spin are `O(N^-5)`. The high-order estimate controls
their error relative to this small scale by `O(N^(1/2-m))`; hence the
nonzero normalized row in (16) survives for any fixed `m≥3` as `N` grows.
The **dimensional** coefficient does not remain finite automatically:
at fixed total material mass density an `O(N³)` parcel number density
gives a first-moment density of order `N^-2`, tending to zero. The
construction uses finite sufficiently large `N` and claims a nonzero
finite response there, or its explicitly normalized asymptotic row. A
nonzero continuum limit requires a separately derived physical scale
choice; it cannot be obtained by forgetting this factor or importing the
disjoint static cage's moments.

### Constructive nonconstant material spin over a core traversal

There is a stronger moment choice than the arbitrary positive covariance
used in 0112. Write `C` for covariance here, distinct from the momentum
moment `Q` in (16). Let `Mp` be the actual **particle** return and choose a real
matrix `S0` with
`Mp=S0 diag(1,Rθp) S0^-1`, where its first column is `u0` and the
periodic covector obeys `k0·u0=1`. Set

```
C0=S0 diag(c_t,c_n,c_n) S0^T,       c_t,c_n>0.          (17)
```

This is positive definite, `Mp C0 Mp^T=C0`, and `C0 k0=c_t u0`.
Consequently its transported covariance `C(t)=M(t)C0M(t)^T` is exactly
periodic at the leading small-parcel order, with `C(t)k(t)=c_t u(t)`.
The full material spin row of 0112 becomes, apart from its known positive
mass, phase and carrier factors,

```
L_t a=c_t[a×A u+u×(Aa+b)],
b=-(k·ω/|k|²)k×a.                                   (18)
```

This row is not zero on `k⊥`. To expose the fact, use physical coordinates
`u=Ue3`, `ω=we3`, `k3≠0`, and write
`a3=-(k1 a1+k2 a2)/k3`. Let `L` be the resulting `3×2` matrix in
`(a1,a2)`, omitting the positive common `c_t`. The actual Beltrami curl
constraints give the exact identity

```
L11+L22-2(k1 L31+k2 L32)/k3
   = U w k3²/|k|² ≠0.                               (19)
```

Therefore at least one real transverse vector `v` has `L0 v≠0`.
The Euler-amplitude return `FE` from 0112 has no eigenvalue one. Choose
the prepared initial polarization

```
a0=(FE-I)^-1 v.
```

This is a physical initial condition chosen from the actual matrices,
not a change to the field, Floquet branch, or action. All coefficients
and the covariance row are periodic after one core traversal. Hence

```
Sleading(T)-Sleading(0)
 = (known nonzero factor) L0(FE-I)a0
 = (known nonzero factor) L0 v ≠0.                   (20)
```

The same prepared actual Euler packet has arbitrarily small error relative
to this fixed nonzero endpoint margin by the high-order bounds of 0112.
Choose the finite observation interval to include a core traversal as
well as the optical time. Its **actual material spin is nonconstant**,
and Cauchy–Schwarz gives
`∫0^T |Sdot|² dt ≥ |S(T)-S(0)|²/T>0`.
No assertion of a generic spectral splitting or absence of all combination
resonances is needed. The invariant particle covariance makes the
observation coefficients periodic, while the nontrivial Euler-amplitude
return supplies the concrete nonconstant response. The exact material
shape/current and ambient traction terms remain those in (15)–(16).
