# Compact-axial actual 3Ω correction: angle, spin, and action cross controls

This is a bounded construction for actual linearized Euler on uniform
rotation. Its physical observables are unchanged. Finite radial
localization leaves an explicitly bounded finite-time error, not an
exact global 3Ω eigenmode. Parent 0125 owns the slow-sideband target and
cancellation algebra; global torus/EPS transfer remains a separate task.

Use `u=Ω(-y,x,0)`, `ω=2Ωez`, Ω>0, `Jv=ez×v`, and the full R3 Leray
projector. Write `ζ=x+iy`, `e±=(1,±i,0)`. All real perturbations below
are real/imaginary parts of the displayed complex fields. The physical
core angle is `Φ=J δω_perp/(2Ω)=J∂zξ_perp` at the axis.

## 1. Exact local mode and its mechanical spin

The parent's local velocity and pressure are

```
v=e^(iNz-3iΩt)[A ζ² e_-+(2A/N²)e_++(4iA/N)ζ ez],
π=e^(iNz-3iΩt) (8iΩA/N²)ζ.                              (1)
```

Direct differentiation gives `div v=0`, the complete linearized Euler
equation including ∇π, and `curl v=-Nv`. Consequently its actual Kelvin
displacement is `ξ=i v/(2Ω)`, not an independently chosen rotor. With
unit core-angle normalization it is

```
ξ_N=e^(iNz-3iΩt)[(N/2)ζ²e_-+(1/N)e_++2iζ ez],
Φ_N(0,t)=e_+e^(-3iΩt).                                  (2)
```

This grows quadratically transversely and is only a local field.
The laboratory frequency is 3Ω. Its rotating-vector-frame inertial
frequency is 2Ω, with one further azimuthal unit; these are not a
chosen Floquet winding or a changed observation axis.

For a centered material cylinder of radius b, area D=πb², and axial
half-length l, the COMPLETE first spin variation is

```
S=ρ∫tag [r×ξ_dot+2ξ×u].                                 (3)
```

Inserting (2), integrating the transverse polynomials and then the
axial exponential gives exactly

```
S_N=j3(N) Φ_N_dot,
j3(N)=2ρD/N³ [sin(Nl)-Nl cos(Nl)
                         -(b²N²/6)sin(Nl)].              (4)
```

The small isotropic-tag limit is the parent's `I_tag/3`, not the
fundamental's `I_tag`. The finite-tag formula also proves independent
angle/spin control: the two same-frequency modes N1=π/l and N2=2π/l
have unit core angle and respectively

```
j3(N1)=2ρD l³/π²,   j3(N2)=-ρD l³/(2π²).                (5)
```

Their angle/spin response matrix is invertible. This is an exact
comparison route, without a numerical frequency splitting. The compact
profile construction below supplies a cleaner localization and the same
independent controls without choosing two axial wave numbers.

## 2. Infrared failure and its finite-return repair

The original compact odd fundamental f of 0123 has
`fhat(k)=-i B_total k+O(k³)` when `B_total=∫zf!=0`.
Thus the displacement-level inverse-square sideband correction
`fhat(k)/k²` is not L2 near k=0. A nonzero axial carrier does not repair
that integral for the complete compact profile.

An actual outside-tag return repairs it. Fix a smooth odd φ supported
in `1<|s|<2`, with `aφ=∫φ²>0`, `bφ=∫sφ>0`, and place

```
f_ret(z)=t_ret R_ret φ(z/R_ret),
t_ret=-B_total/(R_ret³ bφ).                              (6)
```

Its support is outside the observed tag, its first moment is exactly
-B_total, and its squared norm is
`B_total² aφ/(R_ret³ bφ²)`. The physical tag moment and the core slope
are unchanged. The full action includes this norm and the full pressure
return. If the fundamental's inner profile contains the 0123 tail
parameter t, its added norm is proportional to `(B0+tB1)²`: the full
matching equation remains quadratic and its strict simple root persists
for sufficiently distant finite return support, with the pressure row
retuned as in 0123. The return is not free fluid inertia.

Oddness now gives `fhat=O(k³)`. Both the inverse-square displacement
correction and its needed compact primitives are regular. This is the
failure-generated route appended to the initial single-N/envelope
candidates; no empirical selection or numerical floor is involved.

## 3. An exact compact-AXIAL 3Ω profile

More generally choose any real smooth compact EVEN H with `∫H=0`.
Define compact primitives

```
h=H',    K(z)=∫_{-∞}^z H(s)ds,    L(z)=∫_{-∞}^z K(s)ds.
```

Then h is odd with `∫zh=-∫H=0`, K is odd, and L is even. Every
displayed function is compact; a primitive with a nonzero jump at
infinity is not being called a gauge. Define

```
F_H=(h/2)ζ²e_- - K e_+ - 2Hζ ez,
A_H=-i[(H/2)ζ²e_-+L e_+].                               (7)
```

Direct calculus gives

```
curl A_H=F_H,  div F_H=0,
(J-i)F_H=2i∇(Kζ),
[u,F_H]=iΩF_H.                                          (8)
```

The full fixed-Kelvin operator satisfies
`A=-[u,·]-2ΩPJ` on solenoidal fields. Therefore

```
AF_H=-3iΩF_H,
v_H=-2iΩF_H,
π_H=-8Ω² Kζ,                                            (9)
```

with the common factor e^(-3iΩt). Equations (8)--(9) prove the exact
local Euler/Lin mode for ARBITRARY such axial profiles, including the
repaired compact fundamental. No axial Fourier cutoff, zero-k leakage,
or separate carrier-frequency fit remains in this representation.
Transverse polynomial growth is the only remaining localization issue.

For convenience use `Z_H=-iF_H`; then the physical core row is

```
Φ_H(0,t)=H(0)e_+e^(-3iΩt).                              (10)
```

## 4. Actual finite-tag spin and independent controls

The tag is the same invariant cylinder under uniform rotation. Axial
and transverse parity make its first centroid displacement zero.
Its complete spin, including the second term in (3), is

```
S_H=-iρΩ T(H)e_+e^(-3iΩt),
T(H)=3D∫tag zK dz-Ir2∫tag H dz,
Ir2=∫disk r_perp² d²r=D b²/2.                            (11)
```

The symmetric shape row is also retained. For the normalized Z_H its
only nonzero integrated components are

```
δI_(perp,z)=i Dρ[∫tag zK+(b²/2)∫tag H] e_+e^(-3iΩt),   (12)
```

and the transpose. Formula (12) follows from integrating
`ρ(r⊗Z_H+Z_H⊗r)`; it is not set to zero to obtain (11).

Choose even interior profiles supported strictly inside |z|<l. Complete
each by a negative even counterbump outside the tag, normalized to make
its total H integral zero. Inside the tag, `K(z)=∫0^z H`; hence, with
`m=∫tag H`, `m2=∫tag z²H`,

```
T(H)=D/2 [(3l²-b²)m-3m2].                               (13)
```

Take b²>3l², compatible with the original short axial tag inside a
broad transverse plateau. A nonnegative even off-core interior bump
has H(0)=0, m>0 and T(H)<0. Together with a core bump H_A having
H_A(0)=1, it yields the triangular response matrix

```
                H_A        H_S
core angle       1          0
spin factor     T(H_A)     T(H_S)!=0.                    (14)
```

Thus a finite combination supplies any prescribed angle and spin row
at the common 3Ω frequency. An individual bump need not have positive
mechanical spin-rate inertia; (14) is a controllability statement, not
a claim that every tagged subparcel is an independent positive rotor.
All coefficients are actual bump integrals and fixed tag geometry.

## 5. Fully compact radial completion and its finite-time error

Let g(r/a3) be smooth radial, one on |r|≤a3, zero beyond 2a3. Localize
with the ACTUAL vector potential:

```
Z_H^g=-i curl(g A_H).                                   (15)
```

This is compact in all three coordinates and divergence free. The
radial curl completion changes only its axial component, so

```
(Z_H^g)_perp=-ig[(h/2)ζ²e_- - K e_+].                   (16)
```

It equals Z_H throughout the tag when b<a3. In particular the initial
core angle (10) is EXACT. Its future velocity is obtained by actual
Kelvin evolution and full Leray, not by freezing (15).

The residual source can be kept wholly in the radial collar. Before
the harmless -i normalization it is

```
R_H=A curl(gA_H)+3iΩ curl(gA_H)
   =2ΩP[2iKζ∇g-(J-i)(∇g×A_H)].                         (17)
```

Let all axial profiles be supported in |z|<2R_ret, with R_ret≪a3,
and put m_abs=||H||1. Compact primitive bounds give
`||K||1≤4R_ret m_abs`, `||L||1≤16R_ret² m_abs`.
The FULL unprojected source in (17), including its quadratic transverse
growth, therefore obeys

```
||source||1 ≤ CΩ[a3²||K||1+a3³||H||1+a3||L||1]
             ≤ CΩ a3³ m_abs.                            (18)
```

In rotating-vector coordinates the exact Euler propagator is
`exp(-2Ωt P_kJ)P_k`. Its degree-zero smooth angular symbol has an
off-diagonal kernel with bounds `|∂^m K_t(x)|≤C_(m,T)|x|^(-3-m)`;
0123 proves this including the separate local delta term. The collar
is at a fixed fractional distance from the tag and rotates without
changing that gap. Duhamel gives on the tag

```
|∂_x^m δξ|≤C_(m,T) ΩT m_abs/a3^m.                       (19)
```

The transverse return is odd in z and the axial return is an m=1
angular field. Thus the zeroth-order bound is improved INSIDE the tag
to `|δξ_perp|≤C_T ΩT m_abs |z|/a3` and
`|δξ_z|≤C_T ΩT m_abs |r_perp|/a3`. The same statements hold after
fixed numbers of time derivatives with the corresponding Ω powers.
These parity cancellations are why the r² growth in (18) does not
produce an uncontrolled lever-arm error.

Consequently the actual core angle error is at most
`C_T m_abs/a3` on a fixed number of periods. Integrating the FULL
material formula (3), not just its velocity part, gives

```
|δS|≤C_T ρΩ D(b²+l²) m_abs l/a3.                        (20)
```

For the positive off-core bump in (14), its strict spin scale is
`ρΩD b²m` when b² exceeds 3l² by a fixed margin and m_abs=2m.
Its RELATIVE error is O(l/a3). Shape (12) has the analogous bound.
The physical rank in (14) therefore survives at the reference time
and remains uniformly conditioned after normalizing the angle and
spin units on the fixed finite interval. At the reference time one
can use the exact pressure-completed response matrix to match desired
rows exactly. A constant prepared combination then retains an explicit
finite-time error; it is not retuned dynamically into a forced mode.

All fields and their action norms are finite at the chosen microscopic
scales. A correction multiplied by macro wave number squared is thus
O(K_macro²) in a fixed Hilbert norm. Its SELF quadratic energy is
O(K_macro^4). This does not remove the cross terms below.

## 6. Exact symplectic cross integrals and a third control

The radial completion in (15) leaves (16) unchanged. Hence angular
integration gives the exact KKS matrix between two real profiles,
in the bases (Re Z_H, Im Z_H):

```
β_ij=2ρΩ[C0 ∫K_iK_j -(C4/4)∫h_i h_j],
C0=∫g² d²r,  C4=∫r_perp^4 g² d²r.                      (21)
```

Same-parity entries vanish and the real/imaginary cross block is β.
For a finite independent family, compact Poincare bounds give
`Gram(K)≤C R_ret^4 Gram(h)`. With a3/R_ret sufficiently large, (21)
is strictly negative definite as a profile matrix. Its 3Ω fixed
physical registration therefore has positive `-3Ωβ`. The actual
moving-action connection is included in making that statement; it is
not an equality between the raw localized energy and a guessed mass.

Against the actual 0123 fundamental
`Ξ_f(q)=(g0 f q,-(q·∇g0)F)`, F'=f compact, the exact COMPLEX cross is

```
Ω(Ξ_f(q), Z_H^g)
   =-2ρΩ C01 (∫fK)(q_x+i q_y),
C01=∫g0 g d²r.                                         (22)
```

This is generically nonzero. In the real bases its block is
`-2ρΩ C01(∫fK) I2`. It enters the total action at O(K_macro²).

One can prescribe this cross together with angle and spin. Integration
by parts gives `X(H)=∫fK=-∫F H`. Place the fixed outside-tag H
counterbump beyond the complete support of F, so it does not contribute
to X. On interior off-core bumps the spin weight is
`w(z)=D(3l²-b²-3z²)/2`, while X has weight -F(z).

For the actual smooth tapered fundamental of 0123, choose a taper
interval where F/w is not constant. Such an interval exists: constancy
on every taper segment would force f=F' to be proportional to z and
would prevent a nontrivial smooth compact taper. Pick two even bumps
concentrated at points with different F/w values. Their exact integral
ratios X/T remain different for sufficiently small fixed supports.
Their spin-zero linear combination H_X has H_X(0)=0 and X(H_X)!=0.
Together with H_A,H_S from (14), this gives

```
                H_A      H_S      H_X
angle             1        0        0
spin             T_A      T_S       0
KKS cross        X_A      X_S      X_X,
det=T_S X_X !=0.                                       (23)
```

This is a constructive rank argument from the actual taper, not an
assumed numerical rank. If a caller supplies a different f, the exact
condition is independence of the displayed three functionals; a purely
linear profile on every allowed response interval does not meet it.
The return and response supports are finite and fixed before K_macro.

The same rank persists when the reference spin row is replaced by its
exact pressure-completed value, using (20) and sufficiently distant
finite radial collar. One may thus set the leading KKS cross to zero
while matching the actual reference angle and spin. The Hamiltonian
cross still remains. With the full stationary Euler Hessian
`H_E(ξ,η)=-Ω(ξ,Aη)`, equation (17) gives

```
H_E(Ξ_f,Z_H^g)=3iΩ Ω(Ξ_f,Z_H^g)-Ω(Ξ_f,R_H).
```

After setting X=0, this is the actual pressure-residual integral, not
zero by naming the two local frequencies. Choose a3 beyond the entire
radial fundamental support. The same kernel and axial-parity estimate
gives the useful bound

```
|H_cross|≤CρΩ² ||g0||1 (∫|z f|dz) ||H||1/a3.            (24)
```

It can be made small in this explicit geometric hierarchy and is
retained in the action. The fundamental and correction are not being
declared exact global eigenmodes of different frequencies with a
contradictory nonzero conserved symplectic pairing.

## 7. Macro-order and physical-response interface

Let E0 be the full fundamental embedding and E3 the constructed
correction columns, with all actual time evolution and prescribed
physical registration retained. For a finite preparation matrix C,
`E_K=E0+K_macro² E3 C` gives EXACTLY

```
Ω_K=Ω00+K_macro²(Ω03 C+C^T Ω30)+K_macro^4 C^TΩ33 C.       (25)
```

For the complete moving action define
`Q_ab=Ea^T Ω Ebdot` and `H_ab=Ea^T H_E Eb`.
The O(K_macro²) energy cross is the corresponding cross block of
`H+sym Q`, not H alone. Equations (22)--(24) are its actual microscopic
inputs. No unchanged inertia is inferred from the self-energy order.

The parent's unwanted angle/spin sideband supplies the target rows for
(23). Its reference rows can be canceled by a genuine state preparation.
On the finite interval the residual physical error is
`O(K_macro² ε_loc)` with ε_loc given by (19)--(20), plus the parent's
own higher spatial-order error. At fixed finite microscopic scales this
is NOT O(K_macro^4). The meaning is a controllably accurate second-jet
coefficient, not exact cancellation to fourth order at fixed cutoff.
The microscopic error is made smaller than a declared coefficient
accuracy before taking K_macro→0. Achieving an exact autonomous spatial
jet still requires retaining or closing the remaining response/current
error; this attempt does not relabel it.

Established route: compact-axial exact local 3Ω modes, a smooth fully
compact Kelvin preparation, controlled finite-time physical angle/spin
rows, and independent angle/spin/KKS-cross controls with the Hamiltonian
residual kept. The exact single-N local mode remains a comparison;
its separately proposed massive-pressure localization was not needed
and is not claimed executed here. No spectral numerical design, all-time equality, global
torus transfer, or parent constitutive completion is claimed.
