# Periodic coherence, complete-fluid transfer, and a transported-seed join

This is a bounded stronger-Euler construction, not a change to accepted
C-CST-008/009/010 or a derivation of their autonomous coefficients from
an unrestricted Euler solution. Density is constant. Brackets below are
normalized averages over the torus cell; multiplication by density
converts velocity momentum to physical momentum. Complex Bloch notation
denotes the real combination of opposite fibers.

## 1. Same-core periodic background and order of scales

Apply the primary inverse-localization result in `source-receipt.md` to
the strict elliptic-core and fixed-frame-sign construction of 0112/0114.
Choose the local approximation error smaller than its finitely many
strict geometric margins, in the finite differentiability norm required
by the packet order. The resulting exact periodic field has

```
curl u = λ u,   div u = 0,   p = -ρ|u|²/2,   |λ| = Λ,
```

with Λ a sufficiently large odd integer. Its knot lies in a contractible
small region of the cell. Fix this field before taking subsequent limits.
Its periodic lift has literally identical cells, not a distribution of
similar local oscillators. Integer cell translations commute with the
linearized Euler evolution. Thus a rational Bloch wave number is realized
on a finite supercell, and no random cell-frequency dispersion is present.
This assertion does not assert a point eigenvalue of a Bloch generator.

The scale order is: strict local geometry and curl sign; periodic
approximation and a sufficiently small same-shell perturbation from
Section 3; fixed finite observation time T and its smooth norm bounds;
packet order and carrier N; small nonzero rational Bloch wave number;
finally the nonlinear disturbance amplitude. No bound uniform in Λ,
all supercells, or infinite time is used.

The periodic packet is constructed directly. The local WKB recursion
uses the periodic background and compact transported initial support.
Its incompressibility completion can be written as an oscillatory
gradient plus a remainder to arbitrary chosen order. Applying periodic
Leray kills that gradient and contracts the remainder in L². This gives
the same finite-time error estimate as 0112, with constants belonging to
this fixed periodic field. In particular, this argument does not transfer
the whole-space projector by merely invoking local field approximation.
Higher Sobolev estimates follow by differentiating the same construction.

## 2. Exact Bloch mean equation, including all pressure returns

Put `v(x,t)=exp(iκ·x)wκ(x,t)`, with wκ periodic, and
`Dκ=∇+iκ`. The complete linearized Euler equation is

```
wκ_t + u·Dκ wκ + wκ·∇u + Dκ πκ = 0,
Dκ·wκ=0.                                                   (1)
```

Let `mκ=<wκ>` and `Pκ=I-κκ^T/|κ|²` for κ nonzero. Integrating (1),
using `div u=0`, periodic integration by parts and
`div wκ=-iκ·wκ`, gives the exact physical Fourier equation

```
mκ_t = -i Pκ <u⊗wκ+wκ⊗u> κ.                               (2)
```

This is the coefficient of the WHOLE fluid at Fourier wave number κ.
No isolated-tube inertia or truncated pressure Green function enters it.
At κ=0 the mean is conserved. At nonzero κ it can be sourced by the
nonzero microscopic Reynolds stress even when its initial value is zero.

Here is an explicit finite-time small-κ bound. Fix a unit direction n,
put κ=εn, and use the fiber with mean in n-perp also at ε=0. For every
nonzero integer mode l and |ε|≤1/2 the transverse projections obey
`||P_(l+εn)-P_l||≤8|ε|/|l|`. Rotate the transverse planes by their
smallest-angle orthogonal rotations. Their difference from the identity
is an order minus-one Fourier multiplier; on the mean plane use the
identity. After this identification, conjugating the first-order smooth
Euler transport changes its generator by O(ε), bounded on any fixed
Sobolev space. The same holds for the zero-order deformation and pressure
terms. A finite-time energy estimate and Duhamel therefore give

```
sup_[0,T] ||wε-w0||_(Hs) ≤ C_T |ε|                         (3)
```

for an O(ε)-close smooth preparation. Equation (2) then gives, if mε(0)=0,

```
mε(t) = -i ε P_n ∫0^t <u⊗w0(s)+w0(s)⊗u> n ds
        + O(C_T ε² t).                                    (4)
```

The direction-dependent zero-mode projector is essential; replacing it
by the identity is not the ε→0 fiber limit. The exact rational-fiber
test in `verify.py` retains this projector explicitly.

## 3. A same-background compact Kelvin seed with zero mean and positive stress

Take orthogonal directions n and a. Set

```
f=(a·u)n+(n·u)a,  J=curl(ω×P f),  Ti=curl(ω×ei)=∂iω,
ω=curl u,        ξ=curl A,       vξ=P(ξ×ω).                 (5)
```

Here P is the full mean-preserving periodic Leray operator. Two exact
adjoint identities, obtained by integration by parts, are

```
a·<u⊗vξ+vξ⊗u> n = <J,A>,    <vξ_i>=<Ti,A>.               (6)
```

We construct a field for which T1,T2,T3,J are linearly independent on
any open ball in the tube, including a ball off the observed core. The
explicit unit-curl analytic seed is

```
u*=(sin z-sin y, sin x+cos z, cos x+cos y),
n=e1, a=e3.                                               (7)
```

Its three translation derivatives are independent. Its J has a nonzero
mode at `(1,1,0)`, outside the background's unit Fourier shell. Already
the two-wave part `(-sin y,sin x,cos x+cos y)` has

```
P f=(cos y,0,-sin y),
J=(sin x sin y,cos x cos y,0).                             (8)
```

The added z wave cannot remove the displayed off-shell mode. On the
actual periodic EPS field add `b u*(Λx)`, or its matching negative-curl
version, with b arbitrarily small. This remains an exact same-curl field.
The off-shell coefficient in J is a polynomial in b with a nonzero
quadratic coefficient. Each of the three axis-shell coefficients also
has a nonzero leading coefficient. Excluding their finitely many roots
leaves arbitrarily small b preserving the strict core margins. The Ti
are then independent and J is outside their span. All these fields are
analytic, so a relation on any open ball would continue globally and
contradict the shell argument. No sampled rank test is needed.

Choose a smooth nonnegative bump χ supported on such an off-core ball,
positive on an open subset. Define normalized Gram entries and residual

```
Gij=<χ Ti·Tj>,  hi=<χ Ti·J>,
Jperp=J-Σi Ti (G^-1 h)i,
A0=χ Jperp,    ξ0=curl A0.                                (9)
```

G is positive definite by analytic independence. Equations (6) give

```
<v0>=0,
μ := a·<u⊗v0+v0⊗u>n = <χ|Jperp|²> > 0.                 (10)
```

The generator and its vorticity rearrangement are supported inside the
same tube, but v0 need not be compact: its full periodic pressure return
is included in (10). This is a Kelvin-prepared perturbation, not an
independently supplied mean velocity. The future material tags are
advected by the actual velocity, not held fixed by generator support.

The mean preparation extends exactly to small nonzero κ. The compact
potentials `Ai=χ Σj Tj(G^-1)ji` give three unit mean responses at κ=0.
Replace curls by `curlκ`, form `Pκ(curlκ Ai ×ω)`, and take the two mean
components transverse to n. Their response matrix is invertible near
zero. The corresponding seed has mean O(ε), so subtracting an O(ε)
combination of these two compact potentials makes mε(0) exactly zero.
All generators remain Bloch-divergence-free and the pressure response
remains complete. The rational κ=1/7 fixture verifies this preparation
and a nonzero exact acceleration, not an asymptotic fitted coefficient.

## 4. Nonzero physical finite-time transfer and nonlinear remainder

Let w0 be the exact linearized Euler evolution of v0. Its mean remains
zero. Define `F(t)=a·<u⊗w0+w0⊗u>n`, so F(0)=μ. On fixed finite T,
the smooth energy estimate gives an explicit finite bound
`|F'|≤L_T`, obtainable from the background derivatives and the initial
Sobolev norm. Set `τ=min(T,μ/(2L_T))`, with τ=T if L_T=0. Thus
`F(t)≥μ/2` on [0,τ]. By (4), after choosing a sufficiently small
nonzero rational ε,

```
|a·mε(τ)| ≥ |ε| μ τ/4 > 0.                               (11)
```

This is a complete-fluid Fourier transfer from zero initial coarse
momentum. Its source is exact Euler stress, including intercell pressure.
It is not yet a claim that this stress has an optical-frequency pole.

For completeness it also produces actual nonlinear Euler solutions.
On the fixed finite supercell, take real conjugate-fiber data, with
small amplitude d, prepared along the smooth coadjoint curve of the
divergence-free generator. The data are `u+d vε+O(d²)`. For integer
s≥4, the perturbation estimate has the form

```
Y_s' ≤ C_u Y_s+C_s Y_s²,
sup_[0,T] ||u_d-u-d vε||_(H^(s-1)) ≤ C_T d².              (12)
```

The first inequality, with d sufficiently small, supplies the finite
existence interval by the smooth Euler local construction; the second
follows from the difference equation and the quadratic perturbation
source. Choose d last so `C_T d²` is below the strict linear signal in
(11) multiplied by d. The actual Fourier change from its initial value
is then nonzero. An O(d²) initial preparation error is retained in this
bound, not misreported as exactly zero nonlinear initial momentum.
No global regularity or supercell-uniform existence assertion is made.

## 5. What the original optical packet does and does not supply

For a compact smooth packet with nonstationary advected phase, its
complete stress against the smooth coefficient `P f` is O(N^-M) for
each preselected M. Repeated integration by parts along
`∇φ/|∇φ|²` proves this for the WKB part; the exact periodic Euler
remainder is made smaller by selecting a sufficiently high packet order.
This bound is uniform on fixed finite T. Identically repeating cells
does not enlarge this per-volume mean.

Consequently a nonzero subparcel spin alone is not a resolved optical
residue in the COMPLETE point-Fourier momentum. This observation is not
a no-go for the actual centroid-plus-ambient momentum: the latter differs
by the intrinsic spin dipole and symmetric shape-rate multipoles, which
0117 keeps explicitly. Nor does it erase the nonzero transfer (11).

## 6. Bounded transported-seed join to the positive angular packet

The finite-time transfer can be put in the SAME two-column physical
angular family while controlling its effect on that family. This is a
small actual Euler perturbation, not a prescribed static deformation.

Use the unit-angle packet columns E_N(t) of 0114 in their fixed geometric
registration, with `B0=rJ`, and its actual reconstruction generator A:

```
E_N_dot = A E_N-E_N B0+R_N,
||E_N||_2=O(N^-1),  ||R_N||_2≤C_m N^(-m-1),
H_eff,N ≥ c N^-2 I.                                      (13)
```

Let `Y(t)=exp(tA)ξ0` be the exact smooth generator evolution of the
stress-active seed, and choose r_seed>1, for example r_seed=3. Put

```
d_N=N^(-r_seed),
E'_N(t)=E_N(t)+d_N Y(t) e1^T exp(-B0 t).                  (14)
```

The added term satisfies the SAME homogeneous equation in (13), so the
residual is still R_N. On this fixed finite time interval, boundedness
of the KKS pairing gives

```
||Ω(E'_N,E'_N)-Ω(E_N,E_N)||
       ≤ C_T(d_N/N+d_N²)=o(N^-2).                        (15)
```

Use the already derived moving-frame carrier cancellation of 0114:
`H_eff=-sym(Ω_E B0)+sym(E*Ω R_N)`. Its perturbation obeys (15) plus
the smaller residual correction. The positive angular action therefore
survives for finite sufficiently large N. This uses the full moving
action; estimating the uncorrected carrier Hessian separately would
not establish this result. The added displacement and any fixed finite
order core-angle observation change by at most `C_T d_N`; its initial
generator support is off the core. Future observations include Y(t),
not a frozen tag. This is controlled physical-observable error.

For the initial column e1, the solution of the registered amplitude
equation has `z(t)=exp(B0 t)e1`. Hence the added PHYSICAL generator is
exactly `d_N Y(t)`: the two registration factors cancel. Its mean stress
dominates the original O(N^-M) packet stress for M>r_seed and N large.
If necessary remove the packet's superalgebraically small initial mean
using the same compact mean-response directions from Section 3; this
does not affect (13)--(15) at their retained orders. Choose the packet
order so its actual evolution error is below this added signal. Then
choose ε sufficiently small that its O(ε²) Bloch error is below
`|ε|d_N μ`, and the nonlinear amplitude last. Equations (4), (11),
and (12) give a nonzero complete-fluid coarse transfer in this SAME
finite-time angular packet family, with all error scales ordered.

At nonzero ε, continuity of the two-fiber real action and its finite-time
observables preserves the strict margin (15), after ε is chosen for
the already fixed N. No spectral small-splitting approximation is used.

Crucially, (14) does NOT turn the seed into an optical laboratory
eigenmode. Its observable time trace is Y(t), while the packet has the
fixed comoving geometric phase of 0114. The bounded result is positive
angular action, actual same-core angular evolution, and nonzero actual
finite-time coarse response with controlled errors. An autonomous
optical Bloch pole with nonzero residue is a further construction.

## 7. Route verdict and next independently frozen route

Established as stated: periodic same-core coherence; exact complete-fluid
Bloch stress transfer; compact-generator zero-mean stress-active
preparation; finite-time linear and small-amplitude nonlinear transfer;
and the bounded transported-seed join (14) with its explicit physical
observation error. The analytic compact-domain Gram proof is the theorem
oracle. `verify.py` independently checks the finite Fourier pressure,
adjoints, mean corrections and nonzero rational-fiber response; it does
not identify the ABC fixture itself as a knotted tube.

The candidate claiming a resolved autonomous optical band remains
blocked at its route-specific missing construction: a Bloch resolvent
or spectral-subspace residue estimate linking the physical angular
packet to the measured complete-fluid stress at its optical frequency.
The next route is to construct that residue on the fixed periodic
background, or use the actual hybrid centroid observable with its full
spin/shape current. Neither finite-time transfer nor periodicity alone
supplies it. This is not exhaustion and not an unrestricted constitutive
closure claim. This child attempt is frozen here for bounded integration.
