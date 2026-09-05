# Same-EPS positive gradient action by a geometric cage lift

## Object and the additional ensemble premise

Fix an actual smooth constant-lambda EPS field and one positive physical
core-angle pair `(eta_1,eta_2)` from 0045. Its compact potentials are
`A_1,A_2`, and its exact zero-order matrices are

```
H0=[[h,g],[g,t]] > 0,  Omega0=B J,  B>0,
J=[[0,1],[-1,0]],  nu0^2=(ht-g^2)/B^2,
I0=B^2/t.
```

The field, density, cutoffs, and these physical integrals are fixed before
adding a gradient construction. When instantiating 0045 choose its cage bump
in a small part of the available annulus. Three further disjoint closed
cage balls fit within that same nonzero-vorticity neighborhood, away from
the core and original supports. All fields below remain inside the same
material EPS tube and vanish near its boundary.

Fix three geometric neighboring-frame bonds `d_j e_j`, `d_j>0`, and one
smooth nonzero bump `phi_j` in each added cage. The explicit affine-ensemble
premise is that this cage's displacement amplitude is the centered
neighboring angle difference

```
q(X+d_j e_j/2)-q(X-d_j e_j/2).
```

This is a specified, falsifiable relation between actual fluid displacements
and neighboring physical core-angle coordinates. It is not asserted to be
selected by all Euler initial data. It is the displacement-level analogue
of the parent's prescribed affine coarse-graining. No elastic coefficient
is inserted. The stiffness and all reaction terms are evaluated from the
Euler orbit action on this family.

Choose one compact circularly polarized carrier per added cage,

```
p_j=(cos(k_j z),sin(k_j z),0),  k_j/lambda>0,
Z_j=-phi_j p_j/k_j,
zeta_j=curl Z_j.
```

Its large signed carrier is chosen by the finite bounds below, not by a
desired modulus. The exact slowly varying displacement generator is

```
Xi[q,s]=curl(q A_1+s A_2
 +sum_j [q(x+d_j e_j/2)-q(x-d_j e_j/2)] Z_j).
```

It is divergence free for every smooth macrofield. For a Bloch wave with
macro wave vector `kappa`, its q-generator is

```
eta_1(kappa)+sum_j 2i sin(d_j kappa_j/2) zeta_j(kappa),
eta_i(kappa)=curl_kappa A_i, zeta_j(kappa)=curl_kappa Z_j.
```

The s-generator is unchanged. For uniform amplitudes it is exactly the
original 0045 family. Near the physical core point all added potentials
vanish. The original rotation potential is quadratic and has zero value
and first derivative at that point, so its Bloch completion preserves its
physical tilt jet as well. These facts distinguish the construction from
postulating a positive gradient energy for an abstract director.

## Complete matrix-to-action formula

For any unit macro direction `e`, expand the exact Hermitian/skew-Hermitian
forms with real-field symmetry as

```
H(epsilon e)=H0+i epsilon a(e) J+epsilon^2 H2(e)+O(epsilon^3),
Omega(epsilon e)=B J+i epsilon S(e)+epsilon^2 b2(e) J+O(epsilon^3),
```

where `S,H2` are real symmetric. The characteristic polynomial through
second order is

```
det(H-i nu Omega)
 =D-B^2 nu^2+epsilon nu L
    +epsilon^2 [P+nu^2 Q]+O(epsilon^3),
D=ht-g^2,
L=2Ba+tr(adj(H0) S),
P=tr(adj(H0) H2)-a^2,
Q=det(S)-2B b2.
```

The first-derivative drift is `v=L/(2B^2)`. It is retained. An analytic
near-identity Darboux change makes the KKS form constant, followed by
elimination of the conjugate shape and normalization of the angle inertia.
The laboratory-frame gradient coefficient is exactly

```
C(e)=[tr(adj(H0) H2(e))-a(e)^2
       +nu0^2 (det(S(e))-2B b2(e))]/t.
```

Equivalently it is `I0` times the second coefficient of `nu^2` minus
`2 I0 v^2`. This distinction matters: averaging frequency squares would
give the wrong convective stiffness. The scalar local action retains its
mixed time/space drift term; a parity/isotropic action average may cancel
that odd term, but keeps the displayed `C(e)`. Strict positivity of `C`
is sufficient for the corresponding positive gradient action after that
explicit average. No first derivative is silently deleted.

The formula is affine in `H2_11` with coefficient one. That is the useful
structural fact for the cage lift: positive added q-gradient energy cannot
be canceled by a simultaneously growing, omitted momentum-gradient term.

## Exact compact-support facts about the added jets

KKS is a local integral of generator cross products. The added cage supports
are disjoint from both base generators and from one another. Hence their
cross KKS terms vanish exactly, including Bloch returns. A single real
generator has zero self KKS at zero wave number. Its coefficient in the
q-generator starts at first order in `kappa`, so the self contribution to
KKS starts at third order. Thus **all KKS jets through order two are
unchanged** by the lift.

For `Z_e=sum_j d_j e_j zeta_j`, the added energy jets have the form

```
delta H2_11=H(Z_e,Z_e)+r11(e),
delta H2_12=r12(e), delta H2_22=0,
delta a(e)=-H(Z_e,eta_2).
```

Here `|r11|<=2 |partial_e H_kappa(eta_1,Z_e)|` and
`|r12|<=|partial_e H_kappa(Z_e,eta_2)|` at zero. Signs of these remainders
need not be guessed; their exact integrals remain in the final coefficient.
Only bounds are used in the existence proof.

## Derivative estimates on the actual EPS field

These estimates avoid an unjustified transfer of gradient signs from a
frozen constant-vorticity calculation. Let
`F_i=eta_i cross omega` and
`C_i(e)=(e cross A_i) cross omega`. On a fixed compact support of coordinate
radius at most `R`, the phase-conjugated force derivative is

```
partial_e [exp(i kappa.x) F_i(kappa)] at zero
 =i U_i(e),
U_i(e)=(e.x)F_i+C_i(e).
```

Use the exact identity

```
H_kappa(i,j)=rho [<exp(i kappa.x)F_i(kappa),
                         P exp(i kappa.x)F_j(kappa)>
 -lambda^(-1)<curl(exp(i kappa.x)F_i(kappa)),
                         exp(i kappa.x)F_j(kappa)>].
```

P is an L2 contraction. Curl is self-adjoint, so in a base/cage cross term
put it on the fixed base factor. Derivatives may be taken on the compact
phases/amplitudes rather than on a singular shifted Leray symbol. All base
matrix jets are consequently finite and differentiable to every required
order. The first cross derivative obeys

```
|partial_e H_kappa(eta_i,Z_e)|
 <=rho [M_i1 ||F_Z||_2+M_i0 ||U_Z||_2],
M_i0=||F_i||_2+||curl F_i||_2/|lambda|,
M_i1=sup_|e|=1 [||U_i(e)||_2+||curl U_i(e)||_2/|lambda|].
```

An explicit finite upper bound for `M_i1` is

```
R ||F_i||_2+W ||A_i||_2
 +(||F_i||_2+R ||curl F_i||_2
   +sum_j ||curl((e_j cross A_i) cross omega)||_2)/|lambda|,
W=||omega||_infinity on all supports.
```

For every added carrier magnitude `K_j>=1`, set

```
Fbar_j=W (||phi_j||_2+||grad phi_j||_2),
D_F=(sum_j d_j^2 Fbar_j^2)^(1/2),
D_U=R D_F+W (sum_j d_j^2 ||phi_j||_2^2)^(1/2).
```

Then `||F_Z||<=D_F`, `||U_Z||<=D_U`, uniformly in the added carrier
frequencies. Put `N_i=rho(M_i1 D_F+M_i0 D_U)` and
`Lstar=rho M_20 D_F`. If `astar` bounds the norm of the base coefficient
vector `a(e)`, the complete action-jet correction except for `H(Z_e,Z_e)`
is bounded in absolute value by

```
Rstar=2N_1+2|g|N_2/t+(2 astar Lstar+Lstar^2)/t.
```

This follows directly from the matrix-to-action formula above. Crucially,
`Rstar` is independent of the added large frequencies. It retains the full
shape, KKS, and possible chiral terms rather than assuming them small.

## A finite positive threshold, with every mutual term included

For cage j, the actual-background estimate in 0045 gives

```
H(zeta_j,zeta_j)>=rho[(1+K_j/|lambda|) A_j-C_j],
A_j=integral (phi_j omega_z)^2 > 0,
```

where `C_j` is the explicit cutoff/background norm constant from that proof.
Off-diagonal helicity terms between distinct cage supports vanish exactly
when written in force form. Their remaining Leray-energy cross terms obey
`|H(zeta_j,zeta_l)|<=rho Fbar_j Fbar_l`. Define

```
Moff=max_j sum_{l!=j} d_j d_l Fbar_j Fbar_l,
Cstar=sum_{j,l}|(C_base)_{jl}|.
```

`C_base` is the finite exact gradient tensor obtained from the original
0045 pair by the matrix formula; no positivity of it is assumed. Gershgorin
and the preceding bounds give, for every unit direction e,

```
C_lift(e)>=rho min_j d_j^2[(1+K_j/|lambda|) A_j-C_j]
             -rho Moff-Cstar-Rstar.
```

Select finite signed carriers with `k_j/lambda>0`, `K_j>=1`, and, for each j,

```
rho d_j^2[(1+K_j/|lambda|) A_j-C_j]
   >rho Moff+Cstar+Rstar.
```

Every quantity on the right is fixed by the actual EPS field, original
positive angle action, chosen bumps, and geometric bonds, independently of
the new carrier magnitudes. Since `A_j>0`, finite choices exist. The exact
finite-parameter coefficient is the full matrix-jet integral `C_lift`, not
the bound or an assigned line tension. This proves it strictly positive in
all macro directions.

## Result and boundary

The same actual stationary smooth EPS tube therefore supports both the
positive physical core-angle action of 0045 and a strictly positive
long-wave gradient action in this explicit neighboring-frame cage ensemble.
Constant-angle coefficients and the local core-tilt jet are unchanged.
The construction is an exact conditional quadratic Euler-orbit statement,
with finite norm thresholds and an explicit microscopic displacement map.

The additional neighboring-angle/cage tie is a genuine premise, not a proof
that unconstrained Euler selects that deformation family. It belongs in
the final falsifiable ensemble declaration. The calculation does not borrow
0044 coefficients, infer positivity from the constant-jet symbol, ignore
projection tails, or infer parent campaign completion. The common-angle,
material-translation, and promotion closures remain parent work.
