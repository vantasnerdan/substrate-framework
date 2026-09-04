# Full three-dimensional finite-core twist energy

## Exact object and infrared normalization

Let `f=f_1+f_2` be the smooth six-core transverse vorticity from attempt
0036; `f_i` consists of one alternating triangle. Each component core has
circulation `Gamma>0`, centroid `X_a` at radius `d`, and scale `epsilon`.
Each entire triangle has zero first moment, and the supports of `f_1,f_2`
are disjoint. Coordinates are relative to the molecular axis.

Work on the plane times a periodic axial coordinate `z` of period `L`.
The exact three-dimensional Euler kinetic energy is
`rho/2 <omega,(-Delta_3)^(-1)omega>`. Its axial Fourier kernels are

```
G_k=(-Delta_perp+k^2)^(-1),
G_k(x-y)=K_0(|k||x-y|)/(2 pi)  (k nonzero),
G_0(x-y)=-log(|x-y|/ell)/(2 pi).
```

The nonzero-circulation straight state has a logarithmic infrared energy.
Define energy differences at fixed total circulation using the last kernel
for the zero axial mode. Changing `ell` adds a constant depending only on
that fixed circulation and cancels from all variations used below. This
defines the renormalized whole-space Euler energy, not a screened fluid,
wavenumber-dependent core, or frozen external work source. The K0 kernel
also follows directly by Fourier transforming the three-dimensional Newton
kernel in `z`; the cosine integral and its logarithmic expansion are recorded
in [NIST DLMF 10.32.6](https://dlmf.nist.gov/10.32.E6) and
[10.31.2](https://dlmf.nist.gov/10.31.E2).

The chosen twist has a finite local-gradient coefficient without an outer
cutoff because each triangle has zero first moment. The molecular scale `d`
provides a convenient reference inside logarithms, not a dropped region of
the Biot--Savart integral.

## The complete isovortical second variation

On triangle `i`, push forward by the volume-preserving map
`x_perp=R_{q_i(z)} a_perp`, `z=a_z`. For small angles the disjoint motions
extend to an ambient volume-preserving diffeomorphism, as in 0036. The exact
pushed vorticity is

```
omega_i=(q_i' f_i(R_{-q_i}x) Jx, f_i(R_{-q_i}x)),
J(x1,x2)=(-x2,x1).
```

Define `A_i=f_i Jx`, `B_i=L_z f_i`, where `L_z=Jx.grad_perp`. In
particular `B_i=div_perp A_i`. Under `q_i -> t q_i`,

```
delta omega_i=(q_i' A_i,-q_i B_i),
delta^2 omega_i=(-2 q_i q_i' B_i Jx, q_i^2 L_z^2 f_i).
```

Both variations are divergence free. The full second derivative is

```
H''=rho <delta omega,G delta omega>
     +rho <omega_0,G delta^2 omega>.
```

The second term is retained. Its horizontal part pairs to zero with the
vertical base field; its vertical part supplies the angular-stiffness
correction below. Dropping it would incorrectly retain energy for a uniform
common rotation. A different incomplete calculation, retaining only the
horizontal first-variation norm and omitting the longitudinal resolvent
change, would double the singular self contribution to the twist-gradient
coefficient. The two omissions have distinct mechanisms.

Use axial Fourier coefficients normalized so that the quadratic energy per
length is `1/2 sum_k qhat(k)^* D(k) qhat(k)`. Then exactly

```
D_ij(k)=K_ij+k^2 C_ij(k),
K_ij=rho [<B_i,G_0 B_j>+delta_ij <f,G_0 L_z^2 f_i>],
C_ij(k)=rho [<A_i,G_k A_j>-<B_i,G_k G_0 B_j>].
```

Here `K` is exactly the finite-core straight angular Hessian, already
constructed from Euler interaction integrals in 0036. The formula follows
from the resolvent identity `G_k-G_0=-k^2 G_k G_0` on the zero-mean source
space. In particular the angular-impulse term used to pass to the rotating
action has no axial-gradient contribution.

## Exact positivity: the transverse projector

Let `xi` be transverse Fourier momentum and
`P_T(xi)=I-xi xi^T/|xi|^2`. Since `B_i=div A_i`, the full gradient matrix
has the exact Gram representation

```
C_ij(k)=rho integral
  [P_T Ahat_i]^* . [P_T Ahat_j]/(|xi|^2+k^2) d^2xi/(2 pi)^2.
```

It is positive semidefinite. It is strictly positive definite for the two
triangle twists: if a real linear combination had zero norm, the associated
`sum_i c_i A_i` would have zero curl. The triangle supports are disjoint, so
each nonzero `c_i` would require `curl A_i=div(x f_i)=0`. But integration
by parts gives

```
integral |x|^2 curl A_i = -2 integral |x|^2 f_i < 0,
```

a contradiction. The strict statement is not merely a leading-log sign
guess. It holds for each nonzero `k`, and at `k=0` where the norm is finite.

For finiteness at zero, `integral A_i=J integral x f_i=0`, so
`Ahat_i=O(|xi|)` at small momentum. Compact smooth cores control large
momentum. Dominated convergence therefore gives the finite, exact local
coefficient `C_ij=C_ij(0)`. Moreover
`C(k)=C(0)+O(k^2 |log(kd)|)` at fixed core profile: split the projector
integral at `|xi|=|k|` and use `Ahat=O(|xi|)` near zero. Thus the first
neglected energy contribution is of order `k^4 |log(kd)|`, within the
original slow-varying derivative expansion. No equality of a local action
with every finite-wavenumber Euler motion is assumed.

## A real-space finite-core coefficient, including every mutual term

Using the biharmonic kernel
`Q(x-y)=|x-y|^2 log(|x-y|/ell)/(8 pi)` for `G_0^2` and integrating each
angular derivative by parts gives

```
C_ij = rho/(4 pi) integral integral f_i(x) f_j(y)
 [ -(x.y) log(|x-y|/ell) + (x cross y)^2/|x-y|^2 ] dx dy.
```

An intermediate additional `(x.y)/2` integrates to zero because both
triangle first moments vanish. The change under `ell -> ell_new` vanishes
for the same reason. The fraction is bounded by `min(|x|^2,|y|^2)` and
the logarithmic singularity is integrable for smooth cores. The expression
is invariant under a common rigid transverse rotation. A translated axis
uses the corresponding relative coordinates; the absolute coordinate
origin is not a physical director frame.

This is a complete finite-core Biot--Savart coefficient. There is no assigned
Rankine cutoff constant and no omission of interactions between cores or
between triangles.

## Exposed leading logarithm with a bounded remainder

Write each core in dimensionless local coordinates as

```
f_a(X_a+epsilon y)=Gamma epsilon^(-2) g_a(y),
integral g_a=1,  integral y g_a=0,
supp g_a subset B_Lcore,  ||g_a||_infinity <= M0.
```

The construction in 0036 supplies uniform `Lcore,M0` for small concentration.
Choose `epsilon Lcore/d <= 1/4`. In units `d=1`, the distances between
different core supports lie in `[1/2,5/2]`, and support radii from the axis
are at most `5/4`. For a self pair the coefficient of `log(d/epsilon)`
is exactly `Gamma^2 |X_a|^2`: its first moments remove all cross terms.
All remaining self logarithms and all mutual terms are bounded. One explicit
uniform logarithmic-moment bound is

```
Mlog=log(max(2 Lcore,1))+pi M0/2,
integral integral g_a(y)g_a(z) |log|y-z|| dy dz <= Mlog.
```

Consequently

```
C = rho Gamma^2 d^2/(4 pi) [3 log(d/epsilon) I_2+R],
||R||_op <= (25/16)[3(Mlog+1)+15(log(5/2)+1)].
```

The count is three self pairs per diagonal entry, six ordered distinct-core
pairs within a triangle, and nine pairs in its off-diagonal entry. It gives
a direct analytic lower bound on the smallest coefficient, separate from
the stronger exact Gram positivity. In particular the familiar quantity
`rho Gamma^2 log(d/epsilon)/(4 pi)` emerges as the singular line-tension
coefficient; it was not supplied as a premise. Its finite part is the full
integral above. Units are `C = mass length^3/time^2`, so
`C q'^2/2` is energy per axial length.

## Common bend from the same exact variation

For an axial displacement `x_perp=a_perp+U(z)`, use the same construction
with `A_e=f e`, `B_e=e.grad f`. Its exact gradient matrix at nonzero `k` is

```
C_bend,ab(k)=rho integral |fhat(xi)|^2
 [delta_ab-xi_a xi_b/|xi|^2]/(|xi|^2+k^2) d^2xi/(2 pi)^2.
```

It is positive definite. Sixfold symmetry makes it exactly
`rho <f,G_k f> I_2/2`. Common-bend/twist cross entries vanish by the
threefold rotational symmetry: such an entry would be an invariant
transverse vector under rotation by `2 pi/3`.

Unlike twist, common bend has nonzero `integral A_e=6 Gamma e`. Its local
coefficient therefore retains the physical infrared logarithm. For
`epsilon << d << 1/|k|`,

```
C_bend(k) = rho/(4 pi)
 [6 Gamma^2 log(d/epsilon)+(6 Gamma)^2 log(1/(|k|d))] I_2
 + O(rho Gamma^2).
```

A fixed outer-coherence or return-vorticity prescription is needed to turn
that bend response into a wavenumber-independent local modulus. The twist
coefficient has no such obligation. This file does not silently identify an
infrared cutoff with a fitted elastic constant.

## Claim scope

The positive result is the complete finite-core Euler kinetic-energy
second variation for the declared affine twist and bend families. Together
with 0036 it supplies microscopic angular locking and strictly positive
twist wryness from the same physical Hamiltonian. The two-triangle matrix
is retained for any subsequent collective-field transformation. This result
does not by itself construct a stationary EPS ensemble or replace the
declared affine coarse-graining premise with a full nonlinear Euler
invariant-manifold theorem.
