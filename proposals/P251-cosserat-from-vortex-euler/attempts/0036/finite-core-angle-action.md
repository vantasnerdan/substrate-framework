# Exact affine angle action of a smooth six-core Euler molecule

The base state is the exact smooth six-core rotating equilibrium constructed
in smooth-polygon-construction.md. This file derives the pullback of the Euler
Hamiltonian action onto a specified finite-dimensional, isovortical trial
family. The pullback is exact at nonzero core size. Its use as a coarse-grained
dynamics requires the declared affine collective-motion premise; it is not
a proof that the trial family is an invariant manifold of unrestricted Euler.

## Admissible trial variations

Divide the six cores into two alternating triangles. In triangle `i`, move
the three centroids to radius `r_i` and angle `theta_i+2 pi j/3`; rigidly
rotate each core profile by the same angle change. Core size, circulation,
and the entire vorticity distribution within each core remain unchanged.
Supports stay disjoint in a neighborhood of the regular hexagon. Each motion
extends to a smooth area-preserving diffeomorphism of the plane: use local
linear/quadratic streamfunctions equal to the rigid-motion generators on
the supports and smoothly cut off outside larger disjoint disks. Thus these
are actual tangent directions of the Euler vorticity orbit, not a new
director assigned an independent stored energy.

Let `C_epsilon=integral |x-X|^2 omega_core(x) dx`, with the centroid condition
`integral (x-X) omega_core=0`. Pullback of the Euler orbit's symplectic form
under these rigid motions gives the translation circulation form and a
constant intrinsic angular-impulse term. With the conventions of attempt
0032, a potential is

```
sum_{i=1}^2 [-3 rho Gamma r_i^2/2-3 rho C_epsilon/2] d theta_i.
```

The `C_epsilon` term is locally an exact differential; it contributes no
quadratic dynamics. Thus the centroid canonical momenta are exactly
`p_i=-3 rho Gamma r_i^2/2`, not a point-core approximation. The translation-
rotation cross term vanishes because each core's first moment is zero.
This can also be checked directly from the Euler orbit form on the affine
streamfunctions: its integrals depend only on circulation, first moment,
and angular impulse. Signs agree with the point-vortex convention in 0032.

## Energy and its exposing limit

Use the renormalized Euler kinetic Hamiltonian per unit axial length. All
self energies are constant on the trial family. The variable part is exactly

```
H_epsilon = -rho/(2 pi) sum_{a<b}
  integral integral omega_a(x) omega_b(y) log|x-y| dx dy.
```

It is finite and smooth for disjoint supports. Let `epsilon` be the core
scale and `d` the equilibrium centroid radius. Taylor expansion of each
cross kernel, including two derivatives in the dimensionless collective
coordinates, is uniform on a fixed neighborhood of the regular hexagon.
The first-order term vanishes by the centroid constraints. Therefore

```
||H_epsilon-H_point||_{C^2(x,chi)}
  <= C rho Gamma^2 epsilon^2/d^2,
```

where `x=(r_1^2-r_2^2)/(r_1^2+r_2^2)`, `chi=theta_1-theta_2`, and the
constant depends on a bounded rescaled support and the fixed neighborhood,
not on a fitted target. Derivatives of the kernels are bounded because every
inter-core distance there is bounded below independently of `epsilon`.
This is a direct integral estimate, not an inference from velocity Holder
convergence or a discretized soft eigenvalue.

The exact base state's rotating Euler variational equation ensures that
`H_epsilon` is critical on the fixed-total-angular-impulse subfamily.
Since the intrinsic core impulse is constant, this constraint is precisely
fixed `S=r_1^2+r_2^2`. Reflection and triangle-exchange symmetry make the
mixed `x,chi` Hessian vanish at `x=0, chi=pi/3`. Attempt 0032 supplies
the independent exact point limit

```
A=3 rho Gamma^2/(4 pi),
(H_point)_xx=A/2,
(H_point)_{chi chi}=9 A/2.
```

The uniform `C^2` estimate therefore proves that both corresponding exact
finite-core Hessian entries are positive for all sufficiently small nonzero
`epsilon`. No finite maximum radius is asserted without evaluating the
explicit remainder bound.

## The finite-core angle oscillator

Define `P0=3 rho Gamma S/4`, so the relative momentum is `p=-P0 x`, and
write `q=chi-pi/3`. The quadratic restricted Euler action is exactly

```
L2 = p qdot - (H_epsilon)_xx p^2/(2 P0^2)
              - (H_epsilon)_{chi chi} q^2/2.
```

Eliminating its conjugate momentum gives

```
L2 = I_epsilon qdot^2/2-K_epsilon q^2/2,
I_epsilon=P0^2/(H_epsilon)_xx > 0,
K_epsilon=(H_epsilon)_{chi chi} > 0.
```

These coefficients are derived from convergent, finite smooth-core Euler
integrals. Their point limits are the coefficients in 0032, but at finite
core size they are not replaced by those limits. The action contains reaction
of both triangles and no prescribed strain source.

This establishes an exact finite-core Euler-derived affine angle action.
It does not establish a full-PDE isolated optical eigenmode: that different
statement requires a coupled transport/resolvent analysis. Nor does it by
itself identify a triangle with a macroscopic affine cage, derive three-
dimensional wryness, or ensure stationary EPS compatibility. Those remain
distinct parent constructions. The original affine coarse-graining premise
may license the specified variational restriction, but that premise must be
stated where the continuum action is used.
