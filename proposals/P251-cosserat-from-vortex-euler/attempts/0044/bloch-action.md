# Exact slow angle/shape action and its isotropic continuation

## Microscopic object, map, and scope

Use the actual stationary Beltrami field from frozen attempt 0040,

```
u0=b(-sin y,2 sin x,2 cos x+cos y),
curl_physical u0=u0/ell,
```

and its exact incompressible generators

```
xi_q=(-sin y,sin x,0) cos z,
xi_s=(cos x sin y sin z,sin x cos y sin z,-2 sin x sin y cos z).
```

Coordinates have physical period `2 pi ell`, and physical displacement is
`ell xi`. The angle observable is the local rotation in the core-axis
Jacobian at the designated sections; this does not claim that the entire
noncircular cross-section moves rigidly. The section-to-section local angle
is `chi=2q`, as in 0040.

For each nonzero periodic Fourier wave `m`, construct the Coulomb potential
`A_i(m)=i m cross xi_i(m)/|m|^2`. A slow Bloch wave `kappa=ell k` uses

```
xi_i(kappa;m)=i(m+kappa) cross A_i(m),
delta u_i=P_{m+kappa}(xi_i(kappa) cross omega0).
```

Convolution with `omega0` shifts the microscopic wave before the velocity
projection. These are exact divergence-free variations, including the
return terms from differentiating the slow amplitude. The code does not
multiply an incompressible periodic generator by a slow scalar and then
ignore its new divergence.

The return has an explicit local angle map in the two designated core
jets. Explicitly,
`A_q=(sin x sin z,sin y sin z,(cos x+cos y) cos z)/2`, so the added
first-order transverse rotation is `-i kappa.grad A_q,z/2`, which vanishes
at `x=y=0`, `z=0,pi`. Differentiating the full Bloch factor also contributes
at second order: the local rotation amplitudes are
`plus/minus [1+(kappa_x^2+kappa_y^2)/2] q`. This field map is retained
when the physical core angle is needed, rather than identifying the raw
Bloch amplitude with that angle at every wave number. For the shape potential
`A_s=(-sin x cos y cos z,cos x sin y cos z,0)`, that rotation is zero
at both designated core jets, including the second-order Bloch term. Thus
the completion can translate or strain the local axis without relabeling
the shape coordinate as an independent physical angle. Changing to this
physical-angle coordinate shifts both `I2` and `h2` consistently and leaves
`h2-(h0/I0)I2` and the dispersion invariant.
For the literal difference of the two section angles, compare their slow
envelopes at the pair midpoint. The exact additional factor is
`chi=2 cos(pi kappa_z/2)[1+(kappa_x^2+kappa_y^2)/2] q_mid`. Its second-
order expansion has the axial correction `-pi^2 kappa_z^2/8`. The
independent real-space verifier checks this map and the invariance of the
normalized action coefficient under it.

The Euler orbit Hessian and KKS form, per physical volume, are

```
H_ij=rho <delta u_i^* . delta u_j
          -ell delta u_i^* . curl_physical delta u_j>,
Omega_ij=rho <omega0 . (xi_i^* cross xi_j)>.
```

Write `H=rho b^2 Hhat`, `Omega=rho ell b Omegahat`, and
`omega=(b/ell) nu`. The exact characteristic equation, with 0040's KKS
sign convention, is `det(Hhat-i nu Omegahat)=0`.

The full rational matrices are preserved in repaired-stdout.txt and are
recomputed by bloch_sector.py. The following are their exact expansions,
not values from a numerical wave-number fit.

## All raw action entries through second order

For a wave along each of the three coordinate axes,

```
Hhat=diag(h0+h2 kappa^2,s0+s2 kappa^2)+O(kappa^4),
Omegahat=[[0,B0+B2 kappa^2],[-B0-B2 kappa^2,0]],
h0=7/48, s0=25/96, B0=-1/4.
```

The exact matrices have no first derivative along any axis. Since such a
derivative is a linear form in the wave vector, the entire chiral first-
derivative tensor vanishes in this chosen standing-section sector. It was
computed, not discarded by a parity average. The energy off-diagonal and
the KKS diagonal entries vanish along all three axes as well.

| Axis | `h2` | `s2` | `B2` |
|---|---:|---:|---:|
| x | `35693/54000` | `-553/1728` | `1/16` |
| y | `171293/216000` | `653/1728` | `-1/8` |
| z | `13327/43200` | `-1495/1728` | `0` |

These raw entries are the inputs to the parent's common-angle/body lift.
Adding another same-fluid inertia and then reusing the normalized dispersion
below unchanged would be incorrect.

## Elimination retains both varying forms

For this two-coordinate sector alone, eliminate the conjugate shape using
the complete action. Its angular inertia and potential are

```
Ihat(kappa)=B(kappa)^2/s(kappa)=I0+I2 kappa^2+O(kappa^4),
I0=6/25,
I2=I0[2 B2/B0-s2/s0],
Khat(kappa)=h0+h2 kappa^2+O(kappa^4).
```

The exact coefficients are

| Axis | `I2` | `C=h2-(h0/I0)I2` | coefficient of `kappa^2` in `nu^2` |
|---|---:|---:|---:|
| x | `328/1875` | `29953/54000` | `29953/12960` |
| y | `-203/1875` | `185503/216000` | `185503/51840` |
| z | `299/375` | `-7603/43200` | `-7603/10368` |

Here `nu^2(0)=175/288`. The local collective redefinition
`q(kappa)=[1-I2 kappa^2/(2I0)] Q(kappa)+O(kappa^4)` makes its kinetic
term constant and gives the displayed `C`. It agrees with the physical
core-angle coordinate at zero wave number. Equivalently, one can retain
the raw derivative inertia and derive the same dispersion directly.

The physical gradient stiffness is `rho b^2 ell^2 C`; the derivative
inertia coefficient is `rho ell^4 I2`. Thus the negative `C_z` is a
genuine oriented-cell optical curvature, not a unit mismatch or a missed
shape elimination. The full zero-wave-number Hessian remains strictly
positive, so the long-wave branch is stable near its nonzero gap; a negative
curvature is not a claim of immediate instability.

## The original isotropic ensemble supplies the positive continuation

The parent target is an isotropic affine continuum, not a medium in which
every tube axis is parallel. Declare a uniform SO(3) ensemble of the cell
triad, and assign its scalar angle amplitude by `q=n.Q`, where `n` is its
tube axis. The local angle map and this ensemble premise are separate.
The shape is eliminated within the stated internal sector; a later body
constraint must instead use the raw matrices above.

Average the action tensors, not the frequencies. For either the raw
potential-gradient tensor or derivative-inertia tensor, the SO(3) average
depends only on its transverse trace and axial component. Possible mixed
body-frame curvature entries do not survive this average. Its two tensor
invariants are determined by the three axis evaluations above. The common
zero-order inertia and potential average to `I0 I_3/3` and `h0 I_3/3`.
Because these are scalar multiples of the identity, normalizing the averaged
kinetic operator to order two gives exactly the average of
`h2-(h0/I0)I2`. This proves why the displayed averaged coefficient follows
from the action and does not assume commutation of an arbitrary frequency
average with momentum elimination.

Put `Ct=(C_x+C_y)/2=61063/86400`, `Cz=-7603/43200`, and
`G_ij=partial_j Q_i`. Exact spherical fourth moments give the gradient
quadratic form, before an overall physical factor `rho b^2 ell^2/2`,

```
A ||G||^2+B[(tr G)^2+tr(G^2)],
A=(4Ct+Cz)/15=114523/648000,
B=(Cz-Ct)/15=-25423/432000.
```

This density is not pointwise positive on all unconstrained matrices `G`;
in particular its pure-trace coefficient records the negative oriented axial
curvature. That fact is preserved, not hidden by an isotropy claim.

On a periodic macroscopic domain, or for compactly supported macrofields,
the density is equivalent as an action to the strictly positive form

```
A |curl Q|^2+(A+2B)|div Q|^2,
A=114523/648000 > 0,
A+2B=19127/324000 > 0.
```

The exact difference is `(A+B) N`, where

```
N=tr(G^2)-(tr G)^2
 =div F,
F_j=sum_i [Q_i partial_i Q_j-Q_j partial_i Q_i],
A+B=152777/1296000.
```

Thus its integral vanishes under the stated periodic or compact-support
boundary condition. With a physical boundary the displayed flux must be
retained; this calculation does not silently change boundary tractions.
The positive curl/divergence action is coercive for periodic mean-zero
fields and produces positive transverse and longitudinal dispersion. With
the averaged inertia `rho ell^2 I0/3`, their squared-speed coefficients are

```
c_T^2=b^2 (114523/51840),
c_L^2=b^2 (19127/25920).
```

They describe this internal sector before its common-angle/body coupling.
If the collective variable is instead the relative local angle, its zero-
wave-number map is `chi=2q`; all zero-order action coefficients divide by
four together. Retaining the finite-section derivative map above shifts
raw derivative terms consistently and leaves the normalized speeds unchanged.

## Verdicts and continuation boundary

The parallel-oriented-cell route to positive axial optical curvature is
refuted for these exact generators, with its negative coefficient named
above. The original isotropic, periodic affine action route is established:
the complete same-fluid angle/shape sector gives a positive gap and a
strictly positive curl/divergence gradient action, with its boundary term
and local field map explicit. This is a representation/ensemble continuation,
not a fitted change to `a/b` or deletion of a chiral term.

The positive statement does not close the parent campaign by itself. The
common-angle/body sector and any material displacement coupling must be
combined at the raw Hessian/KKS level. The physical angle remains a local
core-axis observable, not a rigid rotation of the entire tube section.
Stationary EPS geometry and claim promotion are parent obligations.
