# Optical histories on the same compensated stationary Euler array

## 1. One background, and the result actually constructed

Use a triangular lattice Lambda with primitive area A and spacing L.
The total physical planar vorticity is
`zeta=q_core-Gamma/A`, not q_core alone. The torus Green function obeys
`Delta G=delta-1/A`, and planar velocity has zero harmonic mean.
The positive excess core has circulation Gamma and characteristic
radius a. All periodic copies belong to this one stationary field.

The result below constructs Kelvin-prepared optical histories on this
same smooth array and its exact Bernoulli force-free lift. Their local
material angle, spin and action approximate the isolated smooth-column
ones in the first two axial-carrier derivatives over a fixed optical
window. The actual transverse Bloch response is retained with a bound;
its sign is not inferred. This is a genuinely uniaxial result, not an
isotropic constitutive theorem or independent-cell superposition.

Fix a small nonzero dimensionless axial carrier q*=k*a, a compact
interval I around it, and an optical window |Omega t|<=T*. The smooth
core taper and tag moment gaps are fixed first. Let

```
ell=L/a, epsilon=1/ell, u_axial=U0/(Omega a).
```

Dilution ell>>1 and large u_axial are then independent choices. Below,
spatial and temporal estimates are stated in core units a and Omega^-1.
Constants may depend on I, T*, the fixed smooth profile and the fixed
nonzero tag-moment determinants, but not on ell or u_axial sufficiently
large. Keeping q* separated from zero is essential for the pressure
locality argument.

## 2. The same optical core is an admissible stationary array seed

0139 uses the particular0036 radial seed. To join0137/0140, extend its
bordered Green construction to the thin smooth Rankine taper used there.
Let Z(r) be exactly 2Omega on r<a, smoothly decreasing to zero across a
thin annulus of thickness eta*a, with eta fixed and sufficiently small.
Let V=r^-1 integral_0^r sZ(s)ds. Its streamfunction is strictly monotone
away from the origin, so Z is a smooth function F of that streamfunction;
F' is supported on the transition annulus. In the sign convention of
0136 the linearized radial potential is `Qpot=-Z'/V>=0`.

The m=1 regular homogeneous solution is V. Its positive ground-state
identity proves that only the translation kernel remains in m=1;
adding `(m²-1)/r²` removes every m>=2 kernel. Sixfold invariance excludes
translations. It remains to check the m=0 mass border rather than
assuming the0036 seed's radial result transfers.

Normalize the regular radial homogeneous solution Y to equal 1 in the
constant-vorticity interior. Across the thin annulus it obeys

```
(rY')'=-r Qpot Y.
```

The integral of r Qpot is 2+O(eta), since V=Omega a[1+O(eta)] there
and the vorticity drops by exactly 2Omega. Its L1 norm is uniformly
bounded. The integral equation gives Y=1+O(eta) across the annulus,
even though Qpot itself is large. Hence its exterior form is

```
Y=A_log+B_log log r,   B_log=-2+O(eta) !=0.       (1)
```

The zero-mass condition in the bordered Green kernel is precisely zero
exterior logarithmic coefficient. It therefore removes the sole regular
m=0 homogeneous direction. The constant multiplier is then fixed by
the Green equation. This proves the needed bordered invertibility for
the selected optical core, without a new assumed general-seed theorem.

Now apply0139's exact periodic Green/IFT construction with this seed.
Its regular part has zero gradient at the lattice point and bounded
second derivatives of order L^-2. The resulting smooth core profile
and local velocity differ from the selected radial core by O(epsilon²)
in each fixed rescaled regularity norm. The physical compensation
`-Gamma/A` is part of that estimate, not removed afterward. The exact
array remains stationary because its excess vorticity is a function
of its actual streamfunction. The core's sixfold deformation and its
material streamline changes are retained in the comparison.

## 3. Exact Bernoulli lift and Kelvin preparation on that field

For this one stationary planar solution v, set
`B=p0+|v|²/2` with one fixed additive pressure registration. Its exact
Bernoulli lift is

```
u=(v,W), W=sqrt(2(C-B)), C>max B.
```

Because v.grad B=0, this is stationary three-dimensional Euler.
With v=J grad psi and grad B=zeta grad psi, direct differentiation
gives `curl u=(zeta/W)u`. It is a single smooth periodic generalized
force-free field, with all compensation and ambient fluid included.
In an axial Galilean frame use h=W-U0, choosing U0 by the same fixed
Bernoulli registration. At large u_axial the normalized h and every
needed fixed derivative are O(u_axial^-1). Neither constant curl nor
global knotted EPS geometry follows from this lift.

Use0137's compact divergence-free isovortical generators for the two
real optical mode columns, supported on the core and one fixed nearby
collar. Apply these generators to the ACTUAL array vorticity, not to a
separately inserted isolated column. Their induced initial velocity is
the complete periodic/Bloch Leray projection of `xi cross omega_array`.
It is an exact tangent to the array's Kelvin leaf. The constant
compensation can contribute to a three-dimensional optical KKS pairing
and is retained; the planar-translation Poisson-bracket cancellation
from0139 is not invoked for these different generators.

For transverse Bloch vector Q, copies of this same initial generator
have the actual phases exp(i Q.R). The linear Euler evolution on the
array, and its material Lin transport, are then solved without restricting
later core shapes or resetting any harmonic mean. Thus inter-core and
ambient reactions remain in the actual history.

## 4. Full pressure has a massive, exponentially local kernel

Fix an axial Fourier k in I/a. Write
`D_k=(partial_x,partial_y,i k)` and
`G_k=(-Delta_perp+k²)^-1`. For an actual divergence-free perturbation w,
the exact linear Euler pressure satisfies

```
(-Delta_perp+k²)p = 2 D_j[(partial_i u_j)w_i].       (2)
```

Only horizontal background derivatives occur; the axial multiplication
in D_j includes the axial shear term. The pressure force is therefore
the full order-zero operator
`-2 D_l G_k D_j[(partial_i u_j)w_i]`. It is not replaced by a local
model. On the plane `G_k(r)=K0(|k|r)/(2pi)`; the periodic Bloch kernel
is its exact image sum with the Bloch phases.

For a horizontal exponential weight exp(alpha.x), with |alpha|<k_min,
conjugation shifts the Fourier variable by i alpha. The denominator has
real part `|p|²+k²-|alpha|²`, strictly positive. Consequently the
conjugated operators `D_l G_k D_j` and their first two k derivatives
are bounded on L² uniformly for k in I/a. One may fix
|alpha|=k_min/2 before taking the lattice limit.

The transport part contributes at most
`|alpha| ||u_perp||_infinity` to weighted energy growth; its ik h term
is imaginary multiplication. The strain and pressure terms are bounded
by `C(k_min,alpha)||grad u||_infinity`. Gronwall therefore gives the
weighted actual Euler propagator estimate

```
||exp(alpha.x) S(t,k) exp(-alpha.x)|| <= exp(C|t|). (3)
```

The same argument applies after two k derivatives, using Duhamel and
the bounded differentiated pressure operators. High fixed Sobolev
versions follow by commuting derivatives through the smooth background,
with the fixed-profile bounds. The large constant U0 is absent because
the axial Galilean frame was taken first. This proves finite-time
off-cell locality for the COMPLETE pressure evolution, not a finite
speed assertion for Euler.

## 5. Actual transverse Bloch terms are retained and bounded

Observe one transported tag in the central core. Decompose its response
to Bloch-phased initial generators into responses from each translated
source core R. Translation covariance of the SINGLE periodic background
and (3) give, for fixed compact observation/core supports,

```
local_response(Q,k,t)=sum_R exp(i Q.R) response_R(k,t),
||response_R||_(C_k², fixed time/Sobolev rows)
    <= C_T exp(-theta |R|), R!=0, theta<k_min.      (4)
```

The sum is the full Euler response, including repeated interactions.
Differentiating the actual phases gives the transverse Bloch jets. In
core units, for j=1,2,

```
||partial_(aQ)^j [local_response(Q)-response_0]||
    <= C_T ell^j exp(-c q_min ell).               (5)
```

The zero-th image contribution has the analogous exponential bound.
This retains the actual transverse coefficient; its sign could be
either sign and has not been selected from a trial energy. It is small
relative to the order-one axial optical second jet at large ell.
Resolving a positive transverse band curvature would require a
separate sign calculation of these interaction terms, not an absolute
O(epsilon²) isolated-core estimate.

Equation(5) is a LOCAL tagged-response bound, not a bound on every
physical coarse mean. The Fourier mean of a Bloch amplitude contains
the cell weight exp(-i Q.x), whose derivatives give algebraic physical
cell moments even when (5) is exponentially small. Those mean/current
rows must be retained explicitly in the parent's joining calculation.
Also this z-independent background preserves the true axial Fourier
sector: an optical carrier k* is not silently turned into point-Euler
axial mean k=0. A material/hybrid observation can mix those scales only
through its actually defined tag/current map.

## 6. Same-array axial optical comparison in C², not absolute-error inference

Compare a central prepared history to the actual smooth isolated-column
history of0137/0140 with the same normalized core and generators.
The image pressure difference on compact supports is exponentially
small by (4). In a central neighborhood, the smooth stationary array
background differs by O(epsilon²), with all required derivatives.
More generally its regular external field is bounded by
`C epsilon²(1+|x|)` in core units before reaching a neighboring core.
This polynomial is integrable against the exponentially localized
column data/propagator. Outside that neighborhood use (3) and its
exponential tail. Full Duhamel, including the actual projector and
background differences, gives

```
||history_array-history_column||_(C_k², local H^s, |t|<=T)
 <= C_T [epsilon²+u_axial^-1+exp(-c q_min ell)] = C_T d.   (6)
```

Parameter derivatives include the initial Leray projection, the exact
Euler evolution and the same material Lin transport. The high Sobolev
version transfers the core-angle jet and genuine moving-tag integrals.
The reference marker's moment gap is fixed and nonzero, so the actual
array marker's angle chart remains nonsingular. The actual steady array
tag is transported by that array, not rigidly rotated by hand; its
background deviation is included in (6).

The KKS form is computed on the actual generators and actual base
vorticity. It differs from the positive column form by O(d), including
the compensation and W' contributions. Actual Euler linear evolution
preserves this form. On the prepared real two-column family, the actual
material-angle observation gives a nonsingular time-dependent row.
Pulling the same phase action through this row and eliminating its
conjugate quadrature is the exact construction of0140. Its coefficients
and the full mechanical spin/current rows inherit the C² error bound
(6), including temporal connections. This is a time-dependent
prepared-solution construction, not an invariant autonomous oscillator
manifold for arbitrary perturbations. The ordinary optical positivity
survives at this stated finite-window/action-jet scope.

To expose the natural axial-gradient scale, take the symmetric physical
carrier difference on THIS same array,

```
Delta_K history=history(k*+K)+history(k*-K)-2history(k*).
```

The integral Taylor remainder applied to (6) proves

```
||Delta_K history_array-Delta_K history_column||
 <= C_T d K².                                    (7)
```

The corresponding physical action/current coefficient errors are also
O(d), not an unscaled absolute residual. Choose dilution and axial
velocity so that C_T d is smaller than the fixed positive axial
propagation margin. This is the differentiated natural-gradient-scale
comparison needed by0138. It does not assert that the absolute
on-site history error O(d) is o(K²) at fixed d. An autonomous model
would have to retain or separately identify its on-site and first-jet
corrections. The same distinction applies to the finite-window group
error already kept in0138.

## 7. Scope, physical mean and activated next route

Established result: a same-background stationary smooth periodic Euler
array and exact Bernoulli lift, with actual Kelvin-prepared optical
histories, positive column-like local material action/axial gradient
response at controlled C² finite-window accuracy, full pressure and
inter-core reactions, and bounded actual transverse Bloch terms.
No independent array copies with different orientations are averaged
into an asserted isotropic medium.

0139's actual planar acoustic histories are separately available on
this same array seed through its unchanged dynamic localization proof.
The Bernoulli lift preserves the planar subsystem at its proved scope.
This common stationary background is genuine joining progress, but
nonzero mixed acoustic/optical coefficients and the physical hybrid
mean rows are not determined merely by coexistence of these histories.
The parent and0144 own that full mean/observation join.

Route verdict: **established as stated** for the controlled same-array
anisotropic finite-window response. A positive transverse optical band,
three-dimensional isotropic same-field architecture, global stationary
knotted geometry and autonomous coupled continuum remain distinct next
constructions. None is ruled out by the present bound. No numerical
soft-mode calculation, canonical change or new claim review was used.
