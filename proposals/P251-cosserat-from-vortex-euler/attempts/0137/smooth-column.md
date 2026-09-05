# An actual smooth-column surface pole and its physical moment match

## 1. Statement and primary-source boundary

Fix m>=2, core radius a and core angular speed O>0. For every sufficiently
small fixed axial carrier k*>0, there are smooth compactly supported
vorticity profiles, exactly 2O on r<=a, whose actual linear Euler velocity
has a simple real surface-mode branch near the Rankine branch on a
neighborhood of k*. Its intrinsic squared frequency has positive first
derivative with respect to k² and positive second derivative with respect
to k. Its rotating-material KKS scalar action is positive. A nonnegative
material marker supported strictly inside r<a can be chosen so that its
full measured spin equals the canonical angular momentum at k*.

This is a smooth ordinary column in R² times an axial periodic direction
(or a fixed axial Fourier fiber), not yet an EPS knotted field. The
eigenmode has finite transverse energy and exponentially decaying exterior
velocity. The statement is linear Euler, not an exact finite-amplitude
nonlinear normal mode. The physical marker match is made at the fixed
carrier; its zero-carrier singularity from0135 remains explicit.

Primary comparison: Gallay--Smets, *Spectral stability of inviscid columnar
vortices*, arXiv:1805.05064v3, Remark1.4 and §6.2, equations(6.10)--(6.15),
<https://arxiv.org/pdf/1805.05064>. Their Rankine matching calculation is
consistent with the imported determinant. Their general spectral theorem
has its own profile hypotheses; it is not a smoothing theorem used here.
Their distinction between velocity and vorticity function spaces at the
patch interface is retained. The proof below is a direct smooth radial
transfer argument, not an assertion of operator-norm convergence of the
vorticity generators across a discontinuity.

## 2. A bounded-coefficient radial system derived from full Euler

Allow initially the more general stationary column
`u=r O(r)e_theta+W(r)e_z`, and write

```
Z(r)=2O+r O',    s(r)=omega-m O(r)-k W(r).
```

Z is axial vorticity, while p below is perturbed pressure divided by
density. In the `exp(i m theta+i k z-i omega t)` sector the three exact
momentum equations and incompressibility are

```
-i s v_r-2O v_theta = -p',
-i s v_theta+Z v_r = -i m p/r,
-i s v_z+W' v_r = -i k p,
v_r'+v_r/r+i m v_theta/r+i k v_z=0.
```

Set `f=v_r/(-i s)`, which is the radial material displacement. Then

```
v_r=-i s f,   v_theta=m p/(r s)-Z f,   v_z=k p/s-W' f,

f' = -[1/r+2m O/(r s)] f +(m²/r²+k²)p/s²,
p' = (s²-2O Z) f +2m O p/(r s).                 (1)
```

The cancellations in the first equation use
`s'=-m O'-k W'` and `Z=2O+r O'`. In particular (1) contains no
derivative of Z, and no W' outside its contribution to s. This matters
because Z' grows when a vortex patch is smoothed. A pressure-only
second-order equation can hide this bounded transfer structure.

The actual material displacement in this Fourier sector is

```
eta_r=f,
eta_theta=i[m p/(r s)-2O f]/s,
eta_z=i k p/s².
```

These expressions satisfy the full Lin relation in the region s!=0;
they follow by including `r O' eta_r` and `W' eta_r` in the tangential
and axial displacement/velocity relations. Thus both material and
pressure rows remain in the construction.

## 3. Smooth profile, nonresonance, and exact exterior matching

Choose once a smooth nonincreasing cutoff c(t), equal to 1 for t<=0,
zero for t>=1, and between 0 and 1 otherwise. For e>0 set

```
Z_e(r)=2O c((r-a)/e),
O_e(r)=r^-2 integral_0^r s Z_e(s) ds.
```

Then the velocity `r O_e e_theta` is smooth everywhere, equals rigid
rotation on r<=a, and has compact vorticity supported in r<=b=a+e.
Outside b, `O_e=C_e/r²`, where
`C_e=O a²+O(e)` with dimensionful constants fixed by a,O. This is an
actual stationary smooth Euler solution with `p0'=r O_e²`.

Fix a small nonzero k* and a compact k-interval around it. The Rankine
surface frequency satisfies `omega_R(k)=O(m-1)+O(k²a² O)` and is simple.
Choose the interval and a fixed omega-neighborhood so that throughout
the annulus, for all sufficiently small e,

```
|omega-m O_e(r)| >= O/2,    a<=r<=a+e.          (2)
```

This follows from the uniform estimate `O_e(r)=O+O(e/a)` there.
The regular interior pressure is exactly the Bessel core solution from
0135. Let `Y_in(omega,k)=(f(a),p(a))`, with the same nonzero pressure
normalization, and let `T_e` be the fundamental transfer matrix of (1)
from a to b for W=0.

Exterior vorticity is exactly zero. With `K(r)=K_m(k r)` and
`s_b=omega-m C_e/b²`, the decaying irrotational response has boundary
vector proportional to

```
Y_out(omega,k,e)=(K'(b), s_b² K(b)).            (3)
```

Here K' differentiates r, not its argument. For example take the
exterior velocity potential `-i s_b K(r)`; its pressure is
`s(r) s_b K(r)`. This is smooth even at the exterior particle critical
radius where s(r)=0. Its material displacement may have a resonant
label component there, but its velocity and pressure do not.

The complete matching function is the real determinant

```
E_e(omega,k)=det(T_e Y_in, Y_out).              (4)
```

At e=0 it is a nonzero smooth multiple of the canonical Rankine
determinant. Both physical pressure and radial material displacement
are matched. This construction does not replace exterior pressure by
a wall or by a local kernel.

## 4. Quantitative transfer and the genuine mode

On the compact annulus, (2) bounds the coefficient matrix in (1) and
its omega/k derivatives of any fixed finite order uniformly in e.
Z_e is bounded by 2O. All derivatives here are parameter derivatives,
not radial derivatives of the thin layer. Peano--Baker iteration gives

```
||T_e-I|| <= exp(M e)-1,
||partial_parameters^alpha(T_e-I)|| <= C_alpha e, |alpha|<=3.
```

Differentiate the integral equation for T_e to obtain the second bound;
each differentiated term still contains one integral across a length-e
interval. The exterior vector (3) converges in the same parameter norms
because b->a, C_e->Oa² and k stays separated from zero. Consequently

```
||E_e-E_0||_(C^3(omega,k)) <= C e.             (5)
```

The Rankine root has nonzero omega derivative. Shrink the compact
parameter neighborhood once, so this derivative stays bounded away
from zero. The ordinary implicit-function theorem, or strict real
monotonicity in omega followed by differentiation of (4), gives a
unique real root omega_e(k) there, and

```
||omega_e-omega_R||_(C²(k)) <= C' e.           (6)
```

The core solution, the annular solution of (1), and the exact exterior
potential now glue to an actual global smooth Euler eigenvelocity.
At b the base vorticity and all its derivatives vanish; continuation
from the irrotational solution therefore gives smooth matching, not a
new sheet. The regular core solution gives smoothness at the axis.
Exponential decay of K_m(k r) gives finite transverse energy.

This is an eigenmode on the dynamically invariant class with perturbed
vorticity supported in r<=b. Indeed the exterior background has zero
vorticity, its perturbation is zero there, and the linearized vorticity
equation preserves that property. A compact divergence-free generator
can realize its isovortical tangent: use the displayed material
displacement throughout the vorticity support and continue it through
a nonresonant exterior collar before the critical radius, then take a
smooth divergence-free extension. The azimuthal normal flux has zero
integral. Exterior modifications act only on relabeling and create no
vorticity. No spectral gap against *all* exterior-vorticity disturbances
is asserted; the root of the stated matching problem is simple.

## 5. Positive actual action, physical curvature, and matched tag

Put `sigma_e(k)=omega_e(k)-m O`. At fixed small k*, (6) preserves
sigma_e<0 and both positive curvature quantities from0135:

```
partial_(k²) sigma_e² > 0,   partial_k² sigma_e² > 0.
```

In the second expression `partial_k²` denotes the second derivative
with respect to k, not differentiation with respect to k². Both are
positive on a sufficiently small nonzero-carrier interval before e is
chosen. This does not discard their generally nonzero group derivative.

For the two real mode columns, define the complete KKS scalar

```
beta_e=2 pi rho Lz integral_0^b Z_e(r) f(r) B(r) r dr,
B=[m p/(r s)-2O_e f]/s.
```

The core integrand converges to the Rankine one. The annular integral
is O(e), because f,B,Z_e remain bounded there. Thus beta_e converges
to the strictly positive Rankine beta, and remains positive. The
actual SO(2) momentum and the mode evolution give exactly
`H_rot''=-beta_e sigma_e I`. For a material-angle row c_e the scalar
inertia is `M_e=-beta_e/(sigma_e c_e²)>0`. Its derivative terms are
retained: the normal-form spatial coefficient is
`M_e partial_(k²)(sigma_e²)`, not a bare Hessian entry.

Place the two signed radial bumps and the tag strictly inside r<a,
where the actual background remains *exactly* rigid rotation. The
complete shape/spin formulas of0135 hold exactly with the new pressure,
frequency and beta_e. Since k*!=0 and the radial Bessel parameter is
nonzero, the two moment rows Q,T remain independent. Their prescribed
ratio, computed with beta_e and sigma_e, gives exactly `j_tag=M_e>0`
at k*. Both radial bump coefficients can be scaled together to maintain
nonnegative tag density. The reference material shape is nonzero.
There is no dephasing inside this ordinary-column core.

For transfer from the patch marker, first fix k* and its two bump
intervals, then choose e smaller than the associated nonzero moment
matrix gap. The gap is O((k*a)²); the order of choices is essential.
Alternatively solve the two moment equations directly for the smooth
column's exact Bessel profile. No measured number or frequency fit
enters either construction.

## 6. Explicit force-free interface and remaining assembly

Equation(1) also applies exactly to the parent's generalized force-free
column with axial W. After subtracting its uniform axial velocity U0,
the only changed coefficient is `s=omega-m O_e-k(W-U0)`.
If W-U0 tends to zero in the required fixed parameter norms, the same
Evans transfer/IFT argument continues the actual pole. Thus no missing
W' term or omitted exterior pressure obstructs the radial mode transfer.

The material observation is more delicate: W varies with r, so a
transported tag samples phases `exp[-i(omega-m O-k W(r))t]`. Its finite
radial moments need not be an exact single-frequency scalar. Retain the
bound `|exp(i k (W(r)-W_ref)t)-1|<=|k t| osc_tag W` in every actual
shape/spin integral; the small shape denominator and moment gap enter
the resulting relative error. This is a finite-time observation error,
not an exact equality to be inferred from a close eigenfrequency.
The parent owns that force-free/local-torus join.

Route verdict: **established as stated** for the actual smooth ordinary
column mode, positive action/curvature and fixed-carrier physical
moment match. The finite-time fallback was unnecessary for this pole.
This result activates, rather than proves, coherent translation and
counterpropagating-carrier assembly, plus the separate force-free and
global geometric transfer. No accepted claim or canonical API changes.
