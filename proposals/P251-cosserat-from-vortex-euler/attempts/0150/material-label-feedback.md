# Actual finite-k core/ambient feedback through an advected material label

The source conventions are0144 and0146. Fix one primitive periodic cell
of area `A_c`, positive smooth excess core `Q0`, circulation
`Gamma=integral Q0`, total vorticity `zeta=Q0-gamma`,
`gamma=Gamma/A_c`, and velocity `v=J grad psi` with
`Delta psi=zeta`, `<v>=0`. All operators below use this same periodic
Green function. The actual smooth stationary field is0139/0141's
compensated core, not a superposition of isolated velocities.

## 1. An exact finite-wave core variable

Consider a full Euler perturbation with axial Fourier factor `exp(ikz)`.
Write its horizontal velocity `w`, vertical velocity `b`, and vertical
material displacement `chi`. Because the stationary base has no axial
velocity, the exact linearized label equation is

```
(partial_t+A)chi=b,       A=v.grad.
```

The physical vertical vorticity perturbation is `eta=curl_h w`.
Its actual stretching equation is

```
(partial_t+A)eta+w.grad zeta=ik zeta b.
```

Since `A zeta=0`, the combination

```
q=eta-ik zeta chi
```

therefore satisfies exactly

```
(partial_t+A)q=-w.grad zeta.                              (1)
```

If q is initially supported in the core, it stays there. The physical
vorticity eta does not: its ambient part is `ik zeta chi`. Thus (1)
is an exact representation change retaining stretching, not an
assumption that ambient Euler vorticity is frozen. It is the linear
potential-vorticity identity for the actual advected vertical label.

Let `G=(-Delta_h)^(-1)` on mean-zero functions, `K=-J grad G`,
and `d=grad G`. Products below have their mean removed when G is
applied. The divergence constraint fixes the complete velocity as

```
w=m+K q+ik K(zeta chi)+ik d b.                            (2)
```

The zero-average vertical vorticity constraint
`integral(q+ik zeta chi)=0` is retained; it is not reset separately.
The physical mean is `m=<w>`, and the exact mean equation remains
`m_dot=-ik <v b>`.

Define the actual pressure map from0146 by `Q_k`, and put

```
L_c q=-A q-(K q).grad zeta,       T m=-m.grad zeta,
B chi=-(K(zeta chi)).grad zeta,   C b=-(d b).grad zeta,
P_k=Q_k K,     R_k=Q_k K zeta,   S_k=Q_k d.
```

The full Euler/label system is

```
q_dot=L_c q+T m+ik B chi+ik C b,
chi_dot=-A chi+b,
b_dot=-A b-ik P_k q+k^2 R_k chi+k^2 S_k b,
m_dot=-ik <v b>.                                        (3)
```

The initial actual common velocity data have
`q(0)=chi(0)=b(0)=0`, `m(0)=V0`. No fixed Kelvin graph or isolated
reaction inverse is appended. The passive label evolves with the same
Euler perturbation, so the stretching and reconstruction constraint
are present at every time.

## 2. The core translation row is an exact complete-fluid current

Use ordinary coordinates x on a neighborhood containing the compact
core, and define its linear centroid row

```
ell q=(1/Gamma) integral x q.
```

The compact support makes this independent of an arbitrary coordinate
cut outside that neighborhood. Antisymmetry of the actual periodic
Biot--Savart kernel gives

```
ell T=I,       ell L_c=0,
ell B chi=-(1/Gamma) integral zeta v chi,
ell C b=-(1/Gamma) J integral v b.                        (4)
```

For example, integrate (1) against x, then use
`integral Q0 K(zeta chi)=-integral zeta chi K Q0`
and `K Q0=v`. Similarly
`integral Q0 grad G b=-integral grad G Q0 b=-J integral v b`.
These identities retain the uniform compensation and the harmonic
mean convention. The negative circulation background is not a second
independently translated fluid.

If the actual bounded transport coordinates satisfy `A r=v`, then
`A(zeta r)=zeta v`. The first ambient row in (4) consequently has the
explicit bounded primitive `zeta r/Gamma`; no abstract forcing-resonance
exclusion or separatrix frequency gap is needed for this row. In
particular, writing `a=ell q`, direct differentiation gives

```
a_corr=a+(ik/Gamma) integral zeta r chi-(A_c/Gamma)J m,
(a_corr)_dot=m+(ik/Gamma) integral zeta r b.               (5)
```

This is an actual finite-k identity, not just its zero-k version. Its
last term is a retained finite-k production current. It is the useful
repair of the previously uncontrolled translation forcing.

## 3. Why the other finite-wave couplings are compact

On the rearrangement support set let `w0=-d zeta/d psi>=0` and use
`q=sqrt(w0) f`. The multiplication by `sqrt(w0)` is bounded for each
fixed smooth finite core. Since `grad zeta=-w0 grad psi`, B and C
in these weighted variables contain `sqrt(w0) grad psi` multiplied
by a periodic inverse-curl or gradient-Green operator. They are compact
from scalar `L2` into the weighted core Hilbert space. The pressure
map `P_k` from that core space into scalar `L2` is also compact;
it has the full inverse-curl and pressure smoothing. `R_k` and `S_k`
are compact. Their operator norms and compactness are for each fixed
core radius; no uniform shrinking-core estimate is hidden here.

Let Pi be the actual orthogonal projection onto `ker A`, including
ambient orbit averages. Global inversion can exchange two different
closed ambient cells, so an odd function can be constant with opposite
values on those cells. Those zero modes are retained. For
`s=epsilon z`, `Re z>=c>0`, `|z|<=C`,

```
epsilon (s+A)^(-1) ->Pi/z strongly,
||epsilon (s+A)^(-1)||<=1/c.                              (6)
```

This is the spectral theorem for the skew-adjoint transport generator,
or its elementary mean-ergodic resolvent proof. It assumes no spectral
gap. Strong convergence and uniform boundedness, followed by a finite-rank
approximation of a compact operator, give operator-norm convergence
after compact composition. The adjoint statement controls composition
on the opposite side as well. In particular a compact C satisfies
`epsilon^2 R_A C R_A -> Pi C Pi/z^2` in operator norm. Treating
this limit as zero without its Pi factors would erase real ambient
responses. The bounded coordinate identity instead supplies the exact
observable fact `Pi v=0`, which is all that the leading mean needs.

## 4. Core reference and the stationary projection actually needed

The radial0036 weighted quadratic form has two translation zero modes
and a strictly positive complement in the odd sector. For a sufficiently
dilute smooth stationary core constructed by the actual Green IFT,
the corresponding compact weighted Green operator changes continuously
in norm after core rescaling. The exact two translation modes remain
in its kernel. The radial positive complement therefore persists for
the same fixed finite core, rather than being inferred from a sampled
spectrum. The common translation projection is the actual row ell in
(4), not an arbitrary orthogonal gauge.

Inside the rearrangement support all regular contours remain small
deformations of the radial core contours, with periods bounded above.
The odd sector has zero orbital average. Thus the core transport
inverse is bounded there. Together with the positive quadratic form
on the translation-quotiented space, this gives a bounded stationary
inverse of `L_c` on that complement. Equivalently, near s=0 its
resolvent has the decomposition

```
(s-L_c)^(-1)=T ell/s+R_f(s),                              (7)
```

where `R_f` is bounded in a fixed neighborhood of zero. This is a
core-only statement: ambient perturbations are still represented by
chi and b in (3), not smuggled into the weighted space. The positive
form also defines the bounded core evolution on the quotient; it does
not assert bounded evolution for the full core/ambient velocity space.

For completeness, the stationary inversion does not require equating
two different projection gauges. In weighted variables write
`H=I-sqrt(w0) G sqrt(w0)`, `L_c=-A H`, and
`t_i=sqrt(w0) partial_i psi`. The source radial operator is exactly
the normalized `I-B_m` in0036. The constant part of the rescaled
periodic Green kernel vanishes on the odd sector; the remaining smooth
lattice correction is `O(core_radius^2)` in Hilbert--Schmidt norm.
Norm continuity of the weights follows from the fixed smooth Green IFT.
These facts give the stated positive complement of H.

Within the core one has explicitly
`A^{-1}t_x=sqrt(w0)y`, `A^{-1}t_y=-sqrt(w0)x`, with zero orbital
mean chosen. Integration by parts against the compact excess core gives
`<A^{-1}t_x,t_y>=Gamma` and its antisymmetric companion. Therefore
`ell g=0` is precisely the solvability condition that
`-A^{-1}g` be orthogonal to `ker H`. Solve
`H f=-A^{-1}g` with the bounded positive-complement inverse, then add
the unique translation making `ell f=0`. This proves a bounded inverse
on the actual invariant complement without replacing ell by ordinary
orthogonality. A Neumann resolvent expansion about that inverse proves
(7). The bound on core periods is only inside the small monotone
rearrangement support, not at the ambient separatrix.

## 5. Full low-frequency Schur feedback

For Re s sufficiently large, Laplace-transform (3), eliminate chi,
and use (7). Write `R_A=(s+A)^(-1)` and `R_c=(s-L_c)^(-1)`.
The exact vertical operator is

```
W=s+A-k^2(R_k R_A+S_k)-k^2 P_k R_c(B R_A+C),              (8)
```

and the exact physical mean row is

```
[s^2 I+k^2 V W^(-1)P_k T] mhat=s V0,   V b=<v b>.         (9)
```

The sign in (8) follows from the `+ik B,+ik C` convention in (3).
Both conventions describe the same operator, but interchanging only
one of their signs would change the feedback.

Set `epsilon=|k|` and restrict z to a compact right-half-plane acoustic
sector as in (6). Multiplication of (8) by R_A gives `I+D_k`.
Its fast-core and stretching parts have the form

```
k^2 R_A (compact, uniformly norm-continuous in s) R_A,
```

or contain only one R_A. The latter tend to zero; the former converge
in norm to the retained ambient averaged block. The translation part
is different and is retained:

```
-(k^2/s) R_A P_k T (ell B R_A+ell C).                     (10)
```

Here `P_0 T=A v`, hence `R_A P_k T` is bounded on the acoustic
sector, using `A r=v` twice and `P_k-P_0=O(k^2)`.
The row `ell B R_A` is bounded by its explicit primitive in (4)-(5):
`A R_A=I-s R_A`. Thus (10) is `O(|k|)`, not an unexamined
zero-frequency inverse. Put `C_0=R_0+P_0 R_f(0)B`. Therefore

```
D_k -> D_0(z)=-Pi C_0 Pi/z^2,
V W^(-1)P_k T -> C_v=<v tensor v>                       (11)
```

uniformly on compact acoustic sectors with `Re z>=c` chosen so
`c^2>||Pi C_0 Pi||`. Indeed the internal averaged block then has a
bounded inverse by its Neumann estimate. Since `Pi v=0`,
`(I+D_0)^(-1)v=v`, and also `V Pi=0`. These are exact complete-fluid
decoupling identities; they do not set the internal averaged motion
itself to zero. The leading
mean resolvent is consequently the actual acoustic resolvent

```
epsilon mhat(epsilon z) -> z(z^2 I+C_v)^(-1)V0.           (12)
```

This is a controlled complete finite-k feedback statement, not the
static k² Taylor coefficient from0144. Every compact operator is the
one from the full periodic Euler/label system. The translation and
ambient currents supplying its apparent singular terms were computed
in (4)-(5). It also indicates exactly which contour geometry and core
stationary estimate an applicability proof consumes.

## 6. The time-domain continuation that (12) does not silently supply

For sufficiently large Re s, (9) is the Laplace transform of the actual
Euler initial-value problem. Equations (8)-(12) additionally give a
low-frequency operator resolvent on the stated acoustic sectors.
Analytic continuation to those sectors does not by itself justify
moving the whole inverse-Laplace contour: higher-frequency poles or
large-response residues may be crossed. Likewise compact-sector
convergence does not bound the high-imaginary-frequency part of that
integral. The fixed-time energy estimate `exp(C_core t)` is insufficient
at `t=T/|k|` and is not used to close this step.

The remaining constructive time-domain task is therefore more precise
than0146's arbitrary pressure forcing: control the complete high-frequency
residue/tail for the actual common-velocity preparation, or construct
a licensed well-prepared actual initial phase that removes that tail
while keeping the physical leading mean and its action. An ordered
fixed-k then shrinking-core construction must estimate these constants
in that order; their dependence on the core cannot be discarded.

No autonomous time-domain isotropic law, unexamined acoustic pole on
the continuous spectrum, or parent completion is asserted by (12).
This new low-frequency result is an actual operator construction;
its exact finite-time Euler observation is the still-active continuation.
