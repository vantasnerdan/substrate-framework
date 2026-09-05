# Full axial Euler response: an exact compensated mean and acoustic-window theorem

This is a new construction, not another review of0144. The base is the same
smooth periodic planar Euler field `u=(v(x,y),0)` used in0144, with normalized
cell average `<.>`, `div v=0`, `<v>=0`, and `v.grad v=-grad p`. Pressure here
is pressure divided by density. The axial Fourier convention is `exp(i k z)`.
The Bernoulli lift with nonconstant axial velocity, and generic nonaxial
wavevectors, are separate transfers; they are not silently included below.

## 1. Exact velocity and Hodge equations

Put `A=v.grad`, `D=(-Delta_h)^(-1)` on mean-zero scalars, and `d=grad_h D`.
Let `P` be the horizontal Leray projection and `P0=P-mean`. For a horizontal
vector field `b` set

```
Q_k b = 2(-Delta_h+k^2)^(-1) sum_ij (partial_i v_j)(partial_j b_i),
L b = P0[-A b-(Dv)b]                  (div b=0, <b>=0),
T c = -(Dv)c                         (c is a constant horizontal vector),
F f = P0[-A-(Dv)] d f,
V f = <v f>.
```

Every inverse on a nonconstant cell mode has its actual periodic Green
meaning. In particular `Q_k-Q_0=O(k^2)` as an operator from horizontal
`L2` to `L2`; the source has zero average, and its derivative can be
placed on the elliptic Green operator. `F`, `Q_k`, and `Q_k d` are bounded
on the indicated `L2` spaces, uniformly for small k. No point pressure
or discarded ambient source occurs here.

For the exact full Euler perturbation write

```
w_h = m+a+i k d w_z,   div a=0, <a>=0, <w_z>=0.
```

Then `div_h w_h+i k w_z=0`, since `div d=-1`. Direct projection of
`w_t+A w+(w.grad)u=-grad_k pi` gives the exact system

```
m_dot = -i k V w_z,
a_dot = L a+T m+i k F w_z,
(w_z)_dot = -A w_z-i k Q_k a+k^2 Q_k d w_z.                 (1)
```

The mean pressure is zero: the source of its Poisson equation has zero
average, also for a horizontal field of nonzero divergence. The horizontal
mean in (1) is the physical Euler velocity, not the mean material
displacement rate. `L T=0` and

```
Q_0 T c = -c.grad p = c.A v.                              (2)
```

These follow by differentiating the stationary Euler equation under an
actual uniform translation. They retain the Galilean Jordan chain at
k=0 instead of replacing it by a fixed-cage potential.

Define the physical integrated mean `X_dot=m`, `X(0)=0`, and set

```
a=T X+y,
Z=w_z+i k v.X,
pi_r=pi+X.grad p
    =(Q_k-Q_0)T X+Q_k y+i k Q_k d w_z.
```

Equations (1)-(2) now give, exactly,

```
y_dot=L y+i k F Z+k^2 F(v.X),
Z_dot+A Z=i k v.m-i k pi_r,
m_dot=-k^2 C_v X-i k <v Z>,       C_v=<v tensor v>.         (3)
```

Thus the cell-pressure remainder is an identified part of the same Euler
operator. Setting `pi_r=0` would be an additional approximation, not a
consequence of closed streamlines.

## 2. What a bounded streamline primitive proves exactly

Assume the actual cell admits bounded unwrapped coordinates `r` with
`A r=v` in the transport domain. For example, choose `r=x-c_a` on each
bounded invariant cell. Its jumps are on invariant separatrices, where
the normal transport flux vanishes, so they contribute no delta term to
`A r`. This construction needs a covering by such invariant cells, not
just the existence of one closed core contour. It does not require a
stationary first integral in a collar or a positive minimum orbit
frequency. On a chosen component one can equivalently use any bounded
transport primitive with this property.

Skew-adjointness of A on periodic `L2` implies

```
<v Z>=<r Z_dot>-i k <r tensor v>m+i k <r pi_r>.
```

Consequently the compensated physical mean

```
p_c=m+i k <r Z>
```

satisfies the exact equation

```
(p_c)_dot=-k^2 C_v X-k^2 C_rv m+k^2 <r pi_r>,
C_rv=<r tensor v>,    C_rv+C_rv^T=0.                       (4)
```

The correction is an explicit observable/current map. It does not rename
`m` as `p_c`: both appear in (4), and `X_dot=m` remains exact. A change
of the primitive by a transport-invariant bounded function changes this
normal-form coordinate, not the physical equation (3).

For the correlation in0144, put `C_r(t)=<r tensor exp(-t A)r>`.
Then `R_v=-C_r''` and

```
integral_0^t (t-s)R_v(s) ds=C_r(0)-C_r(t)+t C_r'(0).       (5)
```

The first two terms are uniformly bounded. `C_r'(0)=-C_rv` is
antisymmetric. Pairing actual backgrounds `v,-v` cancels the last term;
without that pairing it is retained as an odd current contribution.
This establishes the bounded correlation primitive even when periods
diverge at separatrices. Equation (5) alone is not an acoustic-time
Euler approximation; equation (4) displays the missing response row.

## 3. A full-operator acoustic-window theorem, without a transport gap

Here is a sufficient theorem with an independently identifiable operator
hypothesis, rather than a prescribed acoustic evolution. Suppose the
actual group `exp(t L)` on horizontal mean-zero solenoidal velocities
is uniformly bounded for `t>=0`. More generally one may use a norm Y
in which that group is bounded, `F:L2->Y` and `Q_k:Y->L2` are bounded,
and the displayed stationary translation and pressure rows are finite.
The scalar group `exp(-t A)` is unitary by incompressibility; no bound
on its inverse or separatrix period is assumed.

Prepare the actual common velocity data

```
m(0)=V0, a(0)=0, w_z(0)=0, X(0)=0.
```

For every fixed slow time T, let `epsilon=|k|`, `tau=epsilon t`,
and `x=epsilon X`. In the original time, (3) is a pair of bounded
groups for `(y,Z)`, forced and coupled at order epsilon to
`(x,m,y,Z)`. Indeed `w_z=Z-i sign(k)v.x` and

```
||pi_r|| <= C(||y||+epsilon ||Z||+epsilon |x|).
```

Equations (3), variation of constants for the two groups, and the
finite-dimensional equations for `(x,m)` therefore imply by Gronwall
on `0<=tau<=T`

```
sup (|x|+|m|+||y||+||Z||) <= C_T |V0|,                    (6)
```

uniformly for sufficiently small epsilon. The exponent here is `C T`,
not `C T/epsilon`: the rapid groups have been controlled before this
estimate. This is the step which a fixed-time Euler bound `exp(C t)`
cannot replace.

In slow time (4) becomes

```
x_tau=p_c-i k <r Z>,
(p_c)_tau=-C_v x+epsilon[-C_rv m+<r pi_r>].                (7)
```

For negative k only the explicitly imaginary current correction changes
sign; the force remainder in (7) still has epsilon. Equations (6)-(7)
and Duhamel for the finite-dimensional oscillator prove

```
sup_(0<=t<=T/|k|)
 |m(t)-cos(|k| sqrt(C_v)t)V0| <= C_T |k| |V0|.             (8)
```

There is an analogous `O(|k| |V0|)` error for `|k|X` relative to the
matrix-sine solution. The physical mass normalization of the initial
data remains rho, as in0144. The leading speed is calculated from the
same field, not fitted. For the actual sixfold cell, `C_v=c_b^2 I`.
No time-reversal pairing is needed for this leading theorem; its odd
row affects the controlled next-order error and is available if a
paired next-order statement is desired.

This proves a sufficient full Euler operator theorem. It does not yet
assert its horizontal-group hypothesis for the selected compensated
finite-core array. The next section reduces what that array actually
needs, avoiding an unnecessary whole-spectrum requirement.

## 4. The smaller response-specific achievement

Equations (4) and (7) give an a posteriori analytic criterion for the
actual common-velocity solution. Define

```
E_k(T)=|k| sup_(tau<=T)|<r Z>|
       +|k| integral_0^T (|<r Z>|+|C_rv m|+|<r pi_r>|) d tau.
```

Whenever this quantity is finite, the physical acoustic observation
error in (8), with its right side replaced by `C_T E_k(T)`, follows
directly from (7). In particular `E_k(T)->0` suffices. This only asks
for the three indicated finite response rows; it allows an unbounded
full fast semigroup, slow logarithmic response growth, continuous
spectrum down to zero, and unobserved resonances. It is not a request
for unrelated whole-spectrum stability.

The criterion is useful because its rows are exact source operators and
their scaling is fixed before estimating them. It cannot be passed by
using the desired cosine as input or by measuring only `R_v` in (5).

## 5. Executed representation change: the exact resolvent Schur row

For `Re s` sufficiently large all following resolvents exist. Define
`R_L=(s-L)^(-1)` and

```
W(s,k)=s+A-k^2 Q_k[d+R_L F].
```

Laplace transforming (1), with the common-velocity data above, and
using `R_L T=T/s`, gives exactly

```
[s^2 I+k^2 V W(s,k)^(-1)Q_k T] mhat(s)=s V0.              (9)
```

Thus the full pressure, horizontal response and ambient feedback occur
in a named operator, not an isolated-cell inverse appended to the
answer. This identity does not assert resolvent continuation through
the continuous spectrum or an acoustic pole. Those require bounds.

At k=0 the self-energy in (9) has the exact decomposition

```
V(s+A)^(-1)Q_0 T=C_v-s V(s+A)^(-1)v.
```

Since `(s+A)^(-1)v=r-s(s+A)^(-1)r`, the remainder is `O(|s|)`
in any right-half-plane sector `|s|<=C Re s`. This is a genuine
gap-free low-frequency result for the transport row. At nonzero k,
the additional operator `k^2 Q_k R_L F` must still be controlled;
replacing W by `s+A` would omit the precise response in (4).

## 6. Actual core/ambient energy attempt and the remaining construction

For the selected smooth compensated cell, the total stationary vorticity
is constant in the ambient region. The exact zero-k vorticity equation
there is passive transport. Split the initial perturbation by its
transported core and ambient supports. At k=0 the ambient part `b`
obeys `b_t=-A b`; it drives the core through the complete periodic
inverse curl. The core-only isovortical part has the Euler
energy-Casimir form

```
H_core=integral f^2/w - integral f (-Delta)^(-1)f,
w=-d zeta/d psi >=0
```

on its weighted support, with the translation nullspace retained.
The radial0036 form controls the corresponding radial reference on
its nontranslation subspace. Its exact statement does not include
ambient perturbations, and nonradial steady IFT alone does not give a
uniform group estimate. In particular the core equation has a genuine
forcing `B b`, so differentiating its quadratic form produces
`2<H_core f,B b>`, not zero. Even a positive core energy and unitary
ambient transport leave this term. This is the specific mechanism
which prevents importing0036 as a full-cell stability receipt.

The response-specific continuation is now concrete: estimate the
composite rows of `B(s+A)^(-1)` into the translation-quotiented core
resolvent, then their pressure pairing `<r pi_r>`, on the exact odd
sixfold sector excited by a horizontal mean. A finite number of
translation rows is kept explicitly. On centrally symmetric closed
contours an odd scalar has zero orbit average; near a nondegenerate
separatrix its primitive can grow with the logarithmic period. This
suggests a weighted limiting-absorption/cohomological estimate, not a
positive separatrix gap. Establishing the necessary operator bound for
the actual core/ambient couplings is the next uncompleted construction.

### 6.1 A constructive pressure-row primitive, not a whole-spectrum condition

There is a useful further exact reduction of this block. Write its
translation-quotiented core equation as

```
f_dot=L_c f+B b,       b_dot=-A b.
```

For the particular physical pressure functional `ell_c f`, suppose
the stationary adjoint equation `h L_c=ell_c` has a bounded solution.
An inverse of the gapped radial core operator in0036 supplies this
solution on that reference subspace. The corresponding nonradial
core estimate is a transfer to establish, not a new whole-fluid gap.
If the one ambient adjoint equation `d A=h B` has an `L2` solution,
then direct differentiation gives

```
ell_c f=d/dt(h f+d b).                                    (10)
```

Hence the term `k^2 ell_c f` in (4) is removed by the actual additional
current correction `-k^2(h f+d b)`. A resonantly forced core amplitude
may grow linearly while this correction stays `O(k)` on `t=O(1/k)`.
One does not need to show that such a resonance is absent to use (10).
For finite-k equations with forcing additions `F_c,F_b`, the exact
remainder after this correction is `-k^2(h F_c+d F_b)`, with both
terms retained. This gives a concrete quantity for the next estimate.

The individual ambient adjoint equation above has a constructive
gap-free solution under elementary contour hypotheses. Assume the
relevant closed ambient contours are centrally symmetric, their
separatrices are nondegenerate saddles, and `h B` is represented by an
odd `L2` scalar `g` that is bounded near those separatrices. On a
contour with period `T(a)`,
inversion is a half-period shift, so `integral_0^T g=0`. Solve
`A d=-g` by integrating g around the contour and choosing its orbital
mean to be zero. Near a separatrix, `|d|<=2 T(a)||g||_infinity`.
Area coarea in streamfunction/flow-time coordinates gives there

```
||d||_2^2 <=4||g||_infinity^2 integral T(a)^3 da < infinity,
```

since a nondegenerate separatrix has
`T(a)<=C(1+|log|a-a_s||)`. This proves the needed *one functional*
primitive without a bounded inverse of A on all `L2`. On the remaining
contours the periods are bounded, and the mean-zero circle Poincare
estimate gives `||d||_2<=T_max||g||_2`. Values assigned
on the separatrix itself have measure zero. It uses a precise contour
geometry, not an undeclared positive lower bound on orbital frequency.

For the core-to-ambient inverse-curl coupling, the adjoint source is
smooth near the separated ambient separatrices because its source is
supported in the core; the relevant odd representation is preserved by
the centrally symmetric background. An `L2` source in the ambient
collar next to the core suffices because periods there stay bounded.
Verifying the bounded core adjoint row and contour
partition for the chosen cell is the exact applicability step.

Similarly, if a retained core translation coordinate obeys
`q_dot=ell_0 B b` and `d_0 A=ell_0 B`, then `q+d_0 b` is conserved at
k=0. This keeps the ambient impulse in the translation coordinate;
it is not a harmonic-mean reset. Its finite-k production terms must
be included in the eventual acoustic system.

Thus this representation change has already removed the need for a
blanket bounded full-cell group: stationary adjoint rows, their actual
ambient primitives, and the smaller finite-k forcing rows are a
constructive alternative. The remaining step is a uniform estimate
for those full finite-k forcing terms on the acoustic interval. Merely
writing (10), or using its k=0 version at t=1/k, would not supply it.

Primary-source comparison: the radial inviscid-damping theorem of
[Bedrossian, Coti Zelati and Vicol](https://arxiv.org/abs/1711.03668)
concerns strictly monotone radial backgrounds and appropriate regular
data. It is a useful alternative method for the response row, but is
not imported as a theorem about this nonradial periodic compensated
array. No numerical instability, fitted frequency, or finite cutoff is
used as a verdict here.

The strongest established results are the exact normal form, sufficient
full-operator acoustic-window theorem, resolvent representation, and
the individual adjoint pressure-row construction. For the frozen
actual-array target, route A's verdict is `blocked` with the specific
finite-k response estimate missing; route B's verdict is `blocked`
with the corresponding composite resolvent bound missing. These are
route verdicts, not obligation-level blocking or exhaustion.
The parent acoustic objective remains active. Neither fact refutes
that objective or changes the accepted conditional C-CST results.

Whole-field SO(3) averaging is also not a shortcut to an autonomous
isotropic wave: a distribution of direction-dependent cosine responses
generally retains a distribution of frequencies. A generic-direction
same-field closure needs its own actual response/action construction.
