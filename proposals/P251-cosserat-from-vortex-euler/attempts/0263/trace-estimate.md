# Full-operator one-sided trace estimate

This calculation retains the complete normal second derivative in 0261.  It
proves a uniform all-mode trace estimate for the exact radial collar and gives
the corresponding stopped-diffusion representation for every nearby smooth
collar.  The last section constructs the one-sided-to-interior bridge needed
for a nonradial tame inverse.

## 1. Exact generator and boundary orientation

Write the defining-function linearization as

\[
 P_Tv=\sum_{j=1}^2X_j^2v+X_0v+c_Tv,
 \qquad X_j=T\partial_j,                              \tag{1}
\]

\[
 X_0=(2-5T)\nabla T\mathbin\cdot\nabla
      -\frac{\epsilon T^2}{1+\epsilon x}\partial_x. \tag{2}
\]

The Stratonovich diffusion whose generator is `P_T-c_T` is

\[
 dZ_t=X_0(Z_t)\,dt+\sqrt2\sum_{j=1}^2X_j(Z_t)\circ dW_t^j. \tag{3}
\]

The equivalent Itô equation is

\[
 dZ_t=\left[(2-4T)\nabla T
 -\frac{\epsilon T^2}{1+\epsilon x}e_x\right](Z_t)dt
 +\sqrt2T(Z_t)\,dW_t.                                \tag{4}
\]

The `T grad T` difference between (2) and (4) is the Itô correction; omitting
it changes the normal drift at first order.  If `D_t=T(Z_t)`, then

\[
 d\langle D\rangle_t=2T^2|\nabla T|^2(Z_t)dt,         \tag{5}
\]

\[
 (P_T-c_T)T=T^2\Delta_\epsilon^*T
                   +(2-4T)|\nabla T|^2.              \tag{6}
\]

At a solution of the defining-function equation, (6) and the edge eikonal
give

\[
 (P_T-c_T)T=2a_\epsilon+O(T),\qquad
 a_\epsilon=(1+\epsilon x)^2.                        \tag{7}
\]

Thus, after shrinking to `0<=T<=delta_1`, there is a uniform `beta_0>0`
such that the drift of `D_t` is at least `beta_0`, while its quadratic
variation is `O(D_t^2)`.  The edge `T=0` is an entrance boundary: the strong
solution can start there and immediately enters `T>0`, whereas a solution
started at `T>0` cannot reach it.  One direct Lyapunov proof applies Itô's
formula to `1/T` stopped on `eta<T<delta_1`; its drift is negative quadratic
as `T` tends to zero, excluding explosion of `1/T`, and then lets `eta` tend
to zero.

Let `tau` be the first hit of `T=delta_1`.  The barrier
`(delta_1-T)/beta_0` and (7) give

\[
 \sup_{0\le T\le\delta_1}\mathbb E\tau
       \le \delta_1/\beta_0.                          \tag{8}
\]

Iteration of the Markov property gives a uniform exponential moment after
the collar is chosen so that the exponent times the right side of (8) is
less than one.  Consequently a bounded, sign-indefinite `c_T` is harmless.
For `P_Tv=F` and inner exit datum `h` the only possible bounded solution is

\[
 \boxed{
 v(z)=\mathbb E_z\!\left[e^{\int_0^\tau c_T(Z_s)ds}
                   h(Z_\tau)\right]
 -\mathbb E_z\!\left[\int_0^\tau
 e^{\int_0^t c_T(Z_s)ds}F(Z_t)dt\right].}            \tag{9}
\]

The signs follow by applying Itô's formula to
`exp(integral c_T) v(Z_t)`.  Formula (9) prescribes data only at the inner
exit `T=delta_1`; its value at `T=0` is the selected physical trace, not an
arbitrary outer Dirichlet condition.  A constant killing shift may be
inserted and removed, but is unnecessary because of (8).

## 2. Straight collar and exact Fourier killing

For the circular collar `T=d`, `rho=a-d`, the diffusion (4) becomes

\[
 dD_t=b(D_t)dt+\sqrt2D_t\,dW_t^1,
 \qquad
 d\Theta_t=\frac{\sqrt2D_t}{a-D_t}\,dW_t^2,          \tag{10}
\]

\[
 b(d)=2-4d-\frac{d^2}{a-d}.                          \tag{11}
\]

Choose `delta_1` so that `a-d>=a/2` and
`0<b_0<=b(d)<=b_1` there.  The full Fourier operator is

\[
 \mathscr L_m v=d^2v''+b(d)v'
 -\frac{(m^2-1)d^2}{(a-d)^2}v.                       \tag{12}
\]

For `|m|>=2` put

\[
 k=(m^2-1)^{1/2},\qquad M=k^{2/3},\qquad y=Md.       \tag{13}
\]

The following two barriers retain the `d^2v''` term exactly.

### Occupation/source barrier

Set

\[
 W_m(d)=\frac{M^{-1}}{1+Md}.                          \tag{14}
\]

A direct substitution, with no first-order reduction, gives

\[
 -\mathscr L_mW_m=
 \frac{b(d)}{(1+y)^2}
 -\frac{2y^2}{M(1+y)^3}
 +\frac{y^2}{(a-d)^2(1+y)}.                          \tag{15}
\]

For `0<=y<=1` the first positive term is at least `b_0/4`; for
`y>=1` the last is at least `1/(2a^2)`.  The negative term is uniformly
`O(M^{-1})`.  Hence there are constants `m_0,c_*>0`, depending only on the
fixed collar, for which

\[
 -\mathscr L_mW_m\ge c_*\qquad(|m|\ge m_0).          \tag{16}
\]

The finitely many remaining modes are absorbed into the low-mode constant.
By Dynkin's formula or the entrance-boundary maximum principle, the discounted
occupation resolvent `R_m`, defined by
`mathscr L_m R_m=-1`, `R_m(delta_1)=0`, obeys

\[
 0\le R_m(0)\le c_*^{-1}M^{-1}
       \le C|m|^{-2/3}.                               \tag{17}
\]

### Inner-data barrier

Choose a fixed sufficiently small `alpha>0` and put

\[
 p_m(d)=\frac{\alpha k^2d^2}{1+kd^2},\qquad
 B_m(d)=\exp\!\left(-\int_d^{\delta_1}p_m(s)ds\right). \tag{18}
\]

Since `B_m'=p_mB_m`,

\[
 \frac{\mathscr L_mB_m}{B_m}
 =d^2(p_m'+p_m^2)+b p_m-\frac{k^2d^2}{(a-d)^2}.       \tag{19}
\]

Relative to the last term, the three positive ratios are bounded by

\[
 2\alpha a^2\delta_1,qquad
 \alpha^2a^2,qquad
 \alpha b_1a^2.                                      \tag{20}
\]

First choose `alpha` and then `delta_1` so their sum is less than one.
Thus `mathscr L_mB_m<=0`.  Moreover

\[
 \int_0^{\delta_1}p_m(s)ds
 =\alpha\left[k\delta_1-\sqrt{k}
       \arctan(\sqrt{k}\delta_1)\right]
 \ge c k                                             \tag{21}
\]

for all sufficiently large `k`.  If `H_m` is the discounted exit multiplier,
`mathscr L_mH_m=0`, `H_m(delta_1)=1`, comparison gives

\[
 0\le H_m(0)\le B_m(0)\le Ce^{-c|m|}.                \tag{22}
\]

## 3. The actual trace theorem

Let `v_m` be the moderate solution of (12) with source `F_m` and inner
Dirichlet value `h_m`.  Equations (17) and (22) prove

\[
 \boxed{
 |v_m(0)|\le Ce^{-c|m|}|h_m|
       +C\langle m\rangle^{-2/3}
          \sup_{0<d<\delta_1}|F_m(d)|}               \tag{23}
\]

for `|m|>=2`.  The `m=1` block has zero killing and satisfies the same
bounded representation without smoothing; `m=0` is one finite block and is
controlled by the already established 0261 radial spectral gap.  Therefore,
for every real `s`,

\[
 \|\gamma_0v\|_{H^{s+2/3}(S^1)}
 \le C_s\left(
   \|h\|_{H^{-N}(S^1)}
  +\|F\|_{H^1_dH^s_\theta((0,\delta_1)\times S^1)}
  +\|\Pi_{\le1}v\|_{L^2}\right),                    \tag{24}
\]

where any fixed `N` is allowed in the exponentially smoothed inner-data
term and `Pi_{<=1}` denotes the three finite low modes.  This is the
full-operator replacement for the toy kernel in 0261.  Its `2/3` gain is the
same transition scale, but (15) and (19) show that the normal diffusion has
not been discarded.

The inner Neumann value is not independent data.  Once `h` and `F` are
given, (9) or the exact Volterra equation determines it.  Conversely, an
already existing solution may be estimated using both of its inner traces,
as in the frozen 0263 estimate (10), but existence must not prescribe both.

## 4. Moderate branch equals the physical branch

The exact integrating factor from 0261 gives the second homogeneous derivative

\[
 v_m'(d)=C\frac{e^{2/d}d^4}{a-d}+\hbox{moderate part}. \tag{25}
\]

Under the exact conjugacy `u=(Phi/d^2)v`, with `Phi=e^{-1/d}`, the first term
in (25) produces `u` growing like `e^{1/d}` times a power of `d`.  It is not
in the physical quadratic-form domain.  A bounded moderate `v`, on the other
hand, gives `u=e^{-1/d}d^{-2}v`, which is flat and belongs to that domain.
Since the scalar equation has only these two local branches, the physical and
moderate selections coincide mode by mode.  Summability follows from (23).

Thus (9) is not an auxiliary probabilistic boundary condition: it is the
trace of the physical form solution.

## 5. The constructed two-sided extension and collar shrink

Fix `0<delta_0<delta_1`.  A two-sided extension is useful, but it must be
constructed rather than imposed on the physical solution.  Flatten the edge
and smoothly extend the coefficients from `0<=T<delta_1` to
`-eta<T<delta_1`, keeping the normal drift positive and the step-two bracket
uniform.  Extend the source with a fixed bounded linear extension operator and
put any fixed smooth datum, for example zero, on `T=-eta`.  Stop the extended
diffusion when it reaches either `T=-eta` or `T=delta_1`.

The closed half-collar `T>=0` is invariant.  The normal diffusion coefficient
vanishes at `T=0`, the normal drift is strictly positive there, and pathwise
uniqueness holds for the smooth extended SDE.  A path starting at `T>=0`
therefore never samples the arbitrary negative-side coefficients or datum.
Its stopped expectation is exactly (9).  On the other hand, the extended
expectation is a distributional solution on an **open** neighborhood of
`T=0`.  The uniform drift brackets make it an interior hypoelliptic solution
there.  This constructs a smooth extension of the physical one-sided
solution; it does not assume a boundary trace or a Cauchy extension.

The interior drift-Schauder estimate on
`-eta/2<T<delta_0` now yields, uniformly for a small smooth nonradial
coefficient family,

\[
 \|v\|_{\mathcal X^{s+2}(0,\delta_0)}
 +\|\gamma_0v\|_{\mathcal B^{s+2}}
 \le C_s\left(
 \|P_Tv\|_{\mathcal Y^s(0,\delta_1)}
 +\|\gamma_{\delta_1}v\|_{C^0}
 +\|v\|_{C^0(0,\delta_1)}\right).                   \tag{26}
\]

Here `mathcal X^{s+2}` controls `X_iX_jv`, `X_0v`, and all intrinsic
commutations of total degree at most `s`; `mathcal Y^s` is the corresponding
source space, and `mathcal B` is the induced trace scale.  The constants are
uniform because the extension, inward drift, coefficient norms, and bracket
determinant are uniform.  In the exact radial case, (24) identifies the
sharper ordinary boundary gain `H^s -> H^{s+2/3}`.  In ordinary Euclidean
notation the step-two bracket conversion costs finitely many derivatives at
every order; the intrinsic grading is the no-loss tame scale needed for
perturbation.

Estimate (26) is deliberately on `delta_0<delta_1`.  Without a shrink,
ordinary boundary regularity at the inner exit requires high-order inner data.
After the shrink, the inner datum is remote from the estimated region and its
`C^0` norm is enough: hypoelliptic interior regularity supplies the higher
derivatives.  This is the correct version of the `k_0`-only term anticipated
in 0261.

In the global physical problem, uniform elliptic regularity on
`T>=delta_0/2` controls that low inner norm by the global source plus the
physical form norm.  Thus

\[
 \|v\|_{\mathcal X^{s+2}(0,\delta_0)}
 +\|\gamma_0v\|_{\mathcal B^{s+2}}
 \le C_s\left(\|P_Tv\|_{\mathcal Y^s(D)}
  +\left\|\frac{\phi}{T^2}v\right\|_{L^2(d\mu_\epsilon)}
                         \right).                    \tag{27}
\]

The second term is the actual physical low norm under the conjugacy, not an
unweighted surrogate.  The constructed extension proves (26)--(27) for the
nearby nonradial coefficient family.  Differentiating the extended equation
with respect to a coefficient parameter and iterating the same local estimate
gives smooth tame parameter dependence.  The primary sources license the
interior estimate and stochastic uniqueness; invariance and the construction
above are the missing one-sided bridge.
