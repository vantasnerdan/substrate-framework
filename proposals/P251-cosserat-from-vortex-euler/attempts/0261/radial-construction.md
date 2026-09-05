# Constrained radial nondegeneracy on an admissible balanced profile

This file checks and repairs the 0255 Bessel-buffer argument at its stated
scope.  The conclusion is an existence statement in the inverse-designed
straight-profile class.  It is not a finite-radius Euler existence theorem.

## 1. An admissible balanced base with a sacrificial buffer

Choose `b>0` before the first positive zero of both `J_0(lambda rho)` and
`J_1(lambda rho)`, and fix

\[
 \Phi(\rho)=cJ_0(\lambda\rho),\qquad 0\leq\rho\leq b. \tag{1}
\]

Fix the retained 0258 annulus strictly inside `(0,b)`.  Only a still smaller
interval `I_buf compactly contained in (a_core,b)` is used for spectral
perturbations.  Thus `[0,a_core]` is immutable, while `I_buf` is a declared
sacrificial Bessel buffer contiguous with the center at the **base** profile.
This contiguity is load-bearing: it removes the `Y_0` component from a regular
radial zero mode.  Merely inserting a Bessel segment after a general
transition would not do so.

Choose one sufficiently large common support radius `a`, one fixed tiny
strictly decreasing tail, and one common logarithmic collar near `a`.  Outside
`[0,b]`, construct two smooth monotone profiles with these same endpoint
neighborhoods:

- a steep profile drops to the tiny tail in an arbitrarily short interval,
  making `int rho Phi_rho^2` arbitrarily large while its remaining `L^2`
  contribution stays bounded, hence `Q>0`;
- a broad profile stays at a fixed positive fraction of its value over a long
  annulus before joining the same tiny tail, making
  `lambda^2 int rho Phi^2` dominate the bounded transition cost, hence `Q<0`.

The profiles can be joined by a convex interpolation which remains strictly
decreasing and agrees exactly on `[0,b]` and on the outer collar.  Since `Q`
is quadratic on this family, the strict sign change supplies an interior
simple balanced member `Phi`,

\[
 Q(\Phi)=\int_0^a\rho\left(\frac12\Phi_\rho^2-
                        \lambda^2\Phi^2\right)d\rho=0. \tag{2}
\]

Inverse design of `g` is smooth for positive labels, gives `g=0` on (1), and
retains the prescribed logarithmic edge.  Multiplying (2)'s profile equation
by `Phi` shows that `g` is signed somewhere in the outer transition.  The
spectral bumps below are supported in `I_buf`, so this already constructed
signed region and the exact outer collar are untouched.

## 2. Fixed endpoint domain and simplicity of a radial zero

Put

\[
 y=-\Phi'>0,\quad M=\rho y^2,\quad
 V=\frac{y''}{y}+\frac{y'}{\rho y}-\frac1{\rho^2},
 \quad L_0=-\partial_\rho^2-\rho^{-1}\partial_\rho+V. \tag{3}
\]

Near the logarithmic edge, direct differentiation of the 0253 label law gives

\[
 -\lambda^2+g'(\Phi)=d^{-4}-6d^{-3}+O(d^{-2}),
 \qquad d=a-\rho.                                      \tag{4}
\]

Thus the radial realization from the closed quadratic form in
`L^2((0,a),rho d rho)` is semibounded, has the limit-point physical endpoint
at `a`, and has compact resolvent.  A zero eigenvalue is isolated.  Its regular
radial eigenspace is one-dimensional: the weighted Wronskian of two regular
solutions is constant, and regularity at `rho=0` makes that constant zero.

If zero is absent, this profile already has radial nondegeneracy.  Suppose it
is present, with eigenfunction `z`, and set `q=z/y`.  The exact identities are

\[
 (Mq')'=-Mq/\rho^2,                                    \tag{5}
\]

and, under `y_t=y exp(t h)` with
`h in C_c^infinity(I_buf)`,

\[
 \dot V=h''+(M'/M)h',                                  \tag{6}
\]

\[
 C(h):=\langle z,\dot Vz\rangle_{L^2(\rho d\rho)}
 =2\int M\left((q')^2-q^2/\rho^2\right)h\,d\rho.      \tag{7}
\]

The profile reconstructed from the fixed outer edge is

\[
 \Phi_h(\rho)=\int_\rho^a y e^h\,ds.                  \tag{8}
\]

It preserves the entire inner Bessel profile exactly iff

\[
 \mathcal A(h)=\int_{I_{\rm buf}}y(e^h-1)\,d\rho=0.   \tag{9}
\]

At `h=0`, the two equality constraints have linearizations

\[
 A(h)=\int yh\,d\rho,                                  \tag{10}
\]

\[
 B(h)=\int yh\,D(\rho)\,d\rho,
 \quad D(\rho)=\rho y-2\lambda^2
                    \int_0^\rho s\Phi(s)\,ds.         \tag{11}
\]

All endpoint terms in (7), (10), and (11) vanish because the perturbations
are supported strictly inside the buffer.

## 3. The determinant is genuinely nonzero

Because the base profile is Bessel continuously from the center through the
buffer, a regular radial zero has, on the buffer,

\[
 z=A_zJ_0(x),\quad y=c\lambda J_1(x),\quad
 \int_0^\rho s\Phi(s)ds=c\rho J_1(x)/\lambda,
 \qquad x=\lambda\rho.                                 \tag{12}
\]

After division by the positive `A` weight `y`, the `B` and `C` weights are

\[
 D=-cxJ_1(x),                                          \tag{13}
\]

\[
 E=\frac{2M((q')^2-q^2/\rho^2)}{y},\qquad
 q=\frac{A_zJ_0(x)}{c\lambda J_1(x)}.                 \tag{14}
\]

Their convergent Bessel expansions at a positive small buffer are

\[
 D=-\frac c2x^2+\frac c{16}x^4+O(x^6),                \tag{15}
\]

\[
 E=\frac{2A_z^2}{c}-\frac{A_z^2}{12c}x^2
       +\frac{A_z^2}{64c}x^4+O(x^6)
  =\frac{2A_z^2}{c}+\frac{A_z^2}{6c^2}D
       +\frac{A_z^2}{48c^3}D^2+O(D^3).                \tag{16}
\]

The nonzero quadratic coefficient proves that the analytic curve
`rho -> (1,D(rho),E(rho))` is not contained in an affine line on every
sufficiently small positive buffer.  Hence three radii can be chosen there
whose evaluation matrix has nonzero determinant; two can also be chosen with
distinct `D`.  Smooth nonnegative bumps of unit mass concentrated at those
radii preserve both determinants by continuity.  For the resulting
`h_0,h_1,h_2`,

\[
 \det(A(h_i),B(h_i),C(h_i))_{i=0}^2\ne0,
 \qquad
 \det(A(h_i),B(h_i))_{i=1}^2\ne0.                     \tag{17}
\]

This verifies the 0255 expansion and supplies the previously implicit
admissible-base requirement.  It would be false to assert (12) on a detached
Bessel segment, where a `Y_0` component is generally present.

## 4. Exact restoration of core and balance

Apply the ordinary finite-dimensional implicit-function theorem to the exact
map

\[
 (t_0,t_1,t_2)\longmapsto
 \left(\mathcal A(\textstyle\sum t_i h_i),
       Q(\Phi_{\sum t_i h_i})\right).                  \tag{18}
\]

Use `(t_1,t_2)` to solve the two equalities as smooth functions of `t_0`.
The derivative of the simple zero eigenvalue along this exact constraint
curve is a nonzero multiple of

\[
 \frac{\det(A,B,C)}{\det(A,B)}.                        \tag{19}
\]

Therefore one arbitrarily small nonzero `t_0` removes zero from the radial
spectrum while preserving the entire frozen core, the support radius, the
outer logarithmic collar, and `Q=0` exactly.  Positivity of `y exp(h)` is
automatic.  Compact support of the perturbations makes the inverse-designed
label law glue with all positive-label jets at both ends of the buffer.

Finally, rank two of `(A,B)` leaves a core-preserving direction `k` with
`A(k)=0` and `B(k) != 0`.  Solving (9) exactly along that direction gives a
separate local balance coordinate `beta` with

\[
 \left.\partial_\beta Q(\Phi_\beta)\right|_{\beta=0}
 \ne0.                                                 \tag{20}
\]

Radial nondegeneracy is open, so it persists for small `beta`.  Equation (20),
not reuse of the already-spent spectral direction, is the scalar border for
the finite-radius translation equation.

## Verdict and exact scope

`radial_route_verdict: established as stated`.

There exists an admissible inverse-designed straight profile which has the
fixed nonzero Bessel core and 0258 annulus, the fixed logarithmic edge and
support, exact translation balance, a transverse balance coordinate, and no
radial zero eigenvalue.  The only zero modes of the straight scalar
linearization are then the two translations: `m=1` factorization gives those,
`m>=2` is strictly positive, and the radial block excludes zero by the
construction above.

This result does not construct a finite-`R` solution.  Its next consumer must
still build the global range inverse at the characteristic flat edge and close
the nonlinear parameter-dependent problem.  The alternative segmented
Sturm/Wronskian route is not needed for radial nondegeneracy, but remains a
materially different continuation if a later review finds that the declared
buffer cannot be included in the profile class.
