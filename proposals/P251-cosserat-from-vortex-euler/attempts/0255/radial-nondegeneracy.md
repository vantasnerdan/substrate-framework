# Constraint-preserving removal of a radial zero mode

This is a continuation of the balanced profile construction, not a genericity
assertion.  It gives an explicit finite-dimensional condition and a Bessel
buffer in which that condition is nonzero.  The final profile preserves the
fixed inner constant-curl core, the exact logarithmic edge, strict radial
monotonicity, and the 0256 balance.

Fix a balanced straight profile and write

\[
 y=-\Phi'>0,\qquad M=\rho y^2,
 \qquad
 V=\frac{y''}{y}+\frac{y'}{\rho y}-\frac1{\rho^2}.     \tag{1}
\]

The radial linearization is

\[
 L_0=-\partial_\rho^2-\rho^{-1}\partial_\rho+V.       \tag{2}
\]

Equation (1) is just the differentiated-profile identity: the `m=1`
operator has null function `y`, while deleting its angular `rho^-2` term
gives (2).

## 1. Exact zero-mode and potential variations

Suppose first that `L_0 z=0` has a physical radial zero mode.  It is
one-dimensional: the radial Wronskian of two regular solutions is a constant
divided by `rho`, and regularity at the center makes that constant zero.  Put
`z=yq`.  Reduction of order gives

\[
 (Mq')'=-\frac{M}{\rho^2}q.                            \tag{3}
\]

Perturb the monotone derivative, rather than `Phi` directly,

\[
 y_t=y e^{t h},                                       \tag{4}
\]

where `h` is smooth and compactly supported in a transition interval away
from both retained endpoint neighborhoods.  The induced potential satisfies

\[
 \dot V=h''+\frac{M'}M h'.                             \tag{5}
\]

With the radial measure `rho d rho`, the zero-eigenvalue quadratic derivative
is therefore

\[
 \begin{aligned}
 \mathcal C(h)
  &:=\int_0^a\rho z^2\dot V\,d\rho\\
  &=2\int_0^a M\left[(q')^2-\frac{q^2}{\rho^2}\right]
                  h\,d\rho .                          \tag{6}
 \end{aligned}
\]

The second line follows by integrating (5) twice and using (3); no endpoint
term occurs because `h` is compactly supported.  For a normalized isolated
zero eigenfunction, `dot mu_0=C(h)/int rho z^2`.  Thus (6), rather than an
appeal to generic spectra, is the functional that must survive the profile
constraints.

## 2. The two exact profile constraints

Reconstruct the perturbed profile at the same support radius by

\[
 \Phi_t(\rho)=\int_\rho^a y_t(s)\,ds.                 \tag{7}
\]

Because (4) is supported in the transition, the outer edge is unchanged.
For every point interior to the perturbation support, (7) differs from
`Phi` by one common constant.  Hence the entire retained Bessel core, not
only its central value, is unchanged exactly when

\[
 \mathcal A(y_t):=\int_0^a(y_t-y)\,d\rho=0,
 \qquad
 A(h):=D\mathcal A[h]=\int_0^a yh\,d\rho.             \tag{8}
\]

Use the radial normalization of the 0256 balance,

\[
 \widehat Q(\Phi)=\int_0^a\rho
        \left(\frac12y^2-\lambda^2\Phi^2\right)d\rho=0.
                                                                  \tag{9}
\]

Writing `K(s)=int_0^s rho Phi(rho)d rho` and interchanging the order in the
variation of (9) gives the second linear functional

\[
 B(h):=D\widehat Q[h]
 =\int_0^a yh\left[\rho y-2\lambda^2K(\rho)\right]d\rho.          \tag{10}
\]

Equations (8) and (10) display the full core-value and translation-balance
costs.  Neither is assumed independent of (6).

## 3. A checkable three-border determinant

Choose three smooth transition perturbations `h_0,h_1,h_2` and set

\[
 \Delta_{AB}=
 \det\begin{pmatrix}A(h_1)&A(h_2)\\B(h_1)&B(h_2)\end{pmatrix},
\]

\[
 \Delta_{ABC}=
 \det\begin{pmatrix}
 A(h_0)&A(h_1)&A(h_2)\\
 B(h_0)&B(h_1)&B(h_2)\\
 \mathcal C(h_0)&\mathcal C(h_1)&\mathcal C(h_2)
 \end{pmatrix}.                                      \tag{11}
\]

If both determinants are nonzero, apply the ordinary finite-dimensional
implicit-function theorem to the exact nonlinear constraints
`A(y_t)=0` and `Qhat(Phi_t)=0`, using the amplitudes of `h_1,h_2` as borders
for the amplitude of `h_0`.  Along the resulting constraint curve, the radial
eigenvalue derivative is exactly

\[
 \dot\mu_0
 =\frac1{\int_0^a\rho z^2d\rho}
   \frac{\Delta_{ABC}}{\Delta_{AB}}\ne0.              \tag{12}
\]

Thus any radial zero moves away from zero for one sufficiently small nonzero
amplitude.  If the starting operator has no radial zero, no perturbation is
needed.  Positivity of `y_t` is automatic from the exponential
parameterization.  The exact endpoint neighborhoods and support radius do not
move, and the nonlinear constraints preserve the whole fixed core and `Q=0`,
not merely their first variations.  The rank in (11) also leaves an ambient
profile direction transverse to `Q=0`, so the simple translation border
needed at finite `R` is retained after selecting the final balanced profile.

## 4. Why the determinant can be earned in the retained-core design

The attempt may choose a small fixed Bessel core and an additional Bessel
**buffer** immediately outside it.  The buffer belongs to the unperturbed
profile but is not part of the core required to remain common after profile
selection.  The three perturbations in (11) are supported in disjoint small
subintervals of this buffer, so the retained core is untouched.

On the buffer set `x=lambda rho`.  A physical radial zero restricts there to
`z=A_z J_0(x)` for some nonzero constant `A_z`, while

\[
 y=c\lambda J_1(x),\qquad
 K(\rho)=\frac{c\rho J_1(x)}\lambda .                 \tag{13}
\]

The ratios of the three weights in (8), (10), and (6), after dividing by the
positive `A`-weight `y`, are

\[
 D:=\rho y-2\lambda^2K=-\rho y=-cxJ_1(x),             \tag{14}
\]

\[
 E:=\frac{2M[(q')^2-q^2/\rho^2]}y,
 \qquad q=\frac{A_zJ_0(x)}{c\lambda J_1(x)}.           \tag{15}
\]

The Bessel series give

\[
 D=-\frac c2x^2+\frac c{16}x^4+O(x^6),               \tag{16}
\]

\[
 E=\frac{2A_z^2}{c}-\frac{A_z^2}{12c}x^2
       +\frac{A_z^2}{64c}x^4+O(x^6)
  =\frac{2A_z^2}{c}+\frac{A_z^2}{6c^2}D
       +\frac{A_z^2}{48c^3}D^2+O(D^3).               \tag{17}
\]

The nonzero quadratic coefficient proves that the curve `(1,D,E)` is not
contained in an affine line on a sufficiently small positive buffer.  Choose
three distinct sufficiently small buffer radii where their evaluation vectors
have nonzero determinant, and choose nonnegative unit-mass bumps concentrated
at those radii.  Continuity of the three integrals then gives both nonzero
determinants in (11) after relabeling two bumps as the constraint borders.
This is an explicit construction criterion, not an unverified independence
assumption.  It also avoids a monotonicity assumption on `V`, which would be
in tension with the signed transition required by the 0256 balance.

The Bessel buffer costs no new physical parameter.  It may be chosen outside
the fixed annulus used by 0258, and after selection the final inverse-designed
`g` still vanishes on the entire retained core, has the same logarithmic edge,
and contains whatever signed transition is required by `Q=0`.

## Consequence and remaining construction

`radial_route_verdict: established as stated for removal of any physical
radial zero within the expanded constraint-preserving profile family`.
The exact bordered construction supplies a balanced profile whose straight
linearization has only the two translation zero modes; vertical translation
is gauged and radial translation is handled by the simple profile border.

This is not yet the full finite-radius theorem.  The remaining Candidate-I
achievement is to realize the linear result as a bounded right inverse in the
coupled defining-function/interior/boundary-graph spaces and close the
nonlinear estimates.  That theorem must produce one finite large-`R` solution
and `C^4` convergence on the fixed 0258 annulus.  No route-level exhaustion or
parent completion follows from this construction.
