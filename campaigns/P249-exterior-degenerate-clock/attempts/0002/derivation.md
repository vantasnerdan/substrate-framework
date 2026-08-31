# Attempt 0002 exact asymptotic pencil

This attempt derives the full quadratic exterior dynamics from the same action
that earned the compact clock license. No frequency or mass threshold is
imported from M7, M5.32, or a finite-box spectrum.

## Exterior and basis

At the normalized exterior

\[
 M_0=\operatorname{diag}(-4,1,0,0),\qquad
 X_0=\eta^{-1}M_0=\operatorname{diag}(4,1,0,0),
\]

use the ten symmetric covariant channels in order

\[
 (00,01,02,03,11,12,13,22,23,33).
\]

The timelike spectral projector is simple even though the two tangent spatial
eigenvalues coincide. Its exact first variation is unique because the selected
timelike eigenvalue remains separated from every spatial eigenvalue.

## Why the curvature term has no quadratic principal symbol

For a perturbation `M=M0+epsilon W`, every derivative is `O(epsilon)` and

\[
 F_{\mu\nu}=[\partial_\mu M,\partial_\nu M]_\eta=O(\epsilon^2).
\]

The curvature density is quadratic in `F`, hence `O(epsilon^4)`. It contributes
neither kinetic metric nor stiffness/principal symbol at second order. This is
an exact order statement, not a small numerical coefficient.

## Exact stiffness and kinetic matrices

The projected M5.17 potential with `beta=c=1` and unit timelike pin has exact
quadratic coefficient matrix

\[
 K=\operatorname{diag}(1,0,0,0,5/2,0,0,3/2,3,3/2).
\]

The three boost directions are exact Lorentz-orbit tangents and therefore have
zero potential curvature. The spatial block reproduces
`diag(5,3,3,0,0,6)/2` on `(11,22,33,12,13,23)`.

The projector current supplies

\[
 G=\kappa\operatorname{diag}
 (0,1/9,1/16,1/16,0,0,0,0,0,0),\qquad \kappa>0.
\]

Thus `rank K=5`, `rank G=3`. The positive physical kinetic subspace is exactly
the three timelike-mixing channels `(01,02,03)`. The common kernel is the two
director-shear channels `(12,13)`. Five channels are statically stiff but have
no quadratic kinetics.

## SO(2) weights and clock channel

Under the O1 generator in the 2-3 plane:

- `(02,03)` is a weight-one propagating doublet and `(01)` is neutral;
- `(12,13)` is a weight-one inert doublet;
- `(22-33,23)` is the weight-two tangent-split doublet;
- the remaining diagonal combinations are neutral.

The clock split has positive exact stiffness

\[
 K[d(22-33)+b(23)]=3(d^2+b^2),
\]

but lies in `ker G`. A small-amplitude clock is therefore not an ordinary
linear normal mode of this action.

## Principal symbol and permitted verdict

For spatial covector magnitude `k`, the exact generalized exterior pencil is

\[
 H(k)=K+k^2G,\qquad H(k)v=\omega^2Gv.
\]

On `im G`, `K=0`, so all three physical branches obey

\[
 \omega^2=k^2.
\]

The kinetic form is positive there and the essential/radiation threshold is
exactly zero. On `ker G`, the generalized frequency problem is undefined;
stiff directions are algebraic/static at quadratic order and the two common-
kernel directions are inert. Assigning them an infinite or positive mass gap
would be a representation error.

This establishes the action-specific asymptotic dynamics while exposing a
route obstruction: the licensed SO(2) clock must be a nonlinear relative
equilibrium, or the action needs a positive kinetic completion for its
weight-two split. It does not refute nonlinear clocks.
