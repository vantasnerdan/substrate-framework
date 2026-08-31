# Attempt 0003 complete fixed-J calculus

This attempt records the full variational formula before restricting a
background and then uses two independent Derrick directions to decide the
commuting split. It also corrects a scope omission in attempt 0001: a collapse
path must remain on the simple timelike spectral branch.

## Exact inertia functional

For the fixed generator `A`, define

\[
 \zeta(M)=AM+MA^T,
 \qquad U_i(M)=[\zeta(M),\partial_iM]_\eta,
 \qquad V(M)=[A,P_t(M)].
\]

Let `H_2(M)` be the positive two-form metric induced by the spectral Cartan
metric and `q(X,Y)=-Tr(XY)/2` the projector target metric. The coefficient of
`omega^2` in the complete action is

\[
 C[M]=\int d^3x\left(
 2\sum_i\langle U_i,U_i\rangle_{H_2(M)}
 +\kappa q(V,V)\right).
\]

Every dependence is retained: `zeta` varies with `M`, `P_t` and `H_2` are
spectral functions of `M`, and `U_i` contains the full spatial derivative.
The action is parity-even in `omega`, so `B=0` exactly.

## Complete fixed-J variations

With

\[
 E_J[M]=E_0[M]+\frac{J^2}{4C[M]},
 \qquad \omega=\frac{J}{2C[M]},
\]

the first variation is

\[
 \delta E_J=\delta E_0-\omega^2\delta C.
\]

The complete second variation is

\[
 \delta^2E_J=\delta^2E_0-\omega^2\delta^2C
 +\frac{2\omega^2}{C}(\delta C)^2.
\]

The last positive rank-one term is mandatory. For two independent variations
replace `(delta C)^2` by `delta_1 C delta_2 C`. Freezing `zeta(M)`, `P_t(M)`,
or `H_2(M)` deletes load-bearing pieces of `delta C` and `delta^2 C`.

Along a stationary branch, the envelope theorem gives

\[
 \frac{dE}{dJ}=\frac{J}{2C}=\omega.
\]

These are exact identities for any admissible full-field representation.

## Spatial Derrick identity

Under `M_R(x)=M(x/R)`, decompose the static energy and clock inertia as

\[
 E_0(R)=E_4R^{-1}+E_2R+E_0R^3,
 \qquad C(R)=C_1R+C_3R^3.
\]

The complete fixed-J scale condition at `R=1` is

\[
 -E_4+E_2+3E_0-\omega^2(C_1+3C_3)=0.
\]

This condition is necessary, not sufficient, and it retains both curvature
and projector-current inertia classes.

## Branch-preserving commuting split

Attempt 0001 used a fixed time eigenvalue while taking the split amplitude to
infinity. To remain strictly inside the simple-timelike branch, take

\[
 M=\operatorname{diag}[-\tau(x),1,d(x),-d(x)],
 \qquad \tau(x)=4+2d(x),\qquad d\ge0.
\]

Then `tau>1`, `tau>d`, and `tau>-d` pointwise, including the whole collapse
path. The time pin only adds another positive quadratic term. Candidate K's
auxiliary-axis lock also adds a positive quadratic term and does not change the
following homogeneities.

For `d(x)=a f(x/R)` with positive shape integrals, the exact reduced energy has
the form

\[
 E_J(a,R)=U_2+U_4+E_{\rm rot}
 =P_2a^2R^3+P_4a^4R^3+\frac{J^2}{4I_4a^4R},
\]

where `P2,P4,I4>0`. The two logarithmic stationarity equations are

\[
 3(U_2+U_4)-E_{\rm rot}=0,
\]

and

\[
 2U_2+4U_4-4E_{\rm rot}=0.
\]

Eliminating `E_rot` gives

\[
 -10U_2-8U_4=0,
\]

which is impossible for a nonzero core. Therefore the commuting split has no
stationary point at all, not merely no global minimum. The axis lock cannot
repair this because it only increases `U2`.

The same branch-preserving family with `R=a^{-2}` sends
`U2=O(a^{-4})`, `U4=O(a^{-2})`, and `E_rot=O(a^{-2})` to zero. This confirms
the nonattained collapse while respecting the selected timelike spectral
premise.

The verdict is scoped to commuting spatial derivatives. A noncommuting core
adds a positive static curvature term with different scale/amplitude degree and
remains live. The auxiliary-amplitude family also changes the inertia from the
curvature class to a canonical volume class and remains live.
