# Attempt 0006: full-spatial canonical M5 clock

## Why attempt 0005 is not the parent object

Before forming a stability operator, vary the projected spatial M5 tensor
outside the reduced clock ansatz.  For

\[
 S=\operatorname{diag}(1,b,-b)
\]

in the exact projected Landau--de Gennes potential, the director-amplitude and
tangent-trace derivatives are

\[
 \frac{\partial V_{M5}}{\partial e_{11}}=8b^2,
 \qquad
 \frac{\partial V_{M5}}{\partial t}=-6b^2.
\]

They do not vanish for a split core.  Attempt 0005 therefore established a
useful reduced binding mechanism but not a full-field relative equilibrium.
Its provisional O3/O4 licenses are revoked in the campaign manifest.  This
attempt changes representation rather than testing the wrong Hessian.

## Full declared physical field space and action

A constrained auxiliary orthonormal Lorentz frame supplies a unit timelike
vector `U`, a unit spacelike clock axis `N`, and an oriented clock-plane basis.
The frame normalization, orthogonality, and covariant-constancy equations are
local Lorentz-tensor multiplier constraints.  Their directions are frame
gauge, carry no on-shell energy, and are not counted as propagating modes.

The physical M5 variable is now an arbitrary symmetric tensor `S` on the
three-dimensional positive spatial space orthogonal to `U`: all six components
vary.  The exterior is the aligned rank-one projector

\[
 S_\infty=NN^T.
\]

Let `T=I-NN^T` and

\[
 B=TST-\frac12\operatorname{Tr}(TST)T
\]

be the tangent-plane traceless clock tensor.  With the same normalized
`A,H,K` construction as attempt 0005, define

\[
 {cal Q}(\psi)=\operatorname{Re}(\psi^2)H+
                 \operatorname{Im}(\psi^2)K.
\]

The full positive-energy action has canonical Cartan spatial-tensor kinetic
metric and the potential

\[
 V(S,\psi)=V_{M5}(S)+\frac14
 [\operatorname{Tr}S^2-(N^TSN)^2]
 +6\|B-{\cal Q}(\psi)\|^2+W(|\psi|),
\]

\[
 V_{M5}(S)=-\frac12\operatorname{Tr}S^2-
 \operatorname{Tr}S^3+(\operatorname{Tr}S^2)^2+\frac12,
 \qquad
 W(f)=3f^2-4f^3+2f^4.
\]

The kinetic energy is

\[
 \frac14\operatorname{Tr}(\dot S^2+|\nabla S|^2)
 +\frac12(|\dot\psi|^2+|\nabla\psi|^2).
\]

This is the action-specific positive kinetic completion permitted by issue
#189.  It preserves the exact M5 static potential rather than importing the
M7 kinetic pencil.  The axis term vanishes on the aligned uniaxial sector and
prevents the exact director-shear failure with a coefficient fixed at `1/4`.

## Global positivity and vacuum

For real eigenvalues `lambda_i` of `S`, put
`r=(sum lambda_i^2)^(1/2)`.  The norm inequality
`sum lambda_i^3 <= r^3` gives

\[
 V_{M5}\geq r^4-r^3-\frac12r^2+\frac12
 =\frac12(r-1)^2(2r^2+2r+1)\geq0.
\]

Equality requires one positive eigenvalue equal to one and the other two zero:
`S` is a rank-one unit projector.  In an `N=e_1` frame,

\[
 \operatorname{Tr}S^2-(N^TSN)^2
 =2S_{12}^2+2S_{13}^2+S_{22}^2+2S_{23}^2+S_{33}^2.
\]

Thus the axis term selects `S=NN^T` from the M5 vacuum manifold.  The scalar
and phase-lock terms are also nonnegative, with
`W(f)=f^2[2(f-1)^2+1]`.  The unique physical exterior is therefore
`(S,psi)=(NN^T,0)`, fixed pointwise by the diagonal compact action.

## O1 and complete O2 data

The diagonal action is

\[
 \psi\mapsto e^{iq}\psi,\qquad S\mapsto R(q)S R(q)^T,
 \qquad q\in\mathbb R/(2\pi\mathbb Z),
\]

where `R` rotates the tangent plane about `N`.  It is faithful through `psi`,
the action is invariant, and Noether's theorem gives a conserved charge.  On a
relative equilibrium its exact inertia is

\[
 I=\int(|\psi|^2+4\|B\|^2),\qquad Q=\omega I,
 \qquad E_Q=E_0+\frac{Q^2}{2I}.
\]

All six spatial tensor channels are present in the exterior pencil.  In the
adapted coordinates

\[
 S=\begin{pmatrix}
 1+a&u&v\\u&t+p&q\\v&q&t-p
 \end{pmatrix},
\]

the potential Hessian is diagonal at the exterior with entries

\[
 (a,t,p,u,v,q):\quad(5,7,19,1,1,19).
\]

The canonical kinetic metric in these coordinates is
`diag(1/2,1,1,1,1,1)`, so the generalized mass squares are

\[
 (10,7,19,1,1,19).
\]

The two scalar real components have mass square six.  The only charged linear
sectors have weights one (`psi`) and two (`p,q`), hence the exact charged
continuum edge in generator frequency is

\[
 m_*^2=\min(6,19/4)=19/4.
\]

There is no physical negative-kinetic or inert spatial channel.  Auxiliary
frame directions are quotiented constraints rather than zero-norm physical
fields.

## Full-space fixed-charge existence

Use the complete physical Hilbert space

\[
 X=H^1(\mathbb R^3;\operatorname{Sym}(3)\oplus\mathbb C)
\]

relative to the exterior, modulo translations and the compact phase.  Define

\[
 J(S,\psi)=E_0,qquad K=I/2.
\]

The exact positivity above makes `J` coercive relative to the unique exterior.
Every nonlinear interaction has degree at most four, so its derivative is
energy-subcritical in three space dimensions. For the translation group,
energy and charge vanish with their first derivatives at the exterior
(`EC-0`), are invariant (`EC-1`), have the local polynomial splitting property
(`EC-2`), and the positive canonical Hamiltonian is coercive (`EC-3`). These
are the exact hypotheses of Benci--Fortunato arXiv:1212.3236, Theorem 18. The
remaining hylomorphy condition is discharged by the strict plateau comparison
below.

The strict binding witness is the large-ball plateau

\[
 |\psi|=1/2,\qquad S=\operatorname{diag}(1,1/4,-1/4),
\]

smoothed to the exterior over a fixed-thickness shell.  Its volume densities
are

\[
 V=39/64,\qquad I=1/2,
 \qquad J/K\longrightarrow39/16<19/4.
\]

Surface gradients are lower order in the ball radius. After optimizing the
standing frequency, this is precisely the theorem's strict energy/charge
hylomorphy inequality against the lightest charged vanishing sequence.
Vanishing is excluded
by the strict comparison with the charged continuum; dichotomy is excluded on
the theorem's nonempty open charge set by the energy-charge binding
inequality.  Tightness modulo translations and weak lower semicontinuity give
a full-space minimizer, not an ansatz critical point.  It solves every `S` and
`psi` Euler--Lagrange equation.

The clock necessarily lifts the M5 degeneracy.  If `B=0`, all M5 and axis
terms are nonnegative and

\[
 W(f)+6f^4=3f^2-4f^3+8f^4,
 \qquad \frac{2[W(f)+6f^4]}{f^2}
 =16(f-1/4)^2+5\geq5>19/4.
\]

No unsplit state reaches the binding level.  Moreover the full `B` equation at
`B=0` has the nonzero phase-lock source `-12 Q(psi)`, so a nontrivial critical
point cannot stay unsplit.

The complete fixed-charge variations retain

\[
 \delta E_Q=\delta E_0-\frac{\omega^2}{2}\delta I,
\]

\[
 \delta^2E_Q=\delta^2E_0-\frac{\omega^2}{2}\delta^2I
 +\frac{\omega^2}{I}(\delta I)^2,
\]

for arbitrary six-component `delta S` and complex `delta psi`.  The last term
is the positive Legendre rank-one form, and `dE/dQ=omega`.

## O4 localization and O5 stability/radiation

The minimizer lives directly on `R^3`, has finite `E` and `I`, and has
`0<omega^2<19/4`.  Its scalar and charged-tensor homogeneous tail masses obey

\[
 k_\psi^2=6-\omega^2>5/4,
 \qquad k_B^2=19-4\omega^2>0.
\]

Neutral tensor channels have positive exterior mass squares
`10,7,1,1`.  The cumulative Noether charge is integrable, so every fixed
fraction defines a finite median/core radius independent of any box.  The
fixed-charge Derrick identity is

\[
 T+3U-3E_{rot}=0,
\]

which excludes a free dilation modulus and dynamically selects the scale.

Because the solution is a **global constrained energy minimizer in the full
physical space**, its complete fixed-charge Hessian is nonnegative.  Its
kernel is quotiented by the minimizing set's translations and compact phase;
any further degeneracy remains inside that minimizing set rather than becoming
a negative mode.  The canonical semilinear equations are energy-subcritical,
and the positive conserved Hamiltonian controls the energy norm.  Compactness
of the minimizing set modulo translations/phase then gives Lyapunov orbital
stability by the standard energy-charge argument.

The background is static in the co-rotating frame.  Its only active laboratory
harmonics are the weight-one scalar and weight-two clock tensor, and the strict
frequency inequalities put both below their action-derived continua.  Neutral
channels have no periodic source after the **full** Euler equations are
satisfied.  Perturbations may carry away positive neutral radiation, but
energy-charge orbital stability forbids radiative growth of the clock or an
exponentially growing constrained mode.  This is an exact O5 statement for the
declared physical quotient; no soft numerical eigenvalue sign is used.
