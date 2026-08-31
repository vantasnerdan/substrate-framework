# Attempt 0005: phase-locked covariant M5 bridge

## Declared action and the pre-evaluation kinetic repair

Let `A` be a normalized rank-two Lorentz generator on an oriented spacelike
clock plane, `A^2=-T`, and let `H` be a normalized tangent-traceless symmetric
tensor.  Both are auxiliary tensor fields subject to covariant normalization
and covariant-constancy constraints enforced by local multipliers.  They carry
no on-shell energy.  Define

\[
 K=\frac12(AH-HA),\qquad
 R(q)=e^{qA}.
\]

In an adapted clock-plane frame,

\[
 A=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
 H=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
 K=\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\]

Then `R(q) H R(q)^T=cos(2q)H+sin(2q)K`.  If `psi` is a
charge-one complex Lorentz scalar, the tangent tensor

\[
 {cal Q}(\psi)=\operatorname{Re}(\psi^2)H
                +\operatorname{Im}(\psi^2)K
\]

obeys `Q(exp(iq)psi)=R(q)Q(psi)R(q)^T` exactly.

The selected action starts from the reproduced conditional P239 density but
uses the explicitly recorded minimal kinetic replacement.  If `P_A` is the
field-metric orthogonal projection of the curvature two-form onto its
`A`-directed clock-plane component, replace `F` by `F_perp=F-P_A F` in that
positive square and add the positive canonical kinetic square of the
tangent-traceless M5 tensor `B`.  Thus the old clock curvature is not silently
retained alongside the canonical term.  The orthogonal curvature sector and
all static energies are unchanged.  This repair is necessary because the
purely additive choice would contribute
`64 int b^2 |grad b|^2` to the diagonal inertia.

The remaining action delta is

\[
 -\frac12|\partial\psi|^2
 -\frac14\operatorname{tr}|D B|^2
 -W(|\psi|)-6\|B-{\cal Q}(\psi)\|^2,
\quad W(f)=3f^2-4f^3+2f^4,
\]

with `||Z||^2=tr(Z^2)/2` on the clock tensor doublet.  The full diagonal
compact action is

\[
 \psi\mapsto e^{iq}\psi,\qquad B\mapsto R(q)BR(q)^T,
 \qquad q\in\mathbb R/(2\pi\mathbb Z).
\]

It is faithful because of the charge-one scalar, fixes the exterior
`psi=B=0` pointwise, and leaves every term invariant.  All `A,H` indices are
typed Lorentz indices and their constraint contractions are scalars.  The
unchanged uniaxial sector has `psi=B=0`; the altered clock-biaxial kinetic and
potential response is the explicitly declared completion.

## Exact reduced fixed-charge problem

On a diagonal relative equilibrium,

\[
 \psi=e^{i\omega t}f(r),\qquad
 B=b(r)[\cos(2\omega t)H+\sin(2\omega t)K],
\]

the lock is `6(b-f^2)^2`.  The exact static potential, inertia, charge, and
fixed-charge energy are

\[
 V(f,b)=3f^2-4f^3+2f^4+3b^2+4b^4+6(b-f^2)^2,
\]

\[
 I=\int(f^2+4b^2),\qquad Q=\omega I,
\]

\[
 E_Q[f,b]=\int\left[\frac12(|\nabla f|^2+|\nabla b|^2)+V(f,b)\right]
 +\frac{Q^2}{2I}.
\]

The coefficient four in `I` is the exact weight-two contribution of `B`, not
a fitted normalization.  The Legendre map is regular whenever the nonzero
profile has finite `I`.

The complete first variations give

\[
 -\Delta f+(6-\omega^2)f-12f^2+32f^3-24fb=0,
\]

\[
 -\Delta b+(18-4\omega^2)b+16b^3-12f^2=0.
\]

Nothing that vanishes on the background was deleted before variation.  With
variations `(eta,xi)`,

\[
 \delta I=2\int f\eta+8\int b\xi,
 \qquad
 \delta^2I=2\int\eta^2+8\int\xi^2,
\]

and the full Hessian is

\[
 \delta^2E_Q=\delta^2E_0-\frac{\omega^2}{2}\delta^2I
 +\frac{\omega^2}{I}(\delta I)^2.
\]

The potential Hessian blocks are

\[
 V_{ff}=6-24f+96f^2-24b,\quad V_{fb}=-24f,\quad
 V_{bb}=18+48b^2.
\]

The last Legendre term is positive and rank one.  The envelope identity is
`dE/dQ=Q/I=omega`.

## Exact exterior dynamics and existence bridge

At `f=b=0`, the scalar mass is `m_f^2=6`.  The M5 clock doublet has total
quadratic stiffness `9b^2` (the reproduced `3b^2` plus the lock's `6b^2`) and
canonical kinetic metric, hence `m_B^2=18`.  Its diagonal charge weight is two,
so the lightest charged continuum edge in generator frequency is

\[
 m_*^2=\min(6,18/4)=9/2.
\]

The old three boost channels remain massless neutral channels, so the complete
unqualified radiation threshold is still zero.  The two new scalar real
channels and two canonical M5 clock channels enlarge the positive kinetic
rank by four; the diagonal charged threshold above is action-specific and does
not replace that zero neutral threshold.

For the radial Hilbert space set

\[
 J(f,b)=\int\left[\frac12(|\nabla f|^2+|\nabla b|^2)+V(f,b)\right],
 \qquad K=I/2.
\]

The potential is a sum of nonnegative terms because
`W(f)=f^2[2(f-1)^2+1]`, and in particular controls the quadratic norm. Its
expanded nonlinearities have degree at most four: `f^3`, `f^4`, `b f^2`, and
`b^4`. Radial compact embedding for powers strictly between two and six
therefore supplies the compact interaction bridge directly; no endpoint
truncation is imported. Consequently
the abstract hylomorphic minimizer hypotheses `(H1)--(H3)` of the primary
source apply to this two-field `S1` representation.

The binding inequality has an exact plateau witness.  Take
`f=1/2`, `b=1/4` on a ball and smooth them to zero over a fixed-thickness shell.
In the large-radius limit,

\[
 V(1/2,1/4)=37/64,\qquad f^2+4b^2=1/2,
\]

so

\[
 \frac{J}{K}\longrightarrow\frac{2V}{f^2+4b^2}=\frac{37}{16}
 <\frac92=m_*^2.
\]

The theorem therefore supplies a nonempty open charge set of finite-energy
fixed-charge minimizers with `0<omega^2<9/2`.  Palais symmetric criticality for
the spatial rotation fixed set and the diagonal internal action promotes the
radial critical point to the declared full field equations (modulo translation
and global-phase orbits).

The M5 split is necessary, not optional.  If `b=0`, then

\[
 V(f,0)=3f^2-4f^3+8f^4,\qquad
 \frac{J(f,0)}{K(f,0)}\geq5
\]

for every nonzero finite profile, including its positive gradient energy.
Such a state cannot realize the theorem's strict sub-`9/2` binding level.
Every earned hylomorphic minimizer therefore has `b` nonzero and dynamically
lifts the M5 tangent degeneracy in its core.

The tails are unsupported and action-derived:

\[
 f\sim e^{-k_f r}/r,\quad k_f=\sqrt{6-\omega^2},
\]

while `b` has homogeneous exponent
`k_b=sqrt(18-4omega^2)` and a quadratic source from `f^2`, so its actual
leading exponent is `min(k_b,2k_f)` away from resonance.  Both are positive.
The exact fixed-charge scaling identity

\[
 T+3U-3E_{rot}=0
\]

contains competing `R`, `R^3`, and `R^{-3}` powers.  Together with existence
on `R^3`, finite `I`, and the natural tails, this establishes intrinsic rather
than box- or taper-selected localization.

## Scope of the positive result

This attempt establishes the augmented action's exact O1/O2 data and a
theorem-backed O3/O4 phase-locked, necessarily split open-space minimizer.  It
does **not** establish O5.  The neutral massless boost channels, the complete
orthogonal M5 constraint sector, symmetry zero modes, and nonlinear harmonic
radiation still require one frozen full constrained stability transaction.
