# Attempt 0008: complete generator and phase-space repair

## Exact failure of the provisional completion

For

\[
 S=\begin{pmatrix}1+a&u&v\\u&t+p&q\\v&q&t-p\end{pmatrix},\qquad
 A=\begin{pmatrix}0&0&0\\0&0&-1\\0&1&0\end{pmatrix},
\]

the full tensor generator is

\[
 [A,S]=\begin{pmatrix}0&-v&u\\-v&-2q&2p\\u&2p&2q\end{pmatrix}
\]

and therefore

\[
 \frac12\operatorname{Tr}([A,S]^T[A,S])
 =u^2+v^2+4p^2+4q^2.
\]

Attempt 0006 retained only the last two terms.  At its axis-lock coefficient
`zeta=1/4`, `(u,v)` has weight one and mass square one.  The complete charged
edge was therefore one, not `19/4`, and the `39/16` split witness did not bind.
That exact mechanism revokes the provisional O2--O5 verdicts without altering
the valid compact group action or the earlier route evidence.

## Minimal action repair

Keep every field and interaction of attempt 0006, but set the positive
axis-lock coefficient to the exact unit value:

\[
 V_{axis}=\operatorname{Tr}(S^2)-(N^TSN)^2.
\]

For general positive `zeta`, the spatial potential Hessian in
`(a,t,p,u,v,q)` is

\[
 \operatorname{diag}(5,6+4\zeta,18+4\zeta,
 4\zeta,4\zeta,18+4\zeta).
\]

At `zeta=1`, with kinetic metric
`diag(1/2,1,1,1,1,1)`, the generalized mass squares are

\[
 (10,10,22,4,4,22).
\]

The charged sectors are the scalar of weight one and mass square six, the
axis--tangent shear doublet of weight one and mass square four, and the
tangent-traceless tensor doublet of weight two and mass square twenty-two.
The complete generator-frequency edge is consequently

\[
 m_*^2=\min(6,4,22/4)=4.
\]

## Repaired strict hylomorphy and core splitting

At the same exact plateau

\[
 |\psi|=\frac12,\qquad S=\operatorname{diag}(1,1/4,-1/4),
\]

the unit axis lock changes the potential density to `45/64`, while the full
inertia density remains `1/2`.  Thus

\[
 \frac{J}{K}=\frac{2V}{I}=\frac{45}{16}<4,
 \qquad 4-\frac{45}{16}=\frac{19}{16}>0.
\]

The core split remains necessary.  If `B=0`, the scalar contribution obeys

\[
 \frac{2[W(f)+6f^4]}{f^2}=16(f-1/4)^2+5\geq5.
\]

The unit axis lock contributes `2(u^2+v^2)` to the potential while the shear
contributes `u^2+v^2` to the inertia, hence its ratio is four.  All other
static terms are nonnegative, so the full `B=0` ratio is an
inertia-weighted combination bounded below by four.  It cannot attain the
strictly lower split witness.  The `B` equation at `B=0` also retains its
nonzero source `-12 Q(psi)`.

## Canonical phase-space theorem bridge

Use the canonical phase space

\[
 X=[H^1(\mathbb R^3;\operatorname{Sym}(3)\oplus\mathbb C)]
 \times[L^2(\mathbb R^3;\operatorname{Sym}(3)\oplus\mathbb C)]
\]

relative to `S=NN^T`.  The Hamiltonian includes the canonical field
velocities, gradients, and the complete nonnegative potential.  The conserved
Noether charge is

\[
 C=\int\left[\operatorname{Im}(\bar\psi\dot\psi)
 +\frac12\operatorname{Tr}(\dot S[A,S])\right]d^3x.
\]

These `E,C` have zero value and derivative at the exterior, are translation
invariant, have the splitting property, and the positive exterior Hessian plus
quartic growth makes `E` coercive.  A translation-vanishing sequence loses
the subcritical nonlinear terms, so the complete quadratic pencil gives
`Lambda_0=sqrt(4)=2`.  The plateau has optimized hylenic ratio
`sqrt(45/16)<2`.  Benci--Fortunato Theorem 18 therefore supplies a nonempty
translation-compact stable minimizing set.

Velocity variation at a fixed-charge minimizer gives

\[
 \dot\psi=i\omega\psi,\qquad \dot S=\omega[A,S],
 \qquad C=\omega I.
\]

The remaining variation is exactly the full fixed-charge functional.  The
canonical cubic semilinear wave system is energy-subcritical in three spatial
dimensions; its positive conserved Hamiltonian controls the phase-space norm.
Equivariance and uniqueness of this global flow make the minimizing initial
data a relative equilibrium and make the translation/phase minimizing set
Lyapunov orbitally stable.

Finally `0<omega^2<4`, and every active tail is strictly massive:

\[
 6-\omega^2>2,\qquad 4-\omega^2>0,\qquad
 22-4\omega^2>6.
\]

The two neutral tensor channels both have mass square ten.  No numerical
profile, spectral sign, or time integration is needed for this repaired
license chain.
