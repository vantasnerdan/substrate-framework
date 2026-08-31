# P249 theorem bridge

The existence and stability bridge uses Benci--Fortunato,
arXiv:1212.3236, Theorem 18, at its abstract field-theory scope.  The source
allows a finite-dimensional internal field space and translations as the
noncompact group.  P249 supplies that internal space as `Sym(3) + C`.

The theorem is applied on the canonical phase space

\[
 X=\bigl(H^1(\mathbb R^3;\operatorname{Sym}(3)\oplus\mathbb C)\bigr)
 \times\bigl(L^2(\mathbb R^3;\operatorname{Sym}(3)\oplus\mathbb C)\bigr),
\]

relative to `S=NN^T`, with velocities `(dot S,dot psi)`.  Its two conserved
functionals are the full positive Hamiltonian and the diagonal-clock Noether
charge

\[
 C=\int\left[\operatorname{Im}(\bar\psi\dot\psi)
 +\frac12\operatorname{Tr}(\dot S[A,S])\right]d^3x.
\]

For a relative equilibrium, variation in the velocity variables gives
`dot psi=i omega psi`, `dot S=omega[A,S]`, and hence `C=omega I` with

\[
 I=\int\left[|\psi|^2+\frac12
 \operatorname{Tr}([A,S]^T[A,S])\right]d^3x.
\]

This contains both the weight-one axis--tangent shear doublet and the
weight-two tangent-traceless doublet.  Equivariance and uniqueness of the
energy-subcritical canonical flow then turn the constrained critical initial
data into the stated relative equilibrium.

The hypotheses are discharged as follows:

- `EC-0`: energy and charge, together with their first derivatives, vanish at
  the aligned exterior after writing `S-NN^T` as the Hilbert variable.
- `EC-1`: the local action and charge are translation invariant.
- `EC-2`: the quadratic gradients/velocities, the bilinear Noether charge, and
  local subcritical terms through degree four have the splitting property in
  the three-dimensional energy phase space.
- `EC-3`: the exact global M5 lower bound, positive unit-strength axis lock,
  positive phase/scalar locks, canonical velocities and gradients, and the
  positive exterior Hessian make the Hamiltonian positive and coercive
  relative to the unique aligned exterior.
- Strict hylomorphy: translation-vanishing sequences lose the subcritical
  nonlinear terms, so their energy/charge threshold is the complete quadratic
  charged edge `sqrt(4)=2`.  The large-ball split plateau has optimized squared
  ratio `J/K=45/16`, hence `sqrt(45/16)<2` with squared margin `19/16`.

Theorem 18 therefore gives a translation-compact stable energy-at-fixed-charge
minimizing set.  It is applied to the complete six-component spatial tensor,
complex scalar, and their canonical velocities, not to attempt 0005's reduced
two-amplitude ansatz.  The auxiliary Lorentz frame is an explicit
multiplier-constrained gauge quotient and is outside the physical phase-space
tangent.  No theorem statement is imported for particle identity, two-clock
dynamics, gravity, or empirical lifetime.
