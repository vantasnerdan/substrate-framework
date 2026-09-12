# Validation and scope

The repository-interpreter oracle checks algebraic consistency for the
material-current continuity substitution, gauge-variation cancellation,
Gauss-law Coulomb energy and force signs, Maxwell speed, paired
Poynting/fluid-work terms, and paired divergence-constraint terms. The
displayed action and PDE calculation, rather than the placeholder
`x+(-x)=0` regression checks, supplies the actual Poynting identity and
constraint-propagation proof. The oracle does not select the new constants,
prove PDE well-posedness or all-time regularity, or construct a particle.

The local `H^s` statement uses the standard transport/symmetric-hyperbolic
energy estimate on the joint divergence/Gauss constraint kernel. The
finite-rank history map uses a projection reducing that kernel and retains
unresolved initial data. Independent review must inspect this domain, the
material-map Lorentz variation, the gauge boundary class, and the Noether
decay assumptions before the extension theorem is consumed. The six oracle
assertions are algebraic evidence only; they do not prove the PDE or Noether
domain statements.
