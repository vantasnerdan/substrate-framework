# P226 impact analysis

The proposal is an additive framework extension. Existing accepted claims,
registries, release manifests, generated documentation, and induced-gravity
ledgers are unchanged. The new proposed claims remain outside the accepted
registry pending maintainer adjudication.

The implementation follows four ownership seams identified by the bounded
code audit:

- `dimensional_sine_gordon.py` gains the physical Noether stress that belongs
  with its existing physical Lagrangian, Hamiltonian, and residual.
- `sine_gordon_weave.py` owns the tetrahedral collective metric, crossing
  entropy, and composition into a total upstream coupling.
- `local_horizon_gravity.py` owns the explicit Clausius-to-Einstein theorem
  interface; `linearized_einstein.py` owns only weak and radiative consequences
  of that established nonlinear equation.
- `sine_gordon_fiber_source.py` owns the conserved 3+1 embedding,
  `sine_gordon_breather_gravity.py` owns the static breather prediction, and
  `sine_gordon_breather_radiation.py` owns the full-retardation prediction
  while reusing `tt_angular.py` and `numerics.py`.

GitNexus semantic queries located the normalized stress owner, TT projector,
conditional radiation consumers, induced inverse-coupling ledger, and sampled
integration owner before editing. Direct source inspection verified every
selected symbol and exposed one inherited duplicate private vector validator;
the proposal adds no third copy.

Directly affected tests are the dimensional sine-Gordon and numerics suites
plus four new focused suites. Existing consumers continue to receive the same
objects and signatures. Package-level exports make the new APIs discoverable.

The main review risk is scientific, not migration risk: reviewers must decide
whether the explicit topological-doublet, statistically isotropic weave, unique
action identification, and local-equilibrium hypotheses earn the four proposed
claims. The code and evidence keep those hypotheses visible and do not turn them
into accepted authority automatically.
