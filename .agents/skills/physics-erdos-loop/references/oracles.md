# Oracle selection and verifier sensitivity

Use the strongest practical oracle for each claim and state exactly what verdict it earns.

| Claim | Preferred oracle | Maximum verdict |
| --- | --- | --- |
| Exact identity, ansatz residual, algebra, series, exact limit | SymPy or direct symbolic algebra | `symbolic_verified` |
| Finite algebraic, order, combinatorial, or topological theorem | Lean with audited axioms | `formal_verified` |
| Analytically unresolved root, eigenvalue, BVP/IVP, integral, or optimization | SciPy root/integration/optimization/sparse-linear-algebra method plus refinement and an independent route | `numeric_evidence` |
| Time-dependent PDE or nonlinear dynamics | Appropriate spatial discretization plus SciPy time integration, convergence, conservation/stability, and a method cross-check | `simulation_evidence` |
| Figure or visualization | Rendering tool | artifact only |

Separate genuinely different evidence modalities so each receives the right
verdict, while preserving the strongest useful combined statement. Do not split
a claim into trivial fragments merely to simplify review. Use exploratory
numerics to generate hypotheses only when useful; before production numerical
evidence, freeze the analytic specification and work exact algebra, calculus,
scaling, bounds, asymptotics, and applicable theorems as far as the claim
permits. Computation answers the named remainder rather than defining the
question after outputs are visible.

## Audit the objective bridge first

Before selecting or auditing an oracle, record this bridge:

```text
computed predicate
-> mathematical proposition
-> upstream object/ensemble/representation/observable licenses
-> analytic-closure receipt and irreducible numerical remainder
-> maximum scientific verdict
-> exact effect on the parent campaign objective
```

Also name the nearby questions the predicate cannot decide. A perfectly
verified finite inertia, sampled spectrum, residual, fit, or energy difference
does not acquire a missing symmetry, conserved charge, admissible ensemble,
branch identity, or physical observable through numerical accuracy. With an
absent bridge, the result is exploratory or provenance evidence only.

Proof establishes the encoded mathematical statement. Measurement, ordinary
numerics, and simulation test applicability or consequences; record them as
separate evidence scopes rather than treating measurement as the definition of
proof. Use the phrase “computer-assisted proof” only for a rigorous enclosure
of the discretization, truncation, roundoff, and claimed continuum conclusion;
do not silently invent a new verification status outside the registry schema.

## Symbolic checks

- Derive expressions from canonical inputs.
- Substitute the candidate into the governing equation.
- Assert the unsimplified and simplified structure where sign or branch information matters.
- Inspect the evaluated form before assigning an exact verdict. An unevaluated `Integral`, `Sum`, derivative, root object, or unresolved condition is not a symbolic proof merely because it appears inside an equality check. Preserve the failed attempt and change identity, representation, assumptions, or oracle; use independent refined numerics only as a cross-check when the claim remains exact.
- Test domains, branches, dimensions, symmetries, special cases, and limits.
- Use a representative coefficient, sign, branch, or convention mutation when
  it can expose a copied or insensitive residual.

An identity that holds only because both sides contain the same copied literal is not verification.

## Numerical checks

- Link the analytic-closure receipt, state the strongest exact result already
  obtained, and name the one proposition left for computation. “No convenient
  closed form” is not enough; equally, no universal proof of analytic
  impossibility is required.
- Confirm the background is stationary in the claimed admissible space, the
  ensemble and observable are licensed, the representation covers the claimed
  perturbations, and the error enclosure can reach the decision boundary. If
  any of these fail, use `NUMERICALLY_UNRESOLVED` or a representation-scoped
  verdict rather than refuting the physical candidate or question.
- Use the numerical formulation natural to the claim: for example `scipy.integrate.solve_ivp` for an ODE or method-of-lines system, `solve_bvp` for a two-point BVP, `scipy.sparse.linalg` for large sparse spectra, and appropriate SciPy quadrature, root, or optimization routines for those claims.
- For PDEs, state the spatial method (finite difference, finite volume, finite element, spectral, or another justified discretization), boundary implementation, mesh, time integrator, stability restriction, and error norm. A generic integrator does not validate an unspecified PDE discretization.
- Record precision, solver, mesh/domain, timestep, tolerances, and stopping criteria.
- Require the library's success/status result and finite outputs before evaluating physics assertions.
- Refine resolution, timestep, domain, and tolerance independently.
- Compare two methods or an analytically soluble limit.
- Track conserved quantities and discretization error.
- Separate solver convergence from physical-model validity.
- Before FFT differentiation, establish commensurability of every active frequency with the sampled window or measure endpoint closure and use an appropriate nonperiodic derivative. Equality between an FFT coefficient and that same coefficient multiplied by `(i*omega)^n` is representation self-consistency, not independent evidence. For a claimed dominant line, report its fraction of the relevant derivative norm or power and require the preregistered dominance threshold.
- Predeclare tolerances against a dimensional or scale-relative error model. If an absolute near-zero bound fails, preserve the attempt and require refinement, conditioning, or roundoff evidence before replacing it with a justified scale-sensitive bound; keep exact analytic nulls separate from numerical regression.

Numeric evidence does not become exact verification because it has many digits.

Use `substrate_framework.numerics` for shared evidence capture when it fits. Its helpers standardize solver failure, invariant drift, collocation residuals, and empirical order; claim-specific code remains responsible for the governing equations, operators, data, norms, thresholds, and interpretation.

## Formal checks

- Prove the exact statement needed by the physics claim.
- Reject `sorry`, `admit`, `unsafe`, new axioms, proof by importing the target, and silent theorem weakening.
- Inspect the target theorem's axiom footprint.
- Distinguish definitions and decimal attestations from derivations.
- Audit the map between the formal statement and the physical interpretation.

A theorem can be perfectly proved and still encode too weak a proposition.

## Mutation and counterexample audit

Choose the smallest representative set of scientifically meaningful mutations
that covers the load-bearing failure modes of a custom symbolic/numeric verifier
or translation layer:

- coefficient or normalization change;
- sign flip;
- wrong convention;
- parameter outside the claimed regime;
- perturbed imported constant;
- fabricated fitted parameter;
- broken boundary or initial condition.

Each selected mutation must fail a relevant assertion. Do not mutate every
parameter when one mutation exercises the same failure channel. If headline
numbers move while the tally stays green, the verifier does not establish the
headline claim. A kernel-checked Lean theorem instead requires exact statement,
import, proof-escape, axiom-footprint, and physical-encoding audits; do not add a
mutation whose only purpose is ceremony.

## Independent rederivation

Use an independent route for a load-bearing normalization, convention
conversion, or parameter elimination when the primary oracle and existing
evidence do not already close that risk. Reuse a valid recorded rederivation at
the same boundary. Agreement between two calls to the same copied formula is
not independence.

## Comparator gate

Empirical data may test a frozen prediction. It must not:

- choose among concepts before structural criteria are applied;
- appear as a derivation input under a new name;
- determine a tolerance or free parameter after results are seen;
- become a pass condition for a supposedly first-principles derivation.

Open the comparator gate after equations, conventions, selection criteria, and structural tests are frozen. Report disagreement, then continue improving the candidate set; do not refit the framework narrative.
