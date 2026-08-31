# Scientific problem deconstruction

Use this reference before a campaign turns a physical question into symbolic or
numerical predicates. The purpose is to prevent a correct computation on the
wrong object from terminating the scientific search.

## Build an obligation graph, not a milestone list

Separate four levels:

1. **Campaign objective:** the complete positive object requested by the user.
2. **Obligation node:** a mathematical or physical prerequisite whose success
   licenses named downstream work.
3. **Candidate route:** one concept, representation, or mechanism for satisfying
   that obligation.
4. **Attempt:** one reproducible execution of a route and oracle.

For every obligation record:

| Field | Meaning |
| --- | --- |
| `positive_intent` | What must exist or become established |
| `requires` | Upstream licenses consumed |
| `pass_licenses` | Exact downstream propositions enabled |
| `does_not_license` | Nearby interpretations that remain unavailable |
| `maximum_verdict` | Strongest result any oracle at this node may earn |
| `failure_scope` | Candidate class a negative result can actually touch |
| `unlocks` | Next dependency nodes activated by success |

A successful prerequisite activates the next unsatisfied node. A failed route
leaves its obligation active and triggers the continuation ladder. Neither is a
campaign terminal state.

## License chain

Freeze the chain below before production numerics:

```text
mathematical object
-> symmetry/topology/conservation structure
-> ensemble and exact variational functional
-> admissible representation and function space
-> analytic scale, bounds, limits, and asymptotic structure
-> observable
-> irreducible numerical remainder
-> numerical approximation
-> scientific verdict
```

### Object and symmetry

Name the constructed object exactly. Distinguish an arbitrary tangent, compact
group action, normal mode, stationary point, relative equilibrium, and nonlinear
periodic orbit. For a collective clock, fixed-charge language requires a
globally defined compact action with fixed period and generator normalization,
action invariance, vacuum fixation, conserved Noether charge, and a
nondegenerate Legendre map. A weighted, tapered, eigenvector-defined, or
data-dependent tangent without these properties is kinematic only.

### Ensemble and variation

Distinguish fixed charge, fixed frequency, unconstrained statics, and released
dynamics. Derive the transform instead of borrowing its formula. Differentiate
the complete functional—including field-dependent generators, projectors,
inertias, constraints, and terms that vanish on the background—before applying
the background restrictions. Static relaxation followed by fixed-charge
post-processing is not a fixed-charge solution.

### Representation and admissibility

Establish chart smoothness, constraint closure, gauge or bundle descent,
holonomy, phase rank, kinetic-metric sign and rank, boundary conditions, and the
asymptotic principal symbol. A coordinate string, omitted contraction, sign
error, singular phase chart, or undefined kinetic pencil is a representation
verdict. It cannot refute the physical object.

### Observable and localization

State what the observable measures and what could contaminate it. A support
restriction, taper, wall pin, mask, compact extension, or core cutoff used to
define the candidate cannot establish natural localization or open-space box
independence. A quadratic inertia is not a Noether charge; a frozen potential
cross-term is not a relaxed pair law; a negative mode refutes stability rather
than stationary-point existence.

### Analytic remainder and numerical representation

Before designing a production numerical check, freeze an analytic-closure
receipt. It records the exact equations, domains, complete admissible first and
second variations, symmetries, constraints and zero modes; dimensions and
non-dimensional groups; dominant balances and scaling; available identities,
bounds, coercivity or monotonicity results; limiting cases, perturbative and
asymptotic forms; the continuum operator and threshold where applicable; and
the strongest conclusion already fixed by algebra, calculus, or a theorem.
Name one residual proposition that computation must decide, why the current
analytic ladder does not decide it, the remaining design freedoms, and the
maximum verdict. “No convenient closed form” alone is not a receipt, while the
receipt is not required to prove that no unknown analytic method exists.

Only then choose a representation. Require a stationary background in the
claimed admissible class, a conditioning-aware state-error estimate,
branch-identity observables across refinement, and a crossed mesh-by-domain
design: vary mesh scale at fixed domain and domain size at fixed mesh scale,
with basis/order crossed as needed. A single path that co-varies mesh and box
cannot separate truncation from boundary error. Require the full
correct-ensemble Hessian, the constrained tangent space and zero-mode gauge,
and an error budget that reaches the verdict quantity. Derive a fit model and
its free parameters before reading values; a fit with no validation degrees of
freedom is not evidence. Earlier numerical sampling is hypothesis-generation
or debugging only and cannot be promoted later without a fresh frozen design
and execution.

## Exhaustion and PR boundary

The continuation ladder—method repair, representation change, and alternative
concept—is applied after every failed route. It is a minimum response, not proof
of exhaustion. Exhaustion requires the certificate defined in `AGENTS.md`,
including coverage of historical and failure-generated candidates and a
coverage argument for infinite classes.

Keep all obligation nodes, attempts, reusable utilities, and supported
subclaims on the campaign branch. Open one campaign PR only after positive
completion or certified exhaustion. A PR, commit, review, accepted subclaim, or
clean milestone does not change the active obligation graph.
