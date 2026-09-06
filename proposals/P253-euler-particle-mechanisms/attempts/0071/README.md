# P253/0071: independent review of final 0068 Euler--gauge extension

This README preregisters one fixed-boundary, non-author and non-implementer
review of final root-owned attempt `0068` only. The review will adjudicate the
full joint Euler--Maxwell--material-tag constraint reduction, its gauge and
functional domains, the exact forced-PDE and algebra boundary, the static
Coulomb calculation, and the proposed small-coupling continuation. It will
preserve the strongest exact conditional extension while distinguishing it
from an actually constructed stationary charged carrier.

No result from active `0070`, no P4/P5 conclusion, and no electron or particle
interpretation can enter this transaction.

## Independence, ownership, and activation boundary

`particle-balance-review` authored or implemented neither `0068`, its action,
derivation, API, tests, verifier, corrections, nor receipts. It owns only
append-only review artifacts under `attempts/0071`. The reviewer will not edit
the target, any source/API/test, central or generated records, memory, or Git
history.

Before this freeze, no `0068` README, derivation, result, source audit,
validation, pre-review precision record, verifier, captured output, hash
manifest, API, or test body was opened. Inspection was limited to registry and
Git metadata, filenames, content-blind SHA-256 hashes, accepted-authority
metadata, and the user's frozen criteria.

The review opens only after the coordinator at `w4:p2`:

1. centrally registers `0071` with this exact scope and reviewer;
2. records the repository-schema command, stdout, stderr, and exit under
   `attempts/0071`; and
3. makes `attempts/0071/activation-schema.exit` contain exactly `0`.

After activation, the reviewer will recheck every frozen hash, open
`0068/source-audit.md` and `0068/sol-high-pre-review-precision.md` first, pin
each cited source version and exact imported proposition, and only then open
the remaining target bodies, API, tests, and receipts. One substantive pass
may request at most one evidence-driven bounded correction package, followed
by one correction-only check. Final writes are limited to
`0071/review.md` and `0071/verdicts.yaml`, followed by a report to `w4:p2`.

The accepted authority base is release `v0.183.0`.
`governance/releases/current.yaml` and `governance/claims.yaml` have SHA-256
values
`7b889a3c0186c0c22e698fd4ccc69e19066b8262018696cfc41e573ec2f0ffd8`
and
`78dbe694359739ee15f7db2b2cb1ccde7f48b5af34f6d9cd451282732dd08983`.
The observed repository head is
`e72f55aa4f5521586e5af0cbf71b30ebb3948450`. The shared worktree contains
accepted 0069 closure and active 0066/0070 work, including modifications to
`euler_gauge_emergence.py`; all are preserved and excluded. The physics-skill
preflight completed with `ok=7`, `warn=0`.

The target evidence is exact algebraic and focused API evidence. No production
numerics or small-ratio numerical prescription binds this review.

## Frozen final 0068 inventory

The complete current target package is frozen content-blind:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0068/README.md` | `f2e3507f11de1fac39d9665ef07a5139a37b629014b7fcd7e4810971996274b4` |
| `0068/derivation.md` | `26db2c0e6d119644ee7ed9308826e104bf91999ac55d21d1c824d64f9be996bc` |
| `0068/result.yaml` | `ec029045305fee0bfabf53520fad2f4820b7ec775367d630e64848c0bdfe0d5e` |
| `0068/source-audit.md` | `67b8858b9809d3c79b44c6001c59e7459b6d9796fe7379780b82f73935cd8bec` |
| `0068/validation.md` | `c2e573ee17815473d419cc7d28cab40ba2d4d1c621f28f33e408154de95f1c9a` |
| `0068/sol-high-pre-review-precision.md` | `4c826663db94812e565cef9f5dfe9ac0912989207f403fb3a72fb52b4254f523` |
| `0068/verify_gauge_extension.py` | `7ae6492fb77eb983b4a84d2fe05a493037732dabc22b3f80add60e18e24b04cc` |
| `0068/artifact-hashes.sha256` | `716ab1c1dfe94deeeba682e06f51402b827dc74f7e3957bb2818b561bdd8fe18` |
| `0068/activation-schema.command.txt` | `0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8` |
| `0068/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0068/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0068/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0068/exact-check.command.txt` | `0b41e03e87a2e55dd5dc0a06ec657bfb66cf111168ce38165c104e7144653695` |
| `0068/exact-check.stdout` | `0df45bf928031df8ecba6895fb3ce8b7b2f3443eb62dfc6a626303af9e8e6f9c` |
| `0068/exact-check.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0068/exact-check.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0068/focused-pytest.command.txt` | `ece58d2b886eaac150a5adbad86919e6b1690fc9deb265c3dfa054e600769c44` |
| `0068/focused-pytest.stdout` | `9c078973f1624f33eef2fcc94142e43aea74fa0db60028d01031a3f3fda4f835` |
| `0068/focused-pytest.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0068/focused-pytest.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

The importable 0068 API/test pair is frozen separately:

| Artifact | Preactivation SHA-256 |
|---|---|
| `src/substrate_framework/euler_gauge_current.py` | `cea0d100f56542c16868a675fb8dc8f2fb438741ab653542c39219c0c9237e0d` |
| `tests/test_euler_gauge_current.py` | `be984b2344189d61ee0e7029aa7f7a2d32f9ee1c18f914ee7bb19ddc8d61de33` |

`euler_gauge_emergence.py`, its tests, and all `0070` artifacts are outside
the frozen inventory and may not supply definitions, theorems, repairs, or
evidence to this review.

## Fixed claim and verdict boundary

The maximum positive result is an exact, correctly normalized constrained
Euler--Maxwell--tag action/reduction and a conditional local forced-PDE/static
interaction theorem on its stated domain. A static Coulomb Green calculation,
constraint kernel, or formal small-`g` ordering does not itself construct a
stationary charged Euler carrier.

The review will distinguish:

- exact identities and finite-dimensional/kernel algebra;
- local Sobolev well-posedness for prescribed admissible data;
- conditional static field response for a fixed charge profile;
- formal or conditional branch scaling;
- an actual stationary same-field carrier family.

No lower evidence layer inherits the verdict of a higher one.

## Check 1: full joint Euler--Maxwell--tag constraint kernel

1. Write the complete action or Hamiltonian, with Euler velocity/momentum,
   incompressibility multiplier, material tag, electromagnetic potentials or
   fields, charge/current coupling, density, light speed, permittivity and
   permeability conventions. Derive every Euler--Lagrange equation and sign
   from that one functional.
2. Identify all primary constraints and gauge degeneracies before reduction.
   Derive the joint constraint matrix/kernel, including incompressibility,
   Gauss law, longitudinal/transverse Hodge rows, tag transport, harmonic/zero
   modes, and any pressure or scalar-potential multipliers.
3. Prove the stated kernel and quotient are exact in the declared function
   spaces. Check rank, nullity, compatibility/cokernel conditions, residual
   gauge transformations, and whether the reduced form is symplectic or only
   presymplectic.
4. Verify that eliminating pressure, scalar potential, or longitudinal fields
   does not discard a global charge, harmonic, endpoint, or material-history
   row. Algebra at one nonzero Fourier wavevector cannot establish a global
   Hodge reduction.
5. Type any finite-rank projection `Pi` on the mixed state
   `(u,chi,E,B)` before calling it `L2`-orthogonal. These components have
   different physical dimensions, and the material tag has no automatic
   Maxwell/fluid energy weight. Orthogonality therefore requires an explicit
   nondimensionalization and analytic weighted product Hilbert inner product.
   Otherwise retain only the bounded finite-rank projection actually needed
   by the history theorem and drop `orthogonal`.
6. In either formulation, prove that `Pi` is bounded on the declared `H^s`
   mixed-state space, that `Q=I-Pi` is bounded there, and that both subspaces
   reduce the constraint/history operator `K` in the exact sense claimed.
   Smoothing or finite rank alone does not establish physical orthogonality;
   failure of orthogonality does not weaken an exact bounded reducing
   projection or the retained-history map.

**Maximum verdict:** the exact joint constrained kernel/reduced equations on
the fully declared domain.

## Check 2: gauge, spacetime, decay, and first-moment domains

7. Freeze the spacetime domain, time interval, gauge group, gauge condition,
   temporal endpoint variations, spatial boundary conditions, and allowable
   large versus decaying gauge transformations. Check all integration-by-parts
   surface terms and residual gauge modes.
8. For whole-space charge, distinguish potential, electric field, and energy
   spaces: a nonzero total charge has a `1/r` potential and `1/r^2` field.
   State which quantities are in homogeneous versus inhomogeneous Sobolev
   spaces, which have finite field energy, and which cannot be required to
   decay or be square integrable simultaneously.
9. Derive the charge/current continuity equation from material-tag transport
   and confirm it propagates Gauss' constraint. State the regularity and decay
   needed for total charge, first moments, momentum, integration by parts, and
   time differentiation.
10. Check periodic and whole-space zero modes separately. A net periodic charge
   requires compensation; a whole-space Coulomb sector requires its surface
   flux. Neither convention may be silently imported into the other.

**Maximum verdict:** a gauge-consistent domain theorem with every zero mode,
surface term, decay condition, and finite observable typed.

## Check 3: local forced-PDE well-posedness versus algebra

11. State the actual coupled PDE and whether electromagnetic fields are
   prescribed, solved jointly, or eliminated. Verify Euler forcing includes
   the full declared Lorentz/electric term and that incompressibility pressure
   restores the divergence constraint.
12. Match the claimed local theorem to a precise Sobolev index, constraint
    manifold, initial compatibility conditions, gauge, continuation norm, and
    dependence on the material tag. Check derivative loss introduced by Hodge
    or Coulomb inversion and whether the force map closes in the same class.
13. Separate a standard local forced-Euler/Maxwell existence import from what
    is proved in `0068`. The algebra oracle may verify signs, contractions,
    constraint propagation, and scaling; it cannot prove quasilinear PDE
    existence, uniqueness, continuous dependence, or a stationary branch.
14. Confirm that the current used by Maxwell is the same transported current
    entering Euler work/force and the joint Noether identities, rather than a
    prescribed duplicate with no backreaction.

**Maximum verdict:** the strongest local conditional Euler--Maxwell--tag IVP
statement actually supported, with the oracle retained as algebraic evidence.

## Check 4: static Coulomb cross energy versus stationary carrier

15. For a fixed smooth charge profile, derive Gauss' equation, Coulomb
    potential, field energy, and two-profile cross energy with exact sign,
    `4*pi`, coupling, permittivity, density, and double-counting factors.
16. State the support, separation, integrability, first-moment, and
    self-energy hypotheses. Distinguish an exact convolution from its far
    `Q_1 Q_2/(4*pi*epsilon*d)` asymptotic and control the remainder.
17. Check translation differentiation and action--reaction only where the
    charge profiles and field solve the same variational problem. A static
    field cross energy is not automatically a material carrier force.
18. An actual charged stationary carrier must also solve the Euler momentum
    balance with electromagnetic stress/force, pressure, tag stationarity,
    gauge constraints, and finite energy. Record this as a separate theorem,
    not a consequence of the Green kernel.

**Maximum verdict:** exact static Coulomb response/cross energy for fixed
profiles, plus stationary-carrier existence only if the complete coupled
equations are solved.

## Check 5: small-coupling branch and fixed tag normalization

19. Freeze the base `g=0` Euler carrier, its domain and stabilizer, the tag
    profile `chi`, and the fixed nonzero integral `C_chi`. Check the exact
    definition and units of

        Q_g = g C_chi.

20. Derive the small-`g` orders from the coupled equations: charge/current,
    potentials/fields, field energy, Lorentz force, Euler correction, pressure,
    and tag correction. Distinguish parity in `g` from an assumed analytic
    branch and state which remainders are norm-controlled.
21. Identify the nonlinear map, domain/codomain, differentiability, linearized
    inverse or Fredholm decomposition, and every solvability condition needed
    for an implicit-function or Lyapunov--Schmidt theorem. A formal
    `E=O(g)`, force `O(g^2)` ledger is not branch existence.
22. Keep `C_chi` fixed along the family and verify that spatial rescaling,
    carrier deformation, and tag transport do not change it. Opposite sign
    under `g -> -g` must preserve the same base family and all other physical
    parameters.

**Maximum verdict:** exact conditional scaling and charge normalization;
actual branch existence only with the full operator/cokernel closure.

## Check 6: translation cokernel and Noether distinction

23. Compute translation kernel/cokernel modes of the stationary linearization
    in the actual gauge and function spaces. Check adjoint pairings, centering
    conditions, field tails, and whether gauge or pressure modes mix with
    translations.
24. Distinguish three propositions: global translation invariance of the joint
    action, conservation of total Noether momentum for a solved evolution, and
    solvability of the stationary branch equation. The first two do not erase
    a Fredholm cokernel or prove its forcing projection vanishes.
25. If centers are fixed, state the slice and reaction/constraint used. If
    total momentum balance supplies an orthogonality identity, derive it for
    the actual coupled source rather than inferring it from action--reaction
    language.

**Maximum verdict:** exact symmetry/Noether identities and separately proved
stationary solvability conditions.

## Check 7: stabilizer preservation and charge-profile persistence

26. Define the Euler coadjoint stabilizer and its action on `chi`. Verify that
    admissible relabelings or orbit generators preserve the tag profile, its
    integral, support/regularity, and every finite row required by the carrier.
27. Material transport gives `partial_t chi+u dot grad chi=0`; stationarity
    additionally requires `u dot grad chi=0` in the chosen frame. Check which
    condition is proved and on what time interval.
28. Verify that charge density/current built from `chi` persist under the same
    solved flow and deformation. A fixed laboratory profile, an initially
    compatible tag, and a material invariant profile are distinct objects.
29. Determine whether the electromagnetic backreaction preserves the original
    carrier orbit/stabilizer or requires an enlarged phase space and new leaf.
    The branch may not reuse an unperturbed stabilizer without proof.

**Maximum verdict:** exact tag/current conservation and the strongest actual
stabilizer-compatible persistence statement.

## Check 8: parameter-normalization invariants

30. Audit physical dimensions and normalization of `g`, `chi`, `C_chi`,
    charge density/current, potentials, fields, Maxwell kinetic coefficient,
    `c`, permittivity/permeability, energy, action, and force.
31. Identify redundant reparameterizations. Rescaling `chi` and compensating
    `g`, or rescaling the gauge potential and its kinetic/coupling constants,
    must leave physical invariants unchanged. A numerical value or claim that
    changes under such a convention is not selected physics.
32. State which fixed normalization makes “small `g`” meaningful and which
    dimensionless or norm-relative quantity controls the perturbation. Verify
    that `Q_g=g C_chi` and the cross energy are invariant under allowed
    convention changes.
33. Separate a free extension parameter from a value selected by Euler
    topology, action, or another accepted invariant. No electron charge,
    fine-structure constant, or quantum normalization may be fitted or inferred.

**Maximum verdict:** a convention-invariant parameter ledger and exact
physical combinations, not parameter selection.

## Source and oracle audit

34. Open the source audit first after activation and pin every external or
    accepted input by version/path/hash. Match its hypotheses, gauge, domain,
    regularity, constraints, and conclusion to the proposition imported; a
    standard theorem name is not a hypothesis transfer.
35. Map every verifier/API/test predicate through
    `computed algebra -> mathematical proposition -> joint action/domain
    licenses -> maximum verdict`. Inspect copied literals, hard-coded kernels,
    sign/normalization sensitivity, zero modes, invalid domains, and whether
    the tests exercise the complete constraint rather than isolated formulas.
36. Reuse captured receipts only while code and inputs match this frozen
    boundary. No production numerical run is authorized. A focused algebraic
    pass cannot support PDE or branch existence.

## Verdict discipline, exclusions, and deliverables

The review will issue one route-scoped verdict for final `0068`: established,
refuted with an explicit mechanism, or blocked by the named missing
construction. It will separately record the status of constraint algebra,
functional domain, local IVP, static cross energy, small-`g` scaling, actual
branch existence, translation solvability, tag persistence, and parameter
normalization. A declared hypothesis is legitimate scope; an absent
load-bearing bridge is a support gap, not automatically a falsifier.

Excluded from this review are:

- attempt `0070`, its emergence claims, modified API/tests, and receipts;
- P4 or P5 completion, electron/particle identification, quantum statistics,
  magnetic moment, or a quantitative electron relation;
- stability, all-time persistence, or a global no-go not explicitly proved by
  final `0068`;
- target/source/API/test/central edits, commits, production numerics, and
  unrelated downstream replay.

After activation the reviewer will write only `0071/review.md` and
`0071/verdicts.yaml`, preserve the strongest exact joint extension and its
conditional boundary, name at most one bounded correction if concrete, and
report the artifact hashes, verdict, and next dependency to `w4:p2`.
