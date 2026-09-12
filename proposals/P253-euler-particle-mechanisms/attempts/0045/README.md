# P253/0045: independent review of the 0042 retained Euler state and memory theorem

This README preregisters one fixed-boundary, non-author review of root-owned
attempt `0042`. The coherent target is the strongest exact local-time theorem
actually supported for a finite-rank resolved Euler velocity, its full
complement state, a transported material tag, the same-field pressure and
material observables, and the corresponding full-observable Dyson/Mori--
Zwanzig identity. The review will keep the retained-state and retained-memory
representations separate until their equivalence and common domain are proved.

The strongest possible result is an honest exact representation with the
unresolved initial state retained. It is not a finite-dimensional autonomous
particle equation, an existence theorem for a generic infinite-dimensional
orthogonal-dynamics semigroup, or a P2/P4, stability, carrier, particle, or
quantum result.

## Independence, ownership, and activation boundary

`particle-balance-review` did not author or implement `0042` and owns only
append-only review artifacts under `attempts/0045`. Root retains ownership of
`0042`, its importable module/tests, central proposal/registry state,
validation, memory, and commits. The reviewer will not edit `0042`, the
canonical module or tests, active `0043`, any central/generated file, memory,
or Git history.

Before this freeze the reviewer opened only the corrected `0042/README.md`,
SHA-256
`9c37f046edde581edd3e03e0c69fadaffb5d5ef91c7b1001c87bc1797c013a87`,
and inspected repository/Git metadata and content-blind hashes. No `0042`
construction, correction receipt, result, validation, command, stdout, stderr,
module implementation, or test body was deliberately opened.

A broad preactivation symbol inventory inadvertently exposed only these
module-level snippets: the module docstring says it does not construct an
orthogonal-dynamics semigroup or close material moments; the names
`LinearRetainedMemory` and its constructor/return sites; and the test module's
import path. No formula, function body, test predicate, or output was exposed.
The user's prompt and target README had already frozen the review criteria and
the absence-of-semigroup boundary. These snippets are disclosed as
non-authoritative metadata and will receive no evidentiary weight.

The required memory search exposed a high-level supervisor warning that a
formal Dyson identity is not an existence theorem for an orthogonal-dynamics
semigroup, plus summaries of the previously reviewed `0004` material theorem.
Memory is only a pointer; neither statement selects the verdict or substitutes
for the activated bodies.

No further `0042` target body may be opened until root:

1. registers `0045` centrally with this independent owner and exact scope;
2. records the repository-schema invocation, stdout, stderr, and exit under
   `attempts/0045`; and
3. makes `attempts/0045/activation-schema.exit` contain exactly `0`.

After activation, the reviewer will pin this README and the activated target
hashes, perform one substantive pass, and request at most one bounded
evidence-driven correction pass. The accepted authority base is release
`v0.183.0`; the observed preregistration Git head is
`ca929ce584b462c033fde7866e6231076445ecd5`. Shared worktree changes belong to
their existing owners and will be preserved. The physics-skill preflight
completed with `ok=7`, `warn=0`.

## Frozen target inventory

The exact target artifacts are:

- `0042/construction.md`, `result.yaml`, and `validation.md`;
- `0042/readme-correction.md`, solely to audit the corrected Dyson target and
  its provenance;
- `src/substrate_framework/euler_retained_memory.py` and
  `tests/test_euler_retained_memory.py`; and
- the `0042/focused-pytest` exit/stdout/stderr receipt family.

The two activation receipt families are process provenance, not scientific
evidence. The following hashes were calculated without displaying the target
bodies. Central activation must preserve them or explicitly identify and pin
the replacement frozen boundary before review:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0042/construction.md` | `6e9993547fea6617150e4a7f8f72cefa75dc443026b4215f1e173e2767d1d277` |
| `0042/result.yaml` | `d23ebb23b0dc755ecb83926580ab49a932d689bfafa6e78f8037e80e6e8d6424` |
| `0042/validation.md` | `e8efa13afae3422c4e4526394c9dfb35da08fed4709fe66233c97a612758de28` |
| `0042/readme-correction.md` | `3d48d6f4d23cee710c0a61de0646232e787a0f19c225d7f90a8dc60606b8f4f1` |
| `src/substrate_framework/euler_retained_memory.py` | `75ab59d079de41481f1bb85d6f9a97a7d230ab00d3458236e43456f4f9362f0d` |
| `tests/test_euler_retained_memory.py` | `46159071e1c7b90952643c2f24efb626ec3e2e26e63b435abea1f0dddf422e92` |
| `0042/focused-pytest.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0042/focused-pytest.stdout` | `4095c3a2c9393f1fc49c7273a9083197179e91f6bb23de146808ec4b1bdcf57c` |
| `0042/focused-pytest.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |

If activation names another implementation or test as part of the same fixed
transaction, its path and pre-review hash must be appended before opening.
Otherwise no API, test, or later body may be imported into this review by
proximity. Active `0043`, P2/P4 work, and every later edit outside the activated
boundary are excluded.

## Frozen suppliers and non-inheritance

Reviewed `0001/0004` may supply only its exact smooth local Euler material-tag
balances, same-field pressure terms, radial-swirl Poisson source/quadrupole,
harmonic-ball rates, and the initial two-disjoint-swirl/nonclosure examples at
their reviewed domains. Current supplier hashes are:

| Supplier artifact | SHA-256 |
|---|---|
| `0001/material-balances.md` | `9171a613dc7626fa59360d73cda709e17f445bba31422c3ceee12125898dc5dc` |
| `0004/review.md` | `d2b355b213e7e3a4c82f204dc49c2cf7d1f52cde27aa720cb81d04bf31f2d019` |
| `0004/verdict.yaml` | `f8b187a7f5783b3a7961f239668947b7071f368f56557e96a1211b24808af52c` |
| `src/substrate_framework/euler_material_balance.py` | `5ded2b458ccf46d7f32fe87618b96f414477c9069e4115aad342ffec7a84876b` |
| `tests/test_euler_material_balance.py` | `cb2ba5ba2c2a78ace4b3c75bd1891f03e807626b1f1c9500fd8a9cde118b0a51` |

Reviewed `0009/0011` may supply only the exact finite-dimensional linear split,
operator conjugation/Dyson identities, and executed periodic-shear memory
specialization at their stated scope. They do not establish a nonlinear Euler
Liouville or `QL` semigroup, forced-Euler well-posedness, or closure after
discarding the initial unresolved state. Current supplier hashes are:

| Supplier artifact | SHA-256 |
|---|---|
| `0009/retained-euler-memory.md` | `d5b36b813fa3fc73d5dc25a2c0bba8c64b686b2847de1cf27c3a52aa9ee194c6` |
| `0011/review.md` | `5c3a938abfb555c6d84a8ceb8e51bcd62a2896d873e25597a12bc833786940c5` |
| `0011/verdict.yaml` | `9aab64398d4af05175e78d954410edbd63e7cc908f63f5fbba95cfd52d7cb1c1` |
| `src/substrate_framework/euler_shear_memory.py` | `821778ccd55e3f01ba55f769c33b6de7e5f5df745369fd51fe8d7a53e0b4d17f` |
| `tests/test_euler_shear_memory.py` | `ff8183a92d69b367b6c426e37e2245f9917657949572fd772cefbcbf359ca698` |

Standard local well-posedness for smooth incompressible Euler, transport by a
given sufficiently regular divergence-free velocity, Leray/Hodge theory, and
Sobolev product/commutator estimates may be used only with the actual spaces,
regularity, time interval, and dependence on prescribed data stated. A generic
Mori--Zwanzig slogan or formal semigroup notation is not an imported theorem.

## Unit A: pressure convention and exact finite-rank retained-state system

1. Fix physical Euler as

       rho_m(partial_t u+u dot grad u)+grad p=0,
       div u=0,

   with the decaying whole-space normalization. Taking divergence must give

       -Delta p=rho_m partial_i partial_j(u_i u_j),
       p=rho_m(-Delta)^(-1)partial_i partial_j(u_i u_j).

   Audit the repeated-index contraction, derivative order, Fourier sign,
   distinction between physical and specific pressure, and any equivalent
   `partial_i u_j partial_j u_i` form. A sign copied consistently into prose
   and code is not an independent check.
2. State the exact Sobolev/divergence-free space. Require the finite-rank
   projection `Pi` to map to smooth solenoidal fields and to be bounded on all
   spaces used; define `Q=1-Pi`. If kinetic-energy orthogonality or cancellation
   is claimed, verify that `Pi` is orthogonal in the physical `L2` metric rather
   than merely idempotent.
3. Apply `Pi` and `Q` to the full Leray-projected Euler vector field evaluated
   at `u=v+w`. Check signs, pressure equivalence, the constraints
   `Pi v=v`, `Qw=w`, and preservation of divergence. Adding the two projected
   equations must reproduce literal Euler, including its nonlocal pressure;
   no stress, forcing, or closure term may be silently dropped.
4. Derive resolved/complement kinetic-energy exchange and every stated
   momentum/current row. Internal exchange must cancel in the total field,
   with density and boundary decay explicit. A finite-rank state split is an
   exact reparameterization only while `w` and its equation remain part of the
   state.

Representative exposing changes are reversing the pressure sign, contracting
`partial_i partial_j(u_j u_i)` against the wrong free index in a vector row,
dropping the Leray projection from one block, using `Pi+Q!=1`, evaluating
pressure on `v` instead of `v+w`, or asserting an orthogonal energy split for a
nonorthogonal projection.

## Unit B: forced-Euler response map and history substitution

5. For prescribed smooth finite-rank `v(t)`, identify the exact complement
   IVP defining `W[v,w_0]`: equation, initial condition, divergence constraint,
   `Q`-range invariance, pressure reconstruction from `v+w`, and any forcing
   created by the prescribed resolved history. State the Sobolev index and
   local time as functions of norms of both `v` and `w_0`.
6. Audit the local existence proof at the derivative level actually used by
   Euler. A quadratic map `H^s -> H^{s-1}` is not a locally Lipschitz ODE on
   `H^s`; the proof must use the standard transport/energy method, a valid
   forced-Euler theorem, or a precisely weaker solution space. Verify
   uniqueness and continuous/dependent history maps before substituting `W`
   into another equation.
7. Distinguish the complement equation for arbitrary prescribed `v` from the
   coupled exact Euler split. The reconstructed `u=v+W[v,w_0]` is literal
   Euler only when the retained history also satisfies the projected `Pi`
   equation. The response map earns an exact closed history equation for `v`
   with parameter `w_0`; it does not make an arbitrary prescribed `v` an Euler
   solution.
8. Solve the material-tag transport with the same `u=v+W`, and substitute the
   same-field pressure and full material history into every observable. Check
   that centering, moving centroid, pressure force/moment/torque/energy flux,
   ambient momentum, and tag normalization are retained exactly. A velocity-
   only history substitution is not equivalent to the advertised full state.

Representative counterchecks vary `w_0` at fixed prescribed `v`, choose a
history that violates the `Pi` equation, restore a pressure tail outside the
tag support, or lower the Sobolev index until the transport/product estimate
used by the construction fails.

## Unit C: load-bearing unresolved initial state

9. Identify precisely how `w_0` enters the response/noise term and prove that
   it cannot be removed from the claimed state map. Reuse the reviewed
   same-local-tag/different-pressure construction only at its exact initial-
   time scope; verify both full Euler states, their common resolved/tag data,
   and the differing pressure-derived observable rate.
10. Separate non-Markovian dependence from nonuniqueness. Given `(v,w_0)` and
    the required regularity, the local coupled flow should remain unique. The
    obstruction is that resolved/tag data or a resolved history without the
    needed unresolved initial datum does not select the pressure/material
    future.
11. State the strongest implication: the complement state, or an exactly
    equivalent noise/memory datum encoding it, is load-bearing. This does not
    show that every reduced model needs the entire infinite-dimensional field,
    rule out special invariant manifolds, or prove a generic impossibility
    theorem beyond the exhibited class.

An exposing mutation sets `w_0=0` in the noise term while keeping an initial
state whose pressure quadrupole or material acceleration depends on nonzero
unresolved content.

## Unit D: full-observable Dyson/Mori--Zwanzig identity

12. Freeze the Liouville convention for observables: if
    `(e^(tL)A)(U_0)=A(S_t U_0)`, then
    `d/dt e^(tL)A=e^(tL)LA`. Define the observable projection `P`, `Q=1-P`,
    its domain, and whether it is linear, bounded, idempotent, or conditional
    expectation. Do not confuse this `P/Q` with the velocity-space `Pi/Q`
    split.
13. Derive the right-Dyson identity directly from `L=PL+QL` and verify the
    full-observable formula

       d/dt e^(tL)A
        =e^(tL)PLA+e^(tQL)QLA
         +integral_0^t e^((t-s)L)PL e^(sQL)QLA ds.

    Check operator order, all signs, integration endpoints, and the `t=0`
    derivative. The noise term is unprojected on the full left side; applying
    `P` to a term remaining in `Ran Q` would erase it and change the claim.
14. Apply the identity to the complete material-observable vector `R`, not
    only to a scalar or linearized velocity. Trace pressure through `LR`, the
    orthogonal dynamics, and the memory kernel. The initial unresolved state
    must appear in the evaluation of `e^(tQL)QLR`.
15. Audit the exact meaning of every exponential. On an infinite-dimensional
    nonlinear Euler phase space, `L` and especially `QL` are unbounded and
    need not generate global strongly continuous semigroups on the proposed
    observable space. Credit the identity only on a common smooth local-flow
    core, as a rigorously justified finite-time Duhamel formula, or under
    explicit semigroup hypotheses. The declared lack of a generic `QL`
    semigroup is a correct boundary, not a missing assertion to manufacture.
16. Reproduce the finite-dimensional linear matrix specialization and the
    periodic-shear example only as exposing algebra. Agreement there checks
    ordering and signs but cannot establish nonlinear Euler domain or
    semigroup existence.

Representative mutations swap `PL` and `LP`, reverse the Dyson convolution,
put `P` on the full left side/noise, omit the initial noise term, evaluate the
kernel on `e^(sL)` instead of `e^(sQL)`, or silently extend a local flow to all
time.

## Unit E: importable module and tests as bounded algebraic evidence

17. Inspect every public symbol in `euler_retained_memory.py`: inputs, shapes,
    projection identities, sign/order conventions, dtype/domain assumptions,
    error handling, and docstring exclusions. Determine whether it implements
    only finite-matrix Dyson kernels/history terms or claims an Euler PDE
    solver it does not contain.
18. Map each test assertion to one mathematical predicate. Require the tests
    to derive rather than encode the full/Markov/noise/memory decomposition,
    test `t=0` and at least one sign/order mutation, and keep the recorded pass
    count distinct from scientific scope. The focused receipt is regression
    evidence only after its source hashes and predicates are inspected.
19. Verify package/import behavior without assigning scientific authority to
    module presence. A finite-matrix exact oracle can corroborate the Dyson
    algebra; it cannot prove forced-Euler Sobolev existence, material-tag
    transport, pressure reconstruction, or a generic infinite-dimensional
    `QL` semigroup.

No production numerical observable, soft eigenvalue, small force, or energy
splitting is registered, so `small-ratio-numerics` does not bind this review.

## Oracle, verdict, and exclusions

The strongest oracle is a direct pressure/Leray and finite-rank split
derivation, followed by a local forced-Euler Sobolev proof and exact
tag/pressure/history substitution, an exposing initial-state pair, and an
independent Dyson derivation on the declared observable domain. The importable
module and focused tests may corroborate only finite algebra and sign/order
conventions.

The retained-state route and retained-memory route receive separate component
verdicts before one joined representation verdict. `established` requires the
exact local Euler/tag reconstruction and a valid common-domain Dyson identity;
`refuted` requires a counterexample under the stated hypotheses; `blocked`
names the missing projection, Sobolev, pressure, tag, initial-state, or
Liouville-domain construction. A support gap in a generic `QL` semigroup does
not refute a correctly conditional local identity.

After activation, write only `0045/review.md` and `0045/verdicts.yaml`, pin all
inputs and evidence roles, preserve the strongest useful theorem, and request
at most one bounded evidence-driven correction pass. Report final artifacts,
verdict, and next dependency with both
`herdr agent prompt w4:p2 "..."` and
`herdr agent prompt w4:p1 "..."`.

P2/P4, active `0043`, carrier persistence or stability, finite-dimensional
autonomous particle dynamics, reciprocal particle force laws, electromagnetic
or quantum interpretation, electron/neutrino identification, all-time Euler
regularity, and parent-campaign completion are excluded.
