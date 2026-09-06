# P253/0056: independent review of the 0051 Schwinger--Hopf construction

This README preregisters one fixed-boundary, non-author/non-implementer review
of root-owned attempt `0051`. The review seeks the strongest exact classical
Schwinger--Hopf and finite-gate statement supported by the actual physical
Euler supplier. Algebraic `C^2`, `CP^1`, and Pauli identities are valuable
positive results, but they do not inherit an invariant Euler doublet, a
selected action unit, detector law, exchange character, or propagation cone.

The six load-bearing units below receive separate verdicts before the joined
0051 verdict. A support gap is not a refutation: the review will preserve the
strongest exact conditional theorem and name the minimum missing construction
for every stronger physical claim.

## Independence, ownership, and activation boundary

`particle-balance-review` did not author or implement `0051`, its correction,
API, tests, or receipts. It owns only append-only review artifacts under
`attempts/0056`. Root retains ownership of `0051`, its scientific corrections,
source inventory, APIs/tests, central proposal and registry, validation,
memory, and commits. The reviewer will not edit `0051`, any API/test, any
central/generated file, memory, or Git history.

Before this freeze, no `0051` README, construction, result, source audit,
validation, correction receipt, verifier output, API, or test body was opened.
Only filenames, content-blind SHA-256 hashes, central ownership metadata,
accepted authority metadata, Git metadata, and the user's specified review
questions were inspected. Substantive memory and source searches are deferred
until activation so they cannot leak the target argument into the frozen
criteria.

No target body may be opened until root:

1. registers `0056` centrally with this independent owner and exact fixed
   scope;
2. records the repository-schema command, stdout, stderr, and exit under
   `attempts/0056`; and
3. makes `attempts/0056/activation-schema.exit` contain exactly `0`.

After activation, the reviewer will recheck every frozen hash, open the target
and source bodies, and perform one substantive pass. At most one
evidence-driven bounded correction pass may follow. The reviewer will then
write only `0056/review.md` and `0056/verdicts.yaml` and report the result to
`w4:p2`.

The accepted authority base is release `v0.183.0`. At preregistration,
`governance/releases/current.yaml` and `governance/claims.yaml` have SHA-256
values
`7b889a3c0186c0c22e698fd4ccc69e19066b8262018696cfc41e573ec2f0ffd8`
and
`78dbe694359739ee15f7db2b2cb1ccde7f48b5af34f6d9cd451282732dd08983`.
The observed repository head is
`efd222bfbe5b02ac2c501074e45cdb5eafd6903a`. The shared dirty worktree
contains other agents' active work and will be preserved. The physics-skill
preflight completed with `ok=7`, `warn=0`.

This is an exact theorem/source/API review, not a new numerical design. If the
activated package relies on a soft numerical leakage, drift, splitting, or
gate threshold near a discretization or roundoff floor, the small-ratio skill
must be loaded before that evidence is adjudicated. Otherwise no numerical
floor prescription binds.

## Frozen target, API, test, and receipt inventory

The complete current `0051` package is frozen content-blind:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0051/README.md` | `321f8429919a3b6c201693ff2220ef00114fcf44958dc910d0a00be6420c448b` |
| `0051/construction.md` | `8df69d4c88f675af7fae4446656856d1b5f4b5069c544bff64bc3bdf8024d25d` |
| `0051/result.yaml` | `a889c0fb325bdb08db640dcc35dd5244c8a54be8d2db8fd471db48ce145dbe07` |
| `0051/source-audit.md` | `e355dd7046e93595ff1f8062f391bd9f28380e11391ebe240c382408d11db50b` |
| `0051/validation.md` | `3cd41bd877efcdc881be28eca92b1784c10caa53beda8729a181a7fe1be0fae8` |
| `0051/physical-normalization-correction.md` | `046a12c1e05714f3330e187019f97dcec3f92933bc5d4c4f921498c59cf366da` |
| `0051/activation-schema.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0051/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0051/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0051/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0051/focused-pytest.command.txt` | `d8d75be024e6ecd2c6056ce7c1abcbd82e9fa2c7936b1c10172a25edde9a4f0c` |
| `0051/focused-pytest.stdout` | `5fc73e80996797494093e06cbd2b935387b453b382c21a146e148922c2c52ed6` |
| `0051/focused-pytest.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0051/focused-pytest.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0051/normalization-correction.command.txt` | `d8d75be024e6ecd2c6056ce7c1abcbd82e9fa2c7936b1c10172a25edde9a4f0c` |
| `0051/normalization-correction.stdout` | `d1ad42517f60f817268aaf147ece8c68568802a70360e2e8d65f16773f033c31` |
| `0051/normalization-correction.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0051/normalization-correction.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0051/normalization-correction-v2.command.txt` | `d8d75be024e6ecd2c6056ce7c1abcbd82e9fa2c7936b1c10172a25edde9a4f0c` |
| `0051/normalization-correction-v2.stdout` | `ef3f62b56e136ea2c6fe5868700df626f91d87ca9ec74e69f125119686e4eae7` |
| `0051/normalization-correction-v2.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0051/normalization-correction-v2.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

The importable algebraic supplier is frozen separately:

| Supplier artifact | Preactivation SHA-256 |
|---|---|
| `src/substrate_framework/euler_schwinger_hopf.py` | `3e3a9981e8be6ca7f21cc6e30087dca2924c624e40c4ac019f08d7dc5b6f9e4e` |
| `tests/test_euler_schwinger_hopf.py` | `9868f214c652d8629b27ea7fe524e6fdd344e9d7022cd31bb4137c3fac257ee6` |

Hashes and zero exits establish provenance only. After activation, the review
will inspect the exact test predicates, runtime tally, source dependencies,
and correction history. An algebraic API cannot supply an Euler invariant
subspace or physical normalization merely because its tests pass.

## Fixed theorem and non-inheritance boundary

This is one fixed-claim review, not a selection among new mechanisms. The
candidate route is the actual `0051` map from a physically normalized
two-mode classical Euler block to Schwinger--Hopf variables, its reduced
`CP^1` geometry, and its finite-gate estimates. No empirical comparator is
used.

The review will not import a premise merely because it is named in `0051`.
Every claimed Euler supplier will be traced to its accepted or independently
reviewed statement and checked for the exact domain, sign, density, units,
time interval, and invariant-subspace content used. In particular, a KKS
tangent pair, prepared finite-window response, formal two-mode projection, or
positive quadratic block does not automatically provide an invariant
doublet, a globally selectable action level, or physical analyzer controls.

## Check 1: physical energy--action normalization

1. Reconstruct the physical KKS form and quadratic Euler energy on both modes,
   including density, spatial/Fourier volume, real-versus-complex convention,
   coordinate units, and generator sign.
2. Derive each physical action from the energy and frequency rather than
   identifying the KKS scale alone with action. Verify dimensionally that
   `energy = frequency times action` in the declared convention.
3. Treat unequal frequencies explicitly. Separate conservation of individual
   uncoupled actions from proportionality of total energy to one total action,
   and from the stronger `U(2)` invariance needed to rotate an arbitrary
   doublet without changing its Hamiltonian.
4. Check whether rescaling the two modal coordinates to a common canonical
   form changes the physical energy, gate generators, total-action level, or
   Euler embedding. A coordinate normalization cannot select a universal
   numerical action value.

**Maximum verdict:** exact finite physical action and energy normalization for
the stated two-mode block, with unequal-frequency and unit conventions
explicit. It cannot select `hbar`, quantize the action, or prove block
invariance.

## Check 2: exact Schwinger--Hopf reduction

5. Starting from the target's declared symplectic convention on `C^2`, derive
   the total-action moment map and its Hamiltonian `U(1)` phase flow, including
   every factor of two and sign.
6. For each positive regular action level, verify that the level set is the
   appropriate `S^3`, that the phase action is free with the stated period,
   and that the quotient is `CP^1` rather than merely a coordinate chart.
7. Recompute the Hopf/Pauli coordinates and their constraint, Poisson brackets,
   stabilizers, and orientation. Derive the reduced KKS/Fubini--Study area from
   the same normalization; do not import a familiar `2 pi` or `4 pi` factor
   from a different convention.
8. Separate exact classical reduction from any integrality, Chern-number,
   spin, or Hilbert-space interpretation. If quantization is mentioned, state
   the extra prequantization condition and externally selected action scale.

**Maximum verdict:** exact classical `C^2 // U(1) = CP^1 congruent S^2`
reduction with its convention-dependent area. It is not a quantum two-state
system or action-selection theorem.

## Check 3: complex-linear and antilinear operator split

9. For the declared complex structure `J`, expand

       V_C=(V-J V J)/2,    V_A=(V+J V J)/2.

   Verify directly that `V_C J=J V_C` and `V_A J=-J V_A`, that
   `V=V_C+V_A`, and that the split is unique for real-linear `V`.
10. Check all matrix and conjugation formulas against that convention. A map
    commuting with `J` must land entirely in `V_C`; a map anticommuting with
    `J` must land entirely in `V_A`. Audit any sign change induced by using
    `J^2=-1` or by switching temporal-frequency branches.
11. Determine which part is Hamiltonian/symplectic on the physical mode block
    and which part creates leakage, pairing, or Bogoliubov-type mixing. The
    algebraic split alone does not give a closed two-level evolution.

**Maximum verdict:** exact convention-correct real-linear decomposition and
its actual operator meaning on the declared block.

## Check 4: two noncollinear physical Pauli axes

12. Identify every proposed Pauli-axis generator as an actual physical Euler
    perturbation or control, not just a basis matrix on abstract `C^2`.
13. Verify that at least two noncollinear axes survive the same physical
    action level, constraints, stabilizer quotient, and domain, and that their
    flows are independently tunable on the claimed interval.
14. Check detuning and unequal modal frequencies: free phase evolution plus
    one fixed coupling does not supply arbitrary `SU(2)` control unless the
    two generated axes are genuinely noncommuting and the control protocol is
    physically realizable.
15. Distinguish a formal Pauli-coordinate representation from controllability,
    an invariant analyzer, or a measurement interaction.

**Maximum verdict:** two actual noncollinear physical control axes and the
finite-dimensional Lie algebra they generate. One axis, a coordinate choice,
or an unmaterialized coupling earns only a conditional control theorem.

## Check 5: gate-time leakage and action drift

16. Pin the projected block, complement, injection, norm, residual, initial
    data, and exact time interval. Separate leakage from phase error, detuning,
    action drift, and projection error.
17. Derive the gate time from the physical generator with units intact. For a
    coupling-controlled gate, track inverse-coupling growth of the gate time;
    a small instantaneous residual can accumulate to order one on that scale.
18. Audit the Duhamel or energy estimate including every constant and
    quantifier. Check whether the full linear Euler evolution is unitary,
    merely energy bounded, or nonnormal in the chosen norm, and include the
    resulting propagator factor.
19. Verify that action drift is computed from the physical total action and
    not inferred from norm preservation in differently scaled coordinates.
    State whether the estimate is absolute or relative and how unequal
    frequencies enter.
20. Match every API/test predicate to these continuum statements. Exact small
    matrices may validate signs and scaling while leaving the Euler residual,
    complement, and gate interval conditional.

**Maximum verdict:** a finite-time gate approximation with explicit leakage,
phase, and physical-action drift bounds. It cannot imply an all-time invariant
doublet or nonlinear Euler control.

## Check 6: actual Euler supplier and downstream boundary

21. Inventory every accepted claim, reviewed proposal result, and source
    theorem used to realize the two modes. Check same-field status, domain,
    physical normalization, pressure/exterior tails, accessible constraints,
    and the exact local or finite-window quantifiers.
22. Determine whether the supplier proves an invariant two-complex-dimensional
    Euler subspace or only a tangent block, prepared response, approximate
    projection, or conditional reduction. Preserve the strongest actual
    supplier statement and do not infer invariance from a matrix API.
23. Check whether the total action is conserved by the full supplied Euler
    evolution and whether any mechanism selects one action level. Classical
    conservation and reduction do not select an absolute quantum.
24. Keep observer and detector claims separate. No Born rule, projective
    reset/repreparation, two-copy exchange character, or measurement theorem
    follows from the Schwinger--Hopf geometry alone.
25. Keep propagation claims separate. Bare Euler has elliptic pressure and
    Galilean rather than Lorentz kinematics; a finite gate interval does not
    create a Lorentz cone or finite signal speed.

**Maximum verdict:** the strongest exact or conditional classical
Schwinger--Hopf theorem licensed by the actual Euler supplier. Missing
invariance, action selection, detector dynamics, exchange, or Lorentz
structure remains a named frontier rather than a P4 no-go.

## Source, oracle, and correction protocol

After activation, every primary source named by `0051/source-audit.md` will be
pinned by exact version and content hash before use. The review will check the
source theorem's variables, symplectic convention, domains, hypotheses, and
physical-to-reduced map rather than relying on the author's applicability
verdict.

The importable API and focused receipts will be audited only at their actual
predicate. Exact algebra will be independently rederived for the
Schwinger--Hopf identities and the `V_C/V_A` split. A test tally cannot prove
the Euler supplier, physical units, invariant subspace, or long gate-time
estimate. Existing receipts will be reused while their source and inputs stay
frozen; no unchanged oracle will be rerun merely for a count.

If one concrete defect is found, the review will state, in order, the strongest
supported result, the unsupported extension, the minimum correction, and the
evidence that would restore the stronger claim. Root may apply one bounded
correction, followed by one changed-statement/dependency check. No second
substantive pass will be performed.

## Verdict map and exclusions

The final review will assign separate route verdicts to:

1. physical energy--action normalization;
2. exact Schwinger--Hopf quotient geometry and area;
3. the complex-linear/antilinear decomposition;
4. realization of two noncollinear physical Pauli axes;
5. finite-gate leakage and physical-action drift; and
6. the actual Euler supplier and downstream applicability boundary.

Each route ends as `established`, `refuted` with a mechanism, or `blocked`
with the missing construction. The joined verdict will preserve the strongest
useful theorem these units support and leave the larger P4 construction active
unless its full positive contract is actually met.

The frozen review excludes active `0044`, `0046`, `0052`, and `0055` except
where an already reviewed statement is explicitly declared as a dependency.
It also excludes nonlinear branch or stability claims, autonomous particle
dynamics, quantum-state selection, Born/reset/detector dynamics, exchange
statistics, Lorentz covariance, finite propagation speed, electron or
neutrino identification, P4 completion or no-go, registry promotion, source or
central edits, and commits.
