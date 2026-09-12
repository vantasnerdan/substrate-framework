# P253/0050: joined independent review of 0046 and 0049

This README preregisters one fixed-boundary, non-author/non-implementer review
transaction with two separately adjudicated units. Unit A reviews root-owned
attempt `0046`, the compact-pair quantum-two-state test. Unit B reviews
root-owned attempt `0049`, the analyzer and first-event detector construction.
The joined purpose is to preserve the strongest exact or conditional result of
each route without letting either unit inherit evidence from the other.

The useful positive target is precise. Unit A may establish exact classical
KKS/Heisenberg dynamics and a conditional prepared two-state construction even
if no invariant compact quantum sector follows. Unit B may establish an exact
conditional first-event probability and reset theorem even if Euler has not
yet supplied the temporal clock, exclusive capture, action scale, or finite
propagation law. Missing bridges remain named construction frontiers; they are
not a P4 no-go, LP4, or a particle theorem.

## Independence, ownership, and activation boundary

`particle-balance-review` did not author or implement `0046`, `0049`, their
APIs, tests, source audits, or validation receipts. It owns only append-only
review artifacts under `attempts/0050`. Root retains ownership of both target
attempts, canonical modules/tests, central proposal and registry, validation,
memory, and commits. The reviewer will not edit either target, either API/test
pair, active `0048`, any central/generated file, memory, or Git history.

Before this freeze, no `0046` or `0049` README, construction, result,
source-audit, validation, command, stdout, stderr, module, or test body was
opened. Only filenames, content-blind SHA-256 hashes, central registration
metadata, accepted authority, Git metadata, and the user's declared review
scope were inspected. Repository memory exposed only the already declared
high-level route summary; it has no evidentiary role and will not select either
verdict. The local GitNexus index is stale at `e38a8e9` relative to current head
`591a642`, so it supplies no review evidence.

No target body may be opened until root:

1. registers `0050` centrally with this independent owner and exact joined
   scope;
2. records the repository-schema invocation, stdout, stderr, and exit under
   `attempts/0050`; and
3. makes `attempts/0050/activation-schema.exit` contain exactly `0`.

After activation, the reviewer will pin this README, verify every frozen hash,
and conduct one substantive pass. The primary sources declared by the two
source audits will then be inventoried by exact version and local or retrieved
content hash before their propositions are used. At most one bounded,
evidence-driven correction pass may follow. The reviewer will write only
`0050/review.md` and `0050/verdicts.yaml` after activation.

The accepted authority base is release `v0.183.0`, source baseline
`b39f392f8feee21c93a6052a5f971bc2df862779`, and observed preregistration
head `591a642`. The shared dirty worktree belongs to its existing owners and
will be preserved. The physics-skill preflight completed with `ok=7`,
`warn=0`.

## Frozen target inventory

Unit A is frozen content-blind at these hashes:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0046/README.md` | `e36dcc9812f4a6f896d30a2ffe5af28e1f15d23010b7a1befb3c4099dd4e7161` |
| `0046/construction.md` | `744968cbedfaa2e1dfd424eec3e7d2ca8e98d303d481b76a88240bed66608b1f` |
| `0046/result.yaml` | `7bc3acdf32661ea2dfc4eb8a23f10541bf7f3bffe647d8bb9c99724eb829d9ea` |
| `0046/source-audit.md` | `77987e351b920ef74d0a06128f716a1546c07c00cc7c8c882a960ef4883a6147` |
| `0046/validation.md` | `bb85f0e02a4b5d05e1aaadb226505588f894b1f447a1fc1248548932a1e23582` |
| `0046/focused-pytest.command.txt` | `b30b2d5d2169f51708025eb018e6bf8d89c8f80c08be9bad1712527951e1981b` |
| `0046/focused-pytest.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0046/focused-pytest.stdout` | `bb898df9e07ec11703c009f863b6b7daad87d5158e37c1e77a008e0a3e081f64` |
| `0046/focused-pytest.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0046/activation-schema.command.txt` | `0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8` |
| `0046/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0046/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0046/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `src/substrate_framework/euler_quantum_two_state.py` | `fed1fe8a0b181cd00db27be5dfae8da29c00b49fbe052d73992e156f8a24836a` |
| `tests/test_euler_quantum_two_state.py` | `992b45164ce7feed1cbb9d86662c3a04c70e46987ea9070a63521898070ceba0` |

Unit B is frozen content-blind at these hashes:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0049/README.md` | `7cfd498c683b503d025bea07c269b20a7e541764738d316375f705d347a5708a` |
| `0049/construction.md` | `d921309b780f44a8fee72ef723b6f11de9648a1806d71c34086edc68432a9897` |
| `0049/result.yaml` | `228cd30e6ad93eef33d647ef73c1a27498980623cb64d04b38a899dccc2478cc` |
| `0049/source-audit.md` | `3eb37819c0c38767b46f5955a08f531fa19208e989d63cb7266f46f06bcdc970` |
| `0049/validation.md` | `edc3cb99e471cb60e0af308817c2daaa8761e538f3457eae47619c131c82686a` |
| `0049/focused-pytest.command.txt` | `fbe7fef772cf0faff140427185ec8b5454c1c347b77de7ae17d611fe131bb864` |
| `0049/focused-pytest.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0049/focused-pytest.stdout` | `d5579e263d3cf70acb717b91007b4e8f4a22f00eeaa6558ee55446d3ad9d35d1` |
| `0049/focused-pytest.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0049/activation-schema.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0049/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0049/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0049/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `src/substrate_framework/euler_measurement_bridge.py` | `89d459c0121459584a8cab39f0db8e1bec2c20a6ef161c98551e4d51f34f4180` |
| `tests/test_euler_measurement_bridge.py` | `19d3a4723b8bbfb8e8f4150118a25444b225df53d496a8f1f5c052d7dff00b5d` |

Activation and pytest receipts are provenance until their commands, terminal
outputs, and scientific predicates are audited. No test count or API presence
will be allowed to establish a source theorem or a physical Euler bridge.

## Non-inheritance and source applicability

The units receive independent verdicts. Unit A's exact classical structure
does not supply Unit B's random clock, capture, reset, or probability measure.
Unit B's conditional probability theorem does not supply Unit A's KKS orbit,
action, quantization, invariant subspace, or statistics. A joined narrative
will be accepted only as the conjunction of the two separately earned scopes.

After activation, every import named in each target will be reduced to the
exact proposition it supplies. In particular:

- the C-CST primary source inventory will be pinned by exact bibliographic
  version and content hash, then checked for its actual phase space,
  symplectic convention, flow or recurrence hypotheses, and quantization
  scope;
- previously reviewed Euler KKS, action, retained-state, pressure, or carrier
  results may be used only at their individually adjudicated domains;
- a theorem about abstract canonical variables, stochastic processes, or a
  spatial point process is not imported as an actual Euler temporal detector;
  and
- an API returning the desired finite algebra is corroboration, not a proof of
  its physical realization.

## Unit A: compact-pair quantum-two-state test

### Exact classical KKS dynamics

1. Identify the actual compact-pair phase space and the KKS two-form inherited
   from the cited source and Euler construction. Recompute every sign, density,
   circulation, `2 pi`, action, and orientation factor.
2. Derive the Hamiltonian vector field on the claimed plane from the stated
   convention `i_X Omega=dH` or its declared alternative. Check the rotation
   generator, angular frequency, time orientation, conserved radius/action,
   and whether the plane is global, local, reduced, or merely a tangent block.
3. Check that the variables are dynamically accessible Euler directions, not
   only algebraic labels. Record any stabilizer, gauge, domain, or compactness
   quotient needed to make the two-form nondegenerate.

### Heisenberg structure versus compact `su(2)`

4. Derive the Poisson brackets of the affine coordinate and momentum rows.
   Determine whether their bracket is a nonzero central constant, giving the
   Heisenberg algebra, rather than the state-dependent closing brackets of
   compact `su(2)`. Similar dimension or two labels cannot identify the Lie
   algebra.
5. Audit any compact two-level truncation against the exact dynamics. A plane
   oscillator has noncompact phase space and ordinary CCR quantization has an
   infinite ladder unless an additional compact orbit, constraint, quotient,
   or nonlinear closure is constructed.
6. Apply the finite-dimensional trace identity directly: for finite matrices,
   `tr([Q,P])=0`, so exact `[Q,P]=i hbar I` is impossible when `hbar` is
   nonzero. Distinguish an approximate/projection commutator from an exact CCR
   representation and from a compact `su(2)` algebra.

### Dynamics, exchange, and added structure

7. Separate a prepared finite-window two-amplitude evolution from an invariant
   two-dimensional subspace of the autonomous full Euler evolution. Check the
   forcing, preparation, endpoints, leakage/complement, pressure response, and
   time interval actually proved.
8. Check whether two copies form only a direct product or whether the target
   constructs a collision-free same-field exchange path with endpoint
   identification and nontrivial topology. A copy permutation does not by
   itself establish exchange statistics or a fermionic character.
9. Classify stochastic noise, temporal clocks, circulation quantization,
   action units, reset, or capture postulates as explicit augmentations unless
   derived from the literal Euler state and flow. Such conditional additions
   may support a precise theorem but cannot be relabeled as bare Euler.

**Unit A verdicts will be separate:** exact classical KKS oscillator and
frequency; affine Heisenberg rather than compact `su(2)`; exact finite-CCR
obstruction; prepared finite-window versus invariant dynamics; actual
two-copy exchange; and every stochastic/circulation augmentation. The unit's
strongest true positive statement will be retained even if the quantum
two-state interpretation remains conditional or blocked.

## Unit B: analyzer and first-event detector

### Analyzer normalization and first-event law

10. Define the analyzer channels, state norm, unitary basis change, and
    intensities. Verify nonnegativity and exact normalization, including the
    zero-norm input boundary and all zero/nonzero channel cases.
11. Independently derive the exponential race. For independent clocks with
    rates `lambda_a`, compute the survival law, first-event density, and
    `P(a first)=lambda_a/sum_b lambda_b` when the total rate is positive. If
    `lambda_a=gamma I_a`, confirm the normalized intensity fractions and show
    where the common rate factor cancels.
12. Audit zero-rate channels, the all-zero case, finite observation windows,
    no-event probability, conditioning on an event, and ties. Independent
    continuous exponential clocks have zero tie probability; a code tie rule
    is not a derivation for discrete or deterministic clocks.

### Reset, sequential transitions, and deterministic exposure

13. Check the declared post-event reset map for normalization, phase handling,
    idempotence, and its status as hypothesis rather than Euler consequence.
    With fresh independent clocks, derive the next analyzer's transition
    probability from the overlap intensity. Identify any independence or
    Markov assumption required for a sequential product law.
14. Expose deterministic-threshold claims exactly. If a hidden threshold `R`
    with CDF `F` triggers at intensity `I`, then the event probability is
    `F(I)`. The Born value `I` follows only from a specified uniform threshold
    on the normalized interval or another derived law; it is not automatic
    from thresholding.
15. Distinguish an actual deterministic first-event map from ensemble
    frequencies. State the measure over initial unresolved state, thresholds,
    phases, or detector microstates that makes a probability claim meaningful.

### Euler applicability and missing physical bridges

16. Distinguish spatial Poisson marks from temporal exponential clocks. A
    spatial intensity law becomes a temporal hazard only after a derived
    trajectory/flux map, sampling law, stationarity or renewal property, and
    exclusivity rule. Dimensional agreement or the same exponential formula
    does not supply that bridge.
17. Audit the analyzer observable on the literal retained Euler state. Track
    dependence on the unresolved initial datum `w0`, reconstruct pressure from
    the full field, and check whether the analyzer acts on the same evolved
    state. Finite material summaries cannot silently replace `w0` or full
    pressure.
18. Keep the action and signal-speed questions separate from the probability
    algebra. A rate constant does not select `hbar` or a finite action packet;
    an elliptic Euler pressure solve does not establish a finite propagation
    cone. Any locality, detector separation, latency, or finite-speed
    hypothesis must be explicit.
19. Audit the importable API at its actual algebraic boundary: normalization,
    zero-channel handling, finite-window first-event probabilities, reset,
    sequential transitions, deterministic CDF exposure, and invalid-domain
    behavior. Tests must expose sign, normalization, zero-rate, and conditioning
    mistakes without being credited as a physical Euler construction.

**Unit B verdicts will be separate:** normalized analyzer intensities; exact
conditional exponential first-event theorem; finite-window/no-event and tie
boundaries; conditional reset/sequential theorem; deterministic-threshold CDF
exposure; spatial-to-temporal applicability; retained-state/full-pressure
closure; and independent action/finite-speed frontiers. The strongest target
is the exact conditional detector theorem, not an unconditional measurement
or particle claim.

## Oracle and correction policy

The preferred oracle for both units is independent exact derivation plus direct
source-hypothesis checking. Finite API tests are regression evidence for that
algebra. The objective bridge for each predicate will be recorded as

    computed API or symbolic predicate
      -> exact algebraic proposition
      -> source and Euler phase-space licenses
      -> maximum conditional or exact verdict
      -> effect on the P253 measurement/quantum frontier.

The review will check one representative wrong sign or normalization and every
load-bearing domain edge already exercised by the frozen tests. It will reuse a
scientifically adequate recorded receipt at unchanged hashes and will not
re-run tests merely for counts. A rerun is warranted only if a correction
changes a predicate or the receipt cannot support its claimed boundary.

At most one bounded correction may repair a concrete error in signs, factors,
quantifiers, domains, source attribution, API behavior, or evidence role. The
correction check will inspect only those changes and directly affected edges;
it will not become a second substantive pass. The result will preserve the
strongest true claim and name the construction that would restore any stronger
version.

## Frozen joined verdict boundary

The final record will issue separate Unit A and Unit B route verdicts before a
joined disposition. The joined maximum is a conjunction of independently
earned statements, expected to have the form
`EXACT_CLASSICAL_TWO_STATE_BOUNDARY_AND_CONDITIONAL_FIRST_EVENT_DETECTOR` only
if the bodies and sources support it. No unit may upgrade the other.

The review cannot license P4 completion or no-go, LP4, a physical particle,
electron or neutrino identity, a genuine finite-dimensional autonomous quantum
system, an actual Euler exchange-statistics theorem, a Born-rule derivation
without the declared detector measure, finite-speed propagation, active
`0048`, or parent-campaign completion.
