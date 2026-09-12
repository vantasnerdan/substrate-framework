# P253/0047: independent review of the 0043 quantum-sufficiency boundary

This README preregisters one fixed-boundary, non-author/non-implementer review
of root-owned attempt `0043`. The target is the strongest exact mathematical
content actually supported by five separate units: a Koopman representation,
conditional quantization of an Euler KKS sphere, two-carrier and Hopf-map
topology, a positive symplectic-mode complex structure, and bare Euler
kinematics. The units will receive separate verdicts before any joined verdict.

The useful target is a precise sufficiency and non-selection ledger. A
conditional mathematical route may be established even when Euler does not
construct its invariant measure, choose its prequantum integer or statistics
character, import CCR/CAR, or produce a Lorentz cone. Absence of those choices
is not a P4 no-go, a particle theorem, or evidence that no restricted Euler
sector could supply them in a later construction.

## Independence, ownership, and activation boundary

`particle-balance-review` did not author or implement `0043` or its API/tests
and owns only append-only review artifacts under `attempts/0047`. Root retains
ownership of `0043`, the canonical module/tests, central proposal and registry,
validation, memory, and commits. The reviewer will not edit `0043`, its
implementation/tests, active `0044` or `0046`, any central/generated file,
memory, or Git history.

Before this freeze, no `0043` README, construction, result, validation,
source-audit, command, stdout, stderr, module, or test body was opened. Only
filenames, content-blind hashes, central registration metadata, accepted
authority, Git metadata, and the user's stated review units were inspected. A
required memory search exposed only the already supplied high-level statement
that `0043` is an early P4 sufficiency test. Memory has no evidentiary role and
will not select the verdict.

No target body may be opened until root:

1. registers `0047` centrally with this independent owner and exact scope;
2. records the repository-schema invocation, stdout, stderr, and exit under
   `attempts/0047`; and
3. makes `attempts/0047/activation-schema.exit` contain exactly `0`.

After activation, the reviewer will pin this README and verify the target
hashes before one substantive pass. At most one bounded evidence-driven
correction pass may follow. The reviewer will then write only `0047/review.md`
and `0047/verdicts.yaml`. The accepted authority base is release `v0.183.0`;
the observed preregistration Git head is
`ca929ce584b462c033fde7866e6231076445ecd5`. The shared dirty worktree belongs
to its existing owners and will be preserved. The physics-skill preflight
completed with `ok=7`, `warn=0`.

## Frozen target inventory

The exact target attempt contains the following content-blind pinned files:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0043/README.md` | `3e4a8dd16a28a51f40676e148c3a92231338412f6abb30bca4fb20a680e53675` |
| `0043/construction.md` | `828b5d26b12dad5b63956181f608253b4c6c0c933ea0cc112ab962a029fa0e64` |
| `0043/result.yaml` | `5f4cfe8d6b140accb12ad2baa895404c67e8167c3e8bf57ce4b5f6b0c491593f` |
| `0043/source-audit.md` | `4d4a7ab43b057c487d0da788b5badf9ee33b66fb3207b736ca386070eb1e1415` |
| `0043/validation.md` | `88733482406946fee6613fa574b7fe8391f8637613f6e8500cee2ce7f89b7ead` |
| `0043/focused-pytest.command.txt` | `dd0baad1c919492417ce6ed66564aade29c992fe7e430db97bd2f015d13b2f18` |
| `0043/focused-pytest.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0043/focused-pytest.stdout` | `1b564d56b190fb6a550bb038ea4ced8533b5574ea4b3a63e6c41a456b4a2630b` |
| `0043/focused-pytest.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0043/activation-schema.command.txt` | `0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8` |
| `0043/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0043/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0043/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |

The importable algebraic surface is frozen separately:

| Artifact | Preactivation SHA-256 |
|---|---|
| `src/substrate_framework/euler_quantum_bridge.py` | `cc5035538ba233e255fe73b495b6b2ef53b924e69aadc9102c3f451c69708f8e` |
| `tests/test_euler_quantum_bridge.py` | `1e43e97a1b3fdcca8d3d25dc94382f38e9dd24713f02283b3027602682a6f6b6` |

Activation receipts are provenance, not scientific evidence. If activation
names another implementation, test, or source as part of this transaction,
its path and pre-review hash must be appended before it is opened. Otherwise
no later or neighboring artifact may enter by proximity.

## Frozen suppliers and non-inheritance

The following already-reviewed Euler artifacts are potential suppliers only
at their individually adjudicated scopes. Their existence does not establish
any `0043` conclusion, and the target must identify every proposition it
actually imports:

| Supplier artifact | SHA-256 |
|---|---|
| `0005/derivation.md` | `90e17874b928e6ccda19b9e745346836f1c4a2140ced8a3e434fbc44eb70abc2` |
| `0005/result.yaml` | `ae1b9fa001a9dfab37d884cd08a5a0ce39a57b1c33174290db1e0cd5ee02e3e5` |
| `0005/source-audit.md` | `2434db2e0d875cd2c917ecf166cc5c6c59bb9f780609d9764626e5dbf264e85e` |
| `0011/review.md` | `5c3a938abfb555c6d84a8ceb8e51bcd62a2896d873e25597a12bc833786940c5` |
| `0011/verdict.yaml` | `9aab64398d4af05175e78d954410edbd63e7cc908f63f5fbba95cfd52d7cb1c1` |
| `0017/global-carrier.md` | `2206f7ec86f4c4caca0f12bf7850a119d8050c85004bc8ecbbf7a7180fa3b0c3` |
| `0017/result.yaml` | `fa01baffae0d1c6e0d3ca5267e358368d3af128b096285a190b7610e9a8cd3a4` |
| `0018/review.md` | `0e7ef9c5ceb62dcbd562c3a6418ce3396e270127c3e01716868bd7794be03871` |
| `0018/verdict.yaml` | `5f8fd5b16fc18a900fc6eedfeb73e858fa7795a7123fae119c97b7f145909512` |
| `0027/solitary-wave-construction.md` | `c18d394a757f9fd967975ba23a486e62087c285c8fda4f08d8aeaa0159c6f850` |
| `0027/exterior-construction.md` | `0c2cda195dd5bea9163a7a6c582bfcff746a608a2e48d0a9cd7a007c1b260fc3` |
| `0027/result.yaml` | `dcef9066af6c83596af6fcba77d4ed1b0e679fc2b044ddbcc8664308460c0b33` |
| `0031/derivation.md` | `f92079c746fd66a1b182c85fb5da3f1b10dbe9772e2ae22ff6d094f5e778d12a` |
| `0031/hodge-domain-correction.md` | `7706af77edd143e12cc33ba57a5aff882f644bf802f14b0fc21bb8f6c26e30fd` |
| `0031/result.yaml` | `9041e385d75de3fce6a71eeff64c9cd8c5449be5c04c704ac9cb7d3a7411fb42` |

In particular, prior Euler action/KKS work may supply only its actual local or
conditional coadjoint-orbit form, periods, and finite relative observables.
It does not supply an invariant probability measure, geometric quantization,
an Euler-to-Hopf configuration-space inclusion, a statistics character, a
CCR/CAR representation, or relativistic kinematics.

After activation, `0043/source-audit.md` will be checked against the exact
primary statements it cites. Geometric-quantization results may be imported
only for their stated compact integral Kähler orbit. Koopman results require
the claimed measure-preserving flow. The Krusch--Speight/Finkelstein--
Rubinstein result must be checked on its actual based Hopf-map configuration
space, sector, boundary conditions, loop, and parity formula; it cannot be
transferred to Euler merely because both constructions use a sphere or a
linked field. Standard Galilean Euler invariance supplies no Lorentz theorem.

## Unit A: Koopman Hilbert representation and physical observables

1. Identify the exact Euler state space `X`, local or global flow `S_t`, sigma
   algebra, and probability measure `mu`. If `mu` is assumed rather than
   constructed for the actual carrier, the theorem is conditional. Invariance
   must be `mu(S_t^(-1)E)=mu(E)` on the times/domain claimed.
2. On `L2(X,mu)`, check the convention
   `U_t F=F composed with S_t` or its inverse, unitarity, group versus local
   evolution, strong continuity, and the sign/domain of its generator. A
   finite-dimensional or stationary measure example does not establish a
   generic invariant measure for three-dimensional Euler.
3. Define physical observables as real measurable functions and their
   multiplication operators. Verify self-adjointness domains for unbounded
   observables and exact commutativity of multiplication operators. Separate
   commuting classical observables from Koopman evolution operators.
4. State the maximum implication: an invariant measure yields a classical
   Hilbert representation and commuting physical-observable algebra. It does
   not by itself supply a quantum state space, Born rule, noncommutative
   observables, superselection, or a physically selected measure.

**Unit A verdict choices:** established as stated, established conditional on
an explicitly assumed invariant measure, or blocked with the missing actual
Euler measure/flow domain named. None is a P4 no-go.

## Unit B: KKS sphere, period, Chern class, and conditional spin

5. Verify that the reviewed Euler family supplies an actual two-dimensional
   coadjoint orbit with stabilizer and orientation giving `S2`, rather than
   only two tangent vectors or a local chart. Track density, profile,
   amplitude, and physical action units in the pulled-back KKS form.
6. Recompute the global symplectic period, including sign and the full `2pi`
   and `4pi` factors. For the convention
   `Omega=s sin(theta) dtheta wedge dphi`, check

       integral_S2 Omega=4pi s,
       c1=[Omega/(2pi hbar)],
       N=2s/hbar in Z.

   If the target uses another orientation or normalization, derive its exact
   equivalent rather than matching labels.
7. Audit the compact Kähler quantization theorem actually invoked. For
   positive integer `N`, holomorphic quantization of `CP1` with the degree-`N`
   line bundle has dimension `N+1` and carries the irreducible `SU(2)`
   representation of spin `N/2`. Check the metaplectic/half-form convention if
   used; it can shift representation formulas and must not be silent.
8. Distinguish conditional integrality and representation sufficiency from
   dynamical selection. Euler scaling may change the KKS action continuously;
   it must be derived with all density factors. Neither Euler nor an algebraic
   API selects `hbar`, `N=1`, spin one-half, a polarization, or a quantum
   probability interpretation unless an additional construction does so.

**Unit B verdict choices:** exact conditional geometric quantization at the
actual Euler KKS normalization, a narrower abstract `S2` theorem without an
Euler inclusion, or a blocked sphere/global-period bridge. A finite axial
moment alone is not quantum spin.

## Unit C: unordered carriers, Hopf topology, and FR characters

9. Define the ordered two-carrier configuration space, collision exclusion,
   center/scale/orientation variables, and identical-carrier quotient. The
   reduction of a pure nonzero relative displacement by exchange gives
   `(R3 minus {0})/(r equivalent -r)`, homotopy equivalent to `RP2`, with
   fundamental group `Z2`; verify that no extra field variable contracts the
   loop in the space actually claimed.
10. To promote this abstract two-point calculation to Euler, require a
    continuous map into the admissible same-field Euler configuration space,
    collision-free along the exchange, compatible endpoint identification,
    and proof that the exchange loop remains nontrivial after inclusion. A
    sketch of two separated carriers or a quotient of labels is not such a
    proof.
11. Audit the exact Krusch--Speight Hopf-map theorem from the primary source:
    based/unbased mapping space, Hopf sector, compactification and finite-energy
    boundary, fundamental group, rotation/exchange loop class, and its charge
    or parity dependence. Credit only the theorem's actual quantifiers.
12. Require an explicit Euler-to-Hopf map with preserved sector and loop class
    before transferring the FR result. The two one-dimensional characters of
    `Z2` permit bosonic and fermionic signs; topology alone does not select a
    character or prove a relativistic spin-statistics theorem.

**Unit C verdict choices:** established abstract `RP2`/FR sufficiency,
established actual Euler inclusion and nontrivial loop, or blocked with that
inclusion/character selection named. A failed inclusion is not a no-go for a
different carrier or topological encoding.

## Unit D: positive symplectic mode versus CCR/CAR

13. Starting from the actual Euler tangent space, KKS form `Omega`, and
    quadratic Hamiltonian `H_2`, verify nondegeneracy after quotienting all
    stabilizer/neutral directions and positivity on the stated mode. Check the
    Hamiltonian generator, eigenfrequency signs, and normalization.
14. Derive the compatible complex structure rather than naming multiplication
    by `i`. It must satisfy `J^2=-I`, preserve `Omega`, and make
    `g(x,y)=Omega(x,Jy)` positive in the convention used. Check whether the
    construction is global on the claimed positive subspace or only one
    finite-dimensional oscillator block.
15. Separate classical Kähler/one-particle data from quantization. A positive
    symplectic mode supplies complex amplitudes and a compatible inner
    product; CCR, a bosonic Fock representation, CAR, a fermionic Fock
    representation, vacuum selection, and operator domains are additional
    imports. CAR is not obtained by changing a sign in the classical KKS form.
16. If topology and mode quantization are joined, verify their common physical
    state space and character. An `RP2` sign and an oscillator complex
    structure on unrelated spaces do not select fermions or a full quantum
    theory.

**Unit D verdict choices:** established positive classical complex structure,
conditional chosen CCR/CAR representation, or blocked actual-mode positivity.
The strongest supported classical result must be preserved.

## Unit E: Galilean Euler dispersion and absence of a Lorentz cone

17. Fix the Euler domain and background. About rest, the Leray-projected
    transverse symbol has zero temporal frequency; about a uniform background
    `U_0`, it is the Galilean advection branch
    `omega=U_0 dot k`. A nonzero uniform background on whole space has infinite
    absolute energy and must be treated periodically or relative to the
    background.
18. Audit every energy-momentum transformation or dispersion formula with its
    density, volume, relative-energy, and boundary hypotheses. Adding a
    constant Galilean boost to a finite-energy whole-space velocity generally
    changes the asymptotic class, so an absolute finite-energy orbit cannot be
    inferred without a relative formulation.
19. Check the exact API algebra for the Galilean branch and any comparison to
    `omega^2=c^2|k|^2` or `E^2=c^2|P|^2+m^2c^4`. Bare incompressible Euler has
    no finite invariant signal speed or Lorentz mass shell from these
    equations. Do not confuse two transverse spatial polarizations with
    positive/negative temporal branches.
20. Scope the conclusion correctly: bare Euler kinematics does not supply the
    Lorentz cone required by P4. This is a missing sufficiency construction,
    not proof that no special medium, carrier sector, homogenized limit, or
    additional dynamics can have emergent Lorentz symmetry.

**Unit E verdict choices:** established Galilean/no-bare-Lorentz boundary,
qualified relative-background formula, or a correction of sign/domain. It
cannot close P4 negatively.

## Source, API, and oracle audit

21. Inspect every citation in `0043/source-audit.md` after activation. Pin the
    exact version, theorem or proposition, mathematical space, hypotheses,
    conventions, and conclusion. Secondary summaries may locate a result but
    cannot expand it. In particular, do not import Hopf-soliton topology into
    Euler or geometric quantization into a nonintegral/local orbit by analogy.
22. Inspect every public symbol in `euler_quantum_bridge.py`: exact inputs,
    dimensions, positivity and integer domains, sign and `2pi` conventions,
    parity/character output, Galilean convention, error handling, and explicit
    exclusions. Determine whether it implements algebra only or silently
    claims a measure, field-space inclusion, quantization, or Lorentz dynamics.
23. Map every test assertion to one mathematical predicate. Require derived
    checks for KKS period/Chern number, spin and Hilbert dimension,
    `RP2`/`Z2` character algebra, compatible-complex-structure identities, and
    Galilean versus Lorentz dispersion only where the API actually exposes
    them. Representative mutations are a missing factor two, reversed
    symplectic orientation, noninteger/negative Chern input, collapsed parity
    character, wrong `J^2` sign, and a Lorentz cone substituted for advection.
24. The focused receipt is regression evidence only after its source hashes,
    command, assertion predicates, and exit are inspected. Exact API algebra
    cannot prove existence of an invariant Euler measure, a global Euler KKS
    sphere, a nontrivial Euler exchange loop, a selected statistics character,
    CCR/CAR, or emergent relativistic dynamics.

No production numerical quantity, soft eigenvalue, force, energy splitting,
or stability threshold is registered, so `small-ratio-numerics` does not bind.
The strongest oracle is direct functional analysis, differential-form period
and Chern calculation, configuration-space topology with exact source
applicability, symplectic linear algebra, and Euler symbol/boost derivation,
with the importable API credited only as corroborating exact algebra.

## Verdict structure and exclusions

Units A--E receive separate `route_verdict`, `evidence_scope`, strongest
supported statement, unsupported extension, and upgrade path. The joined
result may be `established` as a sufficiency/non-selection ledger even if one
or more stronger Euler realization routes remain conditional or blocked. A
missing selection principle is not a contradiction. `Refuted` is reserved for
an explicit counterexample under the stated hypotheses.

The review will preserve the strongest true statement and request at most one
bounded correction pass supported by concrete evidence. It excludes active
`0044` and `0046`, P2, a completed P4 theorem or P4 no-go, a finite-dimensional
autonomous particle model, quantum probability or measurement dynamics,
universal spin one-half, fermion statistics selection, electron/neutrino
identity, Lorentz-covariant dynamics, stability, interactions, and parent
campaign completion.
