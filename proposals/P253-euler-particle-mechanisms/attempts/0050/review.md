# P253/0050 joined independent review of 0046 and 0049

## Frozen transaction and correction provenance

This transaction independently reviews root-owned attempts `0046` and `0049`
as two separate verdict units. The reviewer authored and implemented neither
target, API, test suite, nor bounded correction. The frozen review README has
SHA-256
`1e0907acdb95e4acb8add50248703efb9104fadee786b7d7bb44557f21c74ccf`.
Central activation exited exactly `0`; its command, stdout, empty stderr, and
exit hashes are respectively
`5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
The accepted authority base is release `v0.183.0`.

The substantive pass requested one bounded precision package. Its receipt is
`0050/target-correction-receipt.md`, SHA-256
`03ed2b9d32e9af51177cb0e381aaf62c9b7af1eaf009d0176311990aa6c156fd`.
The single correction check confirms this complete map:

| Artifact | Frozen SHA-256 | Corrected SHA-256 |
|---|---|---|
| `0046/construction.md` | `744968cbedfaa2e1dfd424eec3e7d2ca8e98d303d481b76a88240bed66608b1f` | `7ee4b9d9679a9e2880b49a1efa4294a46f0c88be573124074dc39b55985a6d63` |
| `0046/result.yaml` | `7bc3acdf32661ea2dfc4eb8a23f10541bf7f3bffe647d8bb9c99724eb829d9ea` | `698841744a5ac7e1ed9a5e8a57fe3669cecbe31fe35b0338416293b303df5f69` |
| `0046/source-audit.md` | `77987e351b920ef74d0a06128f716a1546c07c00cc7c8c882a960ef4883a6147` | `74accea1e7b3a11410c8cf38c536936b36990ff517f61af47d57402cc3d671c2` |
| `0046/validation.md` | `bb85f0e02a4b5d05e1aaadb226505588f894b1f447a1fc1248548932a1e23582` | `c47ca3107b8fe6079f33758b0e3a62a2ddcab9ccb5a2c3619e2b237f15cac86b` |
| `0049/construction.md` | `d921309b780f44a8fee72ef723b6f11de9648a1806d71c34086edc68432a9897` | `b207eac0378998336d26af863b3467c6bbfc2f2cf9c7b6156427889c4b26a06c` |
| `0049/result.yaml` | `228cd30e6ad93eef33d647ef73c1a27498980623cb64d04b38a899dccc2478cc` | `b248708fa624b02bb478ae2ae02161ffc0f1b17a009120ac441dbaac40d3ca27` |
| `0049/source-audit.md` | `3eb37819c0c38767b46f5955a08f531fa19208e989d63cb7266f46f06bcdc970` | `f85f0dcd4aa831d010d05ce4ff0dbab9b7d50ae6974ba7369a9634430338c6f8` |
| `0049/validation.md` | `edc3cb99e471cb60e0af308817c2daaa8761e538f3457eae47619c131c82686a` | `282cc57296e77bbe2c3889a70bfc9311190d43090d1cc9c18a6d0ba32e096ae5` |

The correction changes no API or test. Its full positive-definite hypothesis,
CCR/compact-`su(2)` boundary, threshold convention, all-zero clock case, and
reset/first-exit distinction are all supported. The affected-text stale-claim
scan, target result YAML parse, terminal-newline checks, and scoped
`git diff --check` pass. No second substantive pass or test rerun was made.
Active `0048`, P4/LP4 completion, and particle claims remained outside scope.

## Accepted-source applicability

The source transfer is accurate at the accepted C-CST scopes in release
`v0.183.0`:

- `C-CST-008` supplies an actual compact two-dimensional isovortical tangent
  pair with `Omega(Q,S)=B!=0`, a positive-definite complete Euler Hessian, and
  the stated physical affine angle, `G`, and `L` rows. It expressly does not
  make the local angle a global compact symmetry or the tangent pair an
  invariant nonlinear mode.
- `C-CST-009` supplies a conditional phase Cauchy--Born pullback of the Euler
  material action and canonical optical branch. Its reconstruction residual
  is retained; unrestricted Euler invariance is not supplied.
- `C-CST-011` supplies actual finite-action two-phase Kelvin-prepared histories
  only on fixed compact time windows, with full velocity and pressure tails.
  It is not an exact monochromatic eigenmode or unrestricted rotor.
- `C-CST-017` supplies the periodic prepared nonzero-`K` linear response,
  physical observations, inherited action, and full current through
  `o(|K|^2)` on each fixed finite window. Preparation and source norms retain
  their declared dependence.
- `C-CST-018` transfers those forms and observations to one compact-ring
  same-field ensemble and exactly normalizes the initial physical forms. It
  remains a prepared finite-window linear-response/second-variation theorem,
  not an invariant manifold, all-time evolution, or detector.

The checked internal source anchors include the positive compact-pair proof
(`bd7e9c95...`), its independent C-CST-008 review (`658b8093...`), the local
normalizer (`13a00561...`), exact Gram construction (`dcde8a0f...`), current
transfer (`f33ba31b...`), compact joint theorem (`76db9845...`), and its
independent review (`347967d2...`). None supplies a quantum occupation rule,
temporal random clock, projective reset, exchange loop, universal action unit,
or finite signal speed.

The independently reviewed `0042/0045` theorem supplies only the exact local
retained state `W[v,w0]` and full-field pressure/tag history. It makes `w0`
load-bearing and does not turn it into a stochastic detector law.

## Unit A: compact-pair quantum-two-state test

### Exact KKS oscillator

On the accepted tangent plane, write

    Omega=B dq wedge dp,
    H_2=(Hq q^2+2N q p+P p^2)/2,

with `B!=0`, `Hq>0`, and `Hq P-N^2>0`. Thus the Hessian is positive definite.
For the inherited convention `i_X Omega=dH_2`, direct contraction gives

    A=(1/B) [[N,P],[-Hq,-N]],
    A^2=-(Hq P-N^2) I/B^2,
    omega_*^2=(Hq P-N^2)/B^2>0.

The sign, off-diagonal coefficient, and frequency are exact. This proves one
positive classical oscillator on the tangent block. It does not prove that
the plane is an invariant nonlinear Euler submanifold. On the compact
C-CST-018 field, the canonical physical initial forms are exact, while the
prepared optical histories retain their `o(|K|^2)` fixed-window error.

**KKS-plane verdict: established as an exact positive classical tangent
oscillator; autonomous invariant Euler dynamics is not supplied.**

### Affine Heisenberg algebra and finite CCR boundary

For the accepted rows

    theta=q,    G=(B^2/P)q,    L=Bp,

the inverse KKS form gives

    {theta,L}=1,
    {G,L}=B^2/P,
    {theta,G}=0.

The independent rederivation confirms these are constant central brackets.
The independent generators are the affine Heisenberg translation algebra,
not the state-dependent brackets of compact `su(2)`. The larger C-CST
observation map is simultaneous classical tomography, not a noncommuting
operator algebra.

If an external positive `hbar` and exact CCR are imposed on this noncompact
plane, oscillator quantization gives the infinite ladder
`E_n=hbar nu(n+1/2)`. An exact finite-dimensional CCR truncation is impossible:

    tr([Q,P]-i hbar I_d)=-i hbar d !=0.

The coordinate operators also carry `|1>` outside the span of `|0>,|1>`.
This refutes the exact finite-CCR truncation route for this plane. It does not
refute a different physical reduction to a compact KKS sphere with `su(2)`
brackets and a finite representation. That compact-orbit construction remains
an open positive route and would require its own actual Euler realization.

**Algebra verdict: affine Heisenberg structure and the exact finite-CCR
obstruction are established; a different compact-`su(2)` reduction is not
refuted.**

### Exchange and declared augmentations

The C-CST whole-state translation/O(3)/time-reversal law acts on field, cell,
tags, sources, and observer together. It is an ensemble symmetry, not an Euler
trajectory exchanging two carriers. Two prepared copies therefore supply no
collision-free invariant two-copy phase space, endpoint identification,
exchange holonomy, or selected `Z2` character.

The displayed stochastic variational system is algebraically equivalent to a
Schrodinger equation after adding a configuration probability density,
diffusion `D`, mass `m`, osmotic term, and `hbar_ext=2mD`. These are explicit
new laws, not consequences of deterministic incompressible Euler. Likewise,
Kelvin circulation rescales continuously as `Gamma -> (A/B) Gamma` under
Euler similarity. Conservation alone therefore supplies no universal
circulation/action quantum. This refutes bare topology-only selection across
the similarity class, while leaving background-selected or separately
axiomatized scales open.

**Exchange/augmentation verdict: actual two-copy exchange is blocked; the
stochastic construction is a coherent conditional extension; bare universal
circulation quantization is refuted by continuous scaling.**

### Unit A verdict

Unit A is **established after the bounded precision correction** at scope
`EXACT_CLASSICAL_KKS_HEISENBERG_AND_FINITE_CCR_BOUNDARY`. The accepted compact
pair contains an exact classical oscillator and physical current readout. It
does not by itself realize a quantum two-state system, but the result does not
close a separately constructed compact-`su(2)` carrier route.

## Unit B: analyzer and first-event detector

### Analyzer and exponential race

For a normalized `psi in C^2` and unitary analyzer `U`, unitarity gives the
nonnegative channel fractions

    I_i=|(U psi)_i|^2,    sum_i I_i=1.

This is exact energy-channel algebra conditional on an actual physical
same-carrier analyzer. The current C-CST preparation does not supply that
autonomous analyzer.

For independent exponential waiting times with rates
`lambda_i=kappa I_i`, `kappa>0`, and `Lambda=sum_i lambda_i>0`,

    P(T_i in dt and i first)=lambda_i exp(-Lambda t) dt,
    P(i first)=lambda_i/Lambda.

Hence normalized channels give `P(i first)=I_i`. A zero-rate channel never
wins. Continuous independent positive-rate clocks tie with probability zero.
On a finite window `[0,T]`,

    P(i first by T)=(lambda_i/Lambda)(1-exp(-Lambda T)),
    P(no click by T)=exp(-Lambda T),

so conditioning on a click preserves `lambda_i/Lambda`. If all rates vanish,
no finite event occurs and a conditional winner distribution is undefined;
the API correctly rejects that input rather than fabricating a winner.

**First-event verdict: the normalized infinite- and finite-window winner law
is exactly established conditional on the declared independent temporal
clocks and rate-current rule.**

### Threshold exposure, reset, and sequential transitions

For a detector threshold `lambda` with CDF
`F(x)=P(lambda<=x)`, the corrected event rule `lambda<=I` gives probability
`F(I)`. Obtaining `I` for every normalized intensity therefore requires the
uniform CDF or another law that derives the same function; determinism and
energy conservation alone do not select it. Actual exclusive basins and one
fixed invariant measure for every state and analyzer setting remain missing.

If an event additionally resets the state to output ray `e_i` and the clocks
are freshly and independently prepared, the next unitary analyzer `V` gives

    P(j|i)=|V_ji|^2.

This establishes the sequential transition table under explicit reset and
renewal hypotheses. It is not an autonomous consequence of reversible,
measure-preserving Euler dynamics on a finite invariant phase volume: an
invertible volume-preserving time map cannot collapse positive-volume basins
many-to-one onto projective rays or an attracting reset state. An autonomous
first-exit or scattering construction can preserve reversibility by retaining
escaped/environmental degrees of freedom, but the outcome channel alone does
not reprepare the carrier for a later trial. Sequential reset must therefore
be open-system/coarse-grained, externally prepared, or realized on an enlarged
state space whose discarded degrees remain explicit.

**Reset verdict: the conditional projective sequential law is established;
autonomous finite-volume Euler reset/repreparation is not. First-exit and
scattering remain distinct constructive routes.**

### Euler applicability, retained state, action, and propagation

The accepted C-CST Poisson marking is spatial: it selects good cells in a
stationary random field. It is not a temporal point process or a hitting-time
theorem. A same-substrate clock requires a material trajectory or flux map,
joint rare-event limit, intensity-proportional rates, independence/renewal,
and exclusive capture on the actual carrier. None is supplied.

Any detector coupling must evolve the exact retained state `(v,w,tau)` with
its load-bearing `w0`, compute pressure from the full field `v+w`, and derive
the analyzer intensities and basins on that same state. A finite moment list,
prepared target history, or signed material current cannot replace the
unresolved state or become a nonnegative detection rate by naming it one.

Finally, the common rate scale cancels from all conditional winner
probabilities and selects no action quantum. The probability law contains no
propagation theorem; full Euler pressure is elliptic. A finite event packet,
action scale, latency, or effective finite-speed band therefore remains a
separate physical construction.

**Euler-detector verdict: the deterministic threshold and rare-event Euler
routes are blocked by the named basin-measure and joint hitting/capture
constructions; action and finite-speed selection remain separate frontiers.**

### Unit B verdict

Unit B is **established after the bounded precision correction** at scope
`EXACT_CONDITIONAL_EXPONENTIAL_FIRST_EVENT_AND_RESET_THEOREM`. It preserves a
strong positive detector theorem: normalized channel energy fractions become
the exact first-event and sequential projective probability tables under the
declared analyzer, clocks, exclusive capture, reset, and renewal hypotheses.
The theorem identifies rather than hides every stochastic/open-system input.

## API and oracle scope

`euler_quantum_two_state.py` correctly implements the general KKS generator,
frequency square, affine brackets, finite-CCR trace obstruction, and Euler
circulation scaling. The focused receipt reports `5 passed in 1.61s`, exit
zero. The tests exercise diagonal positive Hessian data and domain mutations;
the independent review supplied the general off-diagonal derivation and full
positive-definite hypothesis. The API does not construct an invariant Euler
mode, a compact orbit, or a quantum representation.

`euler_measurement_bridge.py` correctly implements exactly normalized unitary
channel intensities, conditional exponential-race winner fractions with a
strictly positive total rate, and the assumed-reset transition matrix. Its
all-zero rejection matches the mathematical domain. The focused receipt
reports `5 passed in 1.59s`, exit zero. Finite-window/no-click and tie formulas
are proved analytically in the construction rather than implemented by the
API. The module does not construct temporal clocks, capture, reset, retained
Euler dynamics, or a detector measure.

Both test receipts were reused at unchanged code/test hashes and were not
rerun. They are finite algebraic regression evidence, not source-applicability
or physical-realization proofs.

## Joined verdict and next dependency

The joined transaction is **established after one bounded precision
correction** at scope
`EXACT_CLASSICAL_TWO_STATE_BOUNDARY_AND_CONDITIONAL_FIRST_EVENT_DETECTOR`.
Unit A and Unit B retain separate evidence and verdicts. Their conjunction
shows exactly what a classical Euler carrier already supplies and exactly
which detector hypotheses convert normalized channel energy into projective
probabilities.

The next decisive construction is a deterministic first-exit or scattering
detector on one persistent compact carrier: construct the physical analyzer,
derive an intensity-proportional joint hitting law from the same retained
Euler state including `w0` and full pressure, retain environmental/escaped
degrees for reversibility, and give an explicit repreparation or enlarged-state
renewal mechanism for sequential trials. A separately constructed compact
`su(2)` reduction remains the positive state-space route; it is not ruled out
by the finite-CCR result.

This review does not license P4 or LP4 completion/no-go, a finite-dimensional
autonomous quantum Euler system, Born probabilities without the declared
detector law, an Euler projective reset, exchange statistics, a universal
action quantum, finite propagation speed, a particle/electron/neutrino claim,
active `0048`, or parent-campaign completion.
