# P253/0086: independent review of final 0085 charged-Hessian and Maxwell spectrum

## Frozen review object

This README preregisters one content-blind, fixed-boundary, independent
non-author/non-implementer review of final root-owned `0085`. It reviews the
actual constrained Euler--Maxwell Hessian, tangent, and linear spectral
statements in seven separately adjudicated units:

- **A — field square and constraint Schur:** the exact electromagnetic square,
  Gauss and `div B` constraint reduction, signs, physical coefficients, and
  the field-eliminated Schur remainder;
- **B — infinite negative DA sector:** the claimed infinite-dimensional
  dynamically accessible negative subspace at one common nonzero coupling
  `g`, including finite rows and tag constraints;
- **C — infinite positive Maxwell sector:** the claimed infinite-dimensional
  compactly supported positive electromagnetic subspace on the constrained
  tangent space;
- **D — tag locking:** the exact passive-tag tangent/stabilizer formulas and
  their regular-band, support, and nonisochronous boundaries;
- **E — Maxwell essential spectrum:** the Weyl-sequence argument for the
  `i R` essential inclusion of the full constrained linearized generator;
- **F — Doppler shell and dark current:** the shell/coarea calculation and the
  exact necessity, not sufficiency, of a dark-current condition; and
- **G — interpretation boundary:** what the strong saddle, continuum
  resonance, possible bound state in the continuum, and charged-extension
  status do and do not imply.

Evidence and verdicts do not transfer between units. In particular, an
indefinite constrained Hessian does not alone prove nonlinear instability;
an `i R` essential inclusion does not by itself exclude an embedded
eigenmode; and a radiation-shell compatibility condition does not turn the
charged extension into bare Euler. The review preserves the strongest useful
positive theorem supported by the actual construction and treats every open
statement as a named next construction rather than a parent no-go.

Active `0087` is excluded in full. Its bodies, sources, APIs, tests,
verifiers, receipts, results, and conclusions may not be opened or used.
No P4/P5 completion, electron/neutrino identity, particle theorem, universal
radiation theorem, Born/reset dynamics, or nonlinear instability conclusion
is in scope.

## Independence and activation boundary

`particle-balance-review` authored or implemented none of `0085`. It owns
only append-only review artifacts under `attempts/0086` and will not edit
`0085`, source/API/test files, central/generated records, memory, or Git
history.

Before this freeze no `0085` file body, API body, test body, verifier body,
receipt text, manifest text, or captured output was opened. Inspection was
limited to registry/Git metadata, filenames, and content-blind SHA-256 hashes.
The final `author-completion-receipt.md` and `artifact-hashes.sha256` are the
package selectors; intermediate failed and corrected receipts are preserved
only as provenance until their roles are read after activation. The complete
`0087` tree remained unopened.

The target may be opened only after root:

1. centrally registers `0086` with this exact reviewer and scope;
2. records the repository-schema command, stdout, stderr, and exit under
   `attempts/0086`; and
3. makes `attempts/0086/activation-schema.exit` contain exactly `0`.

After activation the reviewer will recheck every frozen hash, open only the
final 0085 completion receipt, manifest, and source audit before the
derivation, pin every primary-source version and proposition, and then inspect
the frozen bodies/API/tests/verifier. One substantive pass may request at most
one evidence-driven bounded correction package, followed by one
correction-only check. Final writes are limited to `0086/review.md` and
`0086/verdicts.yaml`, followed by a report only to `w4:p2`.

The authority base at freeze is accepted release `v0.183.0`.
`governance/releases/current.yaml`, `governance/claims.yaml`, and the P253
proposal registry have SHA-256 values
`7b889a3c0186c0c22e698fd4ccc69e19066b8262018696cfc41e573ec2f0ffd8`,
`78dbe694359739ee15f7db2b2cb1ccde7f48b5af34f6d9cd451282732dd08983`,
and
`fa80c1cba0dfe5c3be147e012d4a401444302e6334cb44294de46eac3f7c718c`.
The observed repository head is
`1cdf81becee9f7ef828fa651cc363377dfc2932b`. The shared worktree contains
unrelated active changes, all preserved and excluded.

This is an exact variational, constraint, operator, and spectral audit. No
production numerical design or run is frozen. Exact algebra and focused unit
tests can support the field-square, Schur, dispersion, and coarea identities,
but cannot prove global DA admissibility, infinite-dimensional subspaces,
essential spectrum, absence of embedded modes, or nonlinear dynamics by a
passing tally.

## Frozen final 0085 inventory

The content-blind final package selectors are:

- `0085/author-completion-receipt.md` SHA-256
  `755173029492d7b5a6c81ade04460027fd8558c22bbbe3bb329fd0f4c153b69e`;
- `0085/artifact-hashes.sha256` SHA-256
  `09452de8228d977355efb657bbfe707c6a95a221a2fc55f38bdeedb99acdcce5`.

The complete top-level target inventory frozen content-blind is:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0085/README.md` | `874522770f9272d73e538435d016ce7a6f1dae7a60217c3ff77f33f7acc07847` |
| `0085/activation-schema.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0085/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0085/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0085/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0085/artifact-hashes.sha256` | `09452de8228d977355efb657bbfe707c6a95a221a2fc55f38bdeedb99acdcce5` |
| `0085/author-completion-receipt.md` | `755173029492d7b5a6c81ade04460027fd8558c22bbbe3bb329fd0f4c153b69e` |
| `0085/derivation.md` | `da2bd59eb3e26a17fa81955d889bc3e0b2b1f0137e8c6a00b3061e7c173a5082` |
| `0085/exact-check.command.txt` | `b08f9bf1783a1dfd87ca0d362bf929b0bb01f4304b2b30e02fb5de9469c7add3` |
| `0085/exact-check.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0085/exact-check.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0085/exact-check.stdout.txt` | `e8a220226cebc124c59d5d9a47487c5d183a62dc90bd1f7579c89bbe7862ccd9` |
| `0085/exact-v2.command.txt` | `b08f9bf1783a1dfd87ca0d362bf929b0bb01f4304b2b30e02fb5de9469c7add3` |
| `0085/exact-v2.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0085/exact-v2.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0085/exact-v2.stdout.txt` | `da4adb00fe224bd2004e8a18f34013d1e92402954c82567f49dfe270bafb81cb` |
| `0085/first-focused-failure.md` | `0bda46a3c4c75c3bf7936555c529fd887faa3b8fe4de42c0493ba09f8635af35` |
| `0085/first-focused.command.txt` | `c92800bfbcbac78742e81b2a895e7a628c540a374a3f896ee23cb7428586448c` |
| `0085/first-focused.exit` | `4355a46b19d348dc2f57c046f8ef63d4538ebb936000f3c9ee954a27460dd865` |
| `0085/first-focused.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0085/first-focused.stdout.txt` | `e40069fa1a5c5988d1731ded0bd091b5221444c0072df573f923c4f9ec118dfa` |
| `0085/focused-tests.command.txt` | `c92800bfbcbac78742e81b2a895e7a628c540a374a3f896ee23cb7428586448c` |
| `0085/focused-tests.exit` | `4355a46b19d348dc2f57c046f8ef63d4538ebb936000f3c9ee954a27460dd865` |
| `0085/focused-tests.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0085/focused-tests.stdout.txt` | `defc632b9f27aba0c5f87f33e99a42603b6007a49f9ba0e4d9590d3ab5c3233f` |
| `0085/focused-v2.command.txt` | `c92800bfbcbac78742e81b2a895e7a628c540a374a3f896ee23cb7428586448c` |
| `0085/focused-v2.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0085/focused-v2.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0085/focused-v2.stdout.txt` | `f256f2c5e7cd4ab60ca538c5d39efb2f1b67c27714fc8dcba62479825b8f0567` |
| `0085/focused-v3.command.txt` | `c92800bfbcbac78742e81b2a895e7a628c540a374a3f896ee23cb7428586448c` |
| `0085/focused-v3.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0085/focused-v3.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0085/focused-v3.stdout.txt` | `0e0177d95805cd7e9ccb1859b562b7fe1d8ffaed73559d85239e8afb6ed7a331` |
| `0085/postactivation-spectral-scope-receipt.md` | `f7af4be104e354534c01cc62f12d335d5c54f1fa9233fd4e0843f45281fcc34c` |
| `0085/repository-validation.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0085/repository-validation.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0085/repository-validation.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0085/repository-validation.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0085/result.yaml` | `e86304d7682b4fb532b8c3324077f1c7259e1b5a9f7e9841d54f4330a18c295b` |
| `0085/shell-first.command.txt` | `c92800bfbcbac78742e81b2a895e7a628c540a374a3f896ee23cb7428586448c` |
| `0085/shell-first.exit` | `4355a46b19d348dc2f57c046f8ef63d4538ebb936000f3c9ee954a27460dd865` |
| `0085/shell-first.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0085/shell-first.stdout.txt` | `96273a5365be82e920b0de215b1cfefbf60832e8e6436540a219e4688da15d1a` |
| `0085/source-audit.md` | `82f7f130154a7048ef4fb986ff10d4737eb9572e7945056e9e754f15a12c29b5` |
| `0085/validation.md` | `21721bcd9c1dcfb2200840d56bec615def159dbce1280ad5be0eeaf8d4b89579` |
| `0085/verify_charged_hessian.py` | `b055628193247da57958d2e3d5f12beea8cf8443fe920e96eb12d9a47e3f01ea` |

The frozen importable API and focused consumer test are:

| Artifact | Preactivation SHA-256 |
|---|---|
| `src/substrate_framework/euler_charged_hessian.py` | `c17fc51ec7f98051b082432b9dc9f2d807685c1d91c704ae799f45eef037ac30` |
| `tests/test_euler_charged_hessian.py` | `5acb87a2a9140a8b42830a2d3ec625e74b86a9f2a4ecb5e77b9d95b0de9a0fd6` |

The sorted top-level `0085` filename manifest has SHA-256
`a4371016b3ee0f788dec6acd898c27d7d6034aadbeaa1cbf51c26082dc6f5366`.
No differently named module, test, or external source is presumed.

## Source and authority audit

After activation, every primary source actually used will be checked at its
exact version and proposition. The review will record the source operator,
domain, constraint/gauge space, decay and moment hypotheses, essential-spectrum
convention, and exact theorem transferred. Accepted or independently reviewed
suppliers will be reused only at their pinned statements; adjacency to an
Euler, Maxwell, KKS, or Weyl theorem will not fill a missing bridge.

## Unit A — electromagnetic square and Gauss/`div B` Schur

1. Reconstruct the complete physical second variation, including density,
   `epsilon_EM`, `mu_EM`, signed coupling, comoving terms, pressure/constraint
   multipliers, and all fluid--field cross terms.
2. Complete the electromagnetic square directly and identify its exact
   positive field norm. Check signs by expanding the square back to the
   original Hessian rather than accepting a named Schur complement.
3. Freeze the linearized Gauss and `div B` constraints, their function spaces,
   affine Coulomb/finite-energy boundary conditions, gauge quotient, and any
   zero-mode or first-moment rows. Verify that the eliminated potentials solve
   those exact constraints.
4. Derive the residual fluid/tag Schur form after field minimization. State
   whether it is finite rank, compact, or merely bounded, and keep any
   conditional coercivity or sign estimate at its proved scope.

**Unit A verdict:** adjudicate the exact square and constrained Schur identity
separately from any global sign conclusion.

## Unit B — common-`g` infinite DA negative subspace

5. Freeze one sufficiently small but common nonzero `g`; a sequence with
   `g` depending on the mode number does not prove an infinite subspace for
   one charged carrier.
6. Construct every negative direction as a genuine dynamically accessible
   Euler--Maxwell/tag tangent on the same base. Track divergence, compact
   support or decay, circulation, impulse, center, charge, tag, gauge, and
   stabilizer rows.
7. Prove negativity for arbitrary nonzero finite linear combinations, not
   merely one direction at a time. Check support separation, cross terms,
   normalization, closure, and whether a uniform small-`g` bound survives as
   the dimension grows.
8. Distinguish an infinite negative subspace, infinite Morse index, and an
   unbounded-below form. Record exactly which statement is earned.

**Unit B verdict:** decide whether one common charged carrier has the claimed
infinite-dimensional DA negative subspace.

## Unit C — infinite positive compact Maxwell subspace

9. Define compactly supported field variations satisfying linearized Gauss,
   `div B`, gauge, and every finite constraint row, with fluid and tag
   components chosen consistently.
10. Verify strict positivity of the full constrained Hessian on every
    nonzero finite combination, including comoving and background-current
    cross terms. Disjoint supports alone are insufficient if a nonlocal
    constraint solve or Schur term recouples them.
11. Prove infinite dimensionality in the actual energy/domain topology and
    distinguish compact field seeds from compact potentials or compact
    constraint-corrected representatives.

**Unit C verdict:** adjudicate the infinite positive compact Maxwell sector
independently of Unit B.

## Unit D — passive tag-lock formulas

12. Derive the full linearized material-tag equation and the KKS/DA tangent
    formula on the charged base. Track the difference among an arbitrary tag
    variation, an advected variation, and a stabilizer tangent.
13. Recheck every division by the toroidal-vorticity profile or its action
    derivative. State the exact regular band where both are nonzero, and keep
    the center, flat/critical levels, support edge, exterior, and tangential
    stabilizer freedom explicit.
14. Verify that finite charge normalization and the chosen tag profile are
    preserved by the tangent construction. Local stabilizer locking must not
    be promoted to all-time global tag persistence.

**Unit D verdict:** preserve the exact local tag-lock theorem at its genuine
regular-band scope and classify every passive-tag direction left outside it.

## Unit E — Maxwell `i R` essential inclusion

15. Freeze the full linearized constrained generator, Hilbert/graph domain,
    gauge quotient, and essential-spectrum convention. Separate the free
    Maxwell principal part from compactly supported carrier and finite-rank
    constraint terms.
16. Construct graph-domain Weyl singular sequences for every claimed
    `i omega`, satisfying Gauss and `div B` exactly and converging weakly to
    zero. Check normalization, supports escaping the carrier, polarization,
    pressure/fluid/tag components, and the comoving Doppler shift.
17. Prove the perturbation vanishes on the sequence in the operator graph norm
    and that the chosen Weyl/Fredholm criterion yields the stated inclusion.
    Do not infer equality, absence of point spectrum, or a spectral gap unless
    separately proved.

**Unit E verdict:** adjudicate only the claimed `i R` essential inclusion of
the actual constrained generator.

## Unit F — Doppler shell, coarea, and dark-current necessity

18. Reconstruct the Fourier dispersion and the resonant Doppler shell with
    all `c_EM`, comoving-speed, sign, polarization, and measure factors.
19. Apply the coarea formula on the exact shell, checking regularity of its
    defining function, exceptional zero/threshold points, Jacobian, and the
    function space in which the radiation or resolvent norm is evaluated.
20. Identify the exact Fourier transform of the compact effective current
    entering each physical transverse polarization. Prove which vanishing
    trace is necessary for a localized time-harmonic mode.
21. Distinguish vanishing on one shell, all polarizations, and an open family
    of shells. A necessary dark-current condition is not existence of a bound
    state in the continuum, sufficiency for localization, or generic
    nonexistence.

**Unit F verdict:** adjudicate the exact shell/coarea calculation and retain
the dark-current statement only at its necessary-condition scope.

## Unit G — strong saddle, resonance, BIC, and theory boundary

22. State exactly what Units B and C imply for the constrained quadratic form:
    infinite positive and negative subspaces may establish a strong saddle or
    infinite Morse/co-Morse index, but not nonlinear instability, exponential
    growth, or absence of a center manifold without a dynamical theorem.
23. State exactly what Unit E implies for an internal mode. Essential
    continuum at its frequency creates a resonance problem; it does not alone
    remove an eigenfunction or prove radiation damping.
24. Preserve the BIC escape route. Decide whether Unit F gives only a
    necessary dark-current condition or a complete solvability theorem, and
    name the additional outgoing-resolvent/matching construction needed for
    sufficiency.
25. Keep the ontology exact: all charged statements belong to the declared
    transported-tag Euler--Maxwell extension. They neither become bare Euler
    nor establish charge/action selection, persistence, analyzer dynamics,
    P4/P5, or a particle conclusion.

**Unit G verdict:** give the strongest exact strong-saddle and continuum
statement while preserving resonance/BIC and bare-Euler boundaries.

## API, oracle, and final-artifact discipline

After activation, the importable API and focused tests will be mapped to their
actual predicates. The review will inspect whether they encode genuine square,
constraint, Schur, signature, Weyl-dispersion, and shell/coarea algebra, and
whether wrong signs, density factors, coupling parity, polarization, or shell
Jacobians are exposing. Passing algebra cannot prove infinite-dimensional
domain constructions, weak-null graph sequences, or nonlinear consequences.
Historical failed runs will be classified by mechanism rather than counted as
scientific evidence.

The reviewer will preserve the strongest supported theorem and request at most
one bounded evidence-driven correction if a finite repair is exposed. Final
`review.md` and `verdicts.yaml` give Units A--G separate route verdicts and name
the next construction for every open unit. No campaign-level exhaustion or
completion is licensed.
