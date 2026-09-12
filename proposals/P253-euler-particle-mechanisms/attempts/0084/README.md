# P253/0084: independent review of final 0080 Cao HSE and charged continuation

## Frozen review object

This README preregisters one content-blind, fixed-boundary, independent
non-author/non-implementer review of final root-owned `0080`. It reviews only
the exact steady Green-map, bordered-complement, charged-IFT, and conditional
local-parameter claims made by that attempt. It has six separately
adjudicated units:

- **A — global nonlinear map:** whole-space zero extension, the `X/G_1`
  compact-Fredholm construction, `C^2` positive-part/tag regularity, and `C^1`
  Maxwell elimination;
- **B — soft-cell complement:** `Q A Q` on its exact left/right scaled spaces;
- **C — limiting exhaustion and physical matrices:** allowed-growth
  Lane--Emden alternatives and the `q_a/q_t` matrices for both `B_I` and the
  frozen-target `B_R` system;
- **D — bordered isomorphism:** index zero, injectivity and surjectivity,
  including parameter-dominated normalized sequences;
- **E — finite-`L` charged continuation:** the subluminal charged IFT and the
  exact `U(1)`-extension versus bare-Euler boundary; and
- **F — conditional refinements:** the formal/source-derived leading Schur
  jet and a local `epsilon` carrier path, both only at the hypotheses actually
  proved.

Evidence and verdicts do not transfer between units. In particular, failure
of a dynamical Fredholm inverse is not automatically failure of a steady
elliptic bordered map, and a nonzero finite Schur determinant cannot replace
the infinite-dimensional complement theorem. The review will not demand
invertibility with both `(mu,c)` artificially fixed: the physical augmented
unknown/constraint system and its actual bordered inverse are the object.

The quantitative geometric-`delta` coverage and downstream graph/Riesz/
projector continuity needed by a rational-ray IVT are not `0080` claims and
are outside this review. Fixed-`(kappa,R)` `B_R` existence is adjudicated only
through Units C--E at the target's own scope. Analyzer control, P4, particle
mechanics, action or charge quantization, and electron/neutrino conclusions
are excluded. Active `0083` is excluded in full; none of its bodies, sources,
modules, tests, receipts, results, or conclusions may be opened or used.

## Independence and activation boundary

`particle-balance-review` authored or implemented none of `0080`. It owns
only append-only review artifacts under `attempts/0084` and will not edit
`0080`, source/API/test files, central/generated records, memory, or Git
history.

Before this freeze no `0080` file body, API body, test body, verifier body,
receipt text, manifest text, or captured output was opened. Inspection was
limited to registry/Git metadata, filenames, and content-blind SHA-256 hashes.
The final `author-completion-receipt.md` and `artifact-hashes.sha256` are the
sole authoritative package selectors; intermediate failure or execution
receipts are frozen only as provenance and cannot supersede them. The complete
`0083` tree remained unopened.

The target may be opened only after root:

1. centrally registers `0084` with this exact reviewer and scope;
2. records the repository-schema command, stdout, stderr, and exit under
   `attempts/0084`; and
3. makes `attempts/0084/activation-schema.exit` contain exactly `0`.

After activation the reviewer will recheck every frozen hash, open only the
final 0080 completion receipt, manifest, and source audit before the
derivation, pin each primary-source version and theorem, and then inspect the
frozen `0080` bodies/API/tests/verifier. One substantive pass may request at
most one evidence-driven bounded correction package, followed by one
correction-only check. Final writes are limited to `0084/review.md` and
`0084/verdicts.yaml`, followed by a report to `w4:p2`.

The authority base at freeze is accepted release `v0.183.0`.
`governance/releases/current.yaml`, `governance/claims.yaml`, and the P253
proposal registry have SHA-256 values
`7b889a3c0186c0c22e698fd4ccc69e19066b8262018696cfc41e573ec2f0ffd8`,
`78dbe694359739ee15f7db2b2cb1ccde7f48b5af34f6d9cd451282732dd08983`,
and
`17d7c93ba403c76e3392f9ece9866e2d826507baca39bac5af19453d4e3df62a`.
The observed repository head is
`1cdf81becee9f7ef828fa651cc363377dfc2932b`. The shared worktree contains
unrelated active changes, all preserved and excluded.

This is an analytic elliptic/Fredholm/operator review. No production
numerical design or run is frozen. Exact algebra and focused unit tests may
support formulas and API behavior, but they cannot prove Schauder mapping,
Fredholm surjectivity, a nonzero physical determinant, or a conditional local
parameter path.

## Frozen final 0080 inventory

The final package selectors match the coordinator-supplied prefixes:

- `0080/author-completion-receipt.md` SHA-256
  `83d28c1b51e70127bb6fcbdc58a3c37f949133d8b2cd1f35eff84bb07e92a9fc`;
- `0080/artifact-hashes.sha256` SHA-256
  `1b3bbe9593d1a034105487aa70d0e2eb25ddea8b5e0870d29ba27e99bd719faa`.

The complete top-level target inventory frozen content-blind is:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0080/README.md` | `27f71a514462a5eb9ea01cf4b88e919a55b6ea3917858b878da781cf00c4280b` |
| `0080/activation-schema.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0080/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0080/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0080/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0080/artifact-hashes.sha256` | `1b3bbe9593d1a034105487aa70d0e2eb25ddea8b5e0870d29ba27e99bd719faa` |
| `0080/author-completion-receipt.md` | `83d28c1b51e70127bb6fcbdc58a3c37f949133d8b2cd1f35eff84bb07e92a9fc` |
| `0080/derivation.md` | `f36d577639f07fbd649e9b97f450dafe4d67c667a82d554d62f5668dba4958f2` |
| `0080/first-green-schur-failure.md` | `724d52e2262788ea3c24db5555086543b6e182942b0dc0684e9feaee2ae28c64` |
| `0080/focused-tests.command.txt` | `a1eb29d425dc207eacb92511c372fd99fee8da895bcdbcfc2fc96dee467808ca` |
| `0080/focused-tests.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0080/focused-tests.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0080/focused-tests.stdout.txt` | `e0b7eca228d15acdf8482428f5f26811f29a6401362ef4ab8aa3c4a14ffe5bf2` |
| `0080/green-bordered.command.txt` | `d4f043a245184cdb95e3990be464b22b863b3a53afb297dcadc4d29e21747477` |
| `0080/green-bordered.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0080/green-bordered.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0080/green-bordered.stdout.txt` | `601816ad3de2218e31d33f9fe9fcca191e8a65be06b5a25a6cd3b0a8893f3c83` |
| `0080/green-schur-corrected.command.txt` | `d4f043a245184cdb95e3990be464b22b863b3a53afb297dcadc4d29e21747477` |
| `0080/green-schur-corrected.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0080/green-schur-corrected.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0080/green-schur-corrected.stdout.txt` | `8cb7befc5a2fc6e0e55c4262a64abf8b82c5d915ad457eb5ea37ff93f26413d4` |
| `0080/green-schur.command.txt` | `d4f043a245184cdb95e3990be464b22b863b3a53afb297dcadc4d29e21747477` |
| `0080/green-schur.exit` | `4355a46b19d348dc2f57c046f8ef63d4538ebb936000f3c9ee954a27460dd865` |
| `0080/green-schur.stderr.txt` | `ebfa7b1950804398ca5ed258749ce348001c9b00b8a0913a223bc9862578ecf5` |
| `0080/green-schur.stdout.txt` | `d028e45fb23ed27071378269478cd7721316ddacf2ad4fce576b2cfa2a9159ca` |
| `0080/result.yaml` | `b4e9f4d2dd69bc322182ed3f564ce5190b85b1da9ac6019d810cc1ee114ef8e3` |
| `0080/source-audit.md` | `0f5aff991e941de0a4bc7a9e9b4638f9f57dbec71084b50d55a8bf0c4c88d4d7` |
| `0080/validation.md` | `e4fbdfa1fd59bb66656ad1ffb061f0325812ee419dd60c4899e33e9828e72352` |
| `0080/verify_green_schur.py` | `575554224922341826c08138d7a4bcc1b4a91d5a8baad5dca9adee67ac19f3fa` |

The importable API and focused consumer test frozen content-blind are:

| Artifact | Preactivation SHA-256 |
|---|---|
| `src/substrate_framework/euler_cao_schur.py` | `92a621bb6537fbb1178062a26b91309fbc1583c82aefa6486ce6cb045d9d8f44` |
| `tests/test_euler_cao_schur.py` | `feaaae535265d5b8456bbb52e379f2f012f32b3252d2bc72739481be65ab4495` |

The sorted top-level `0080` filename manifest has SHA-256
`34ac1c97fc04122c349c107fd6bae7929e235cd3d98a108ddf20b25f4044c5ab`.
No differently named source module or test is presumed.

## Source and authority audit

After activation, every primary source actually used will be checked at its
exact version and proposition. The review will record the source equation,
domain/codomain, boundary and decay hypotheses, parameter set, interface
regularity, normalization, and exact conclusion transferred. In particular,
Cao pointwise existence, uniqueness, limiting-kernel nondegeneracy, and an
auxiliary asymptotic parameter system are distinct statements; none may be
silently promoted into an exact augmented inverse, a nonzero physical Schur
determinant, or smooth fixed-radius coverage.

## Unit A — global zero extension, compact Fredholm map, and regularity

1. Freeze the exact whole-space axisymmetric source space `X`, the Green
   operator `G_1`, the zero-extension convention across the positive-part
   interface, all axis/infinity weights, and the finite moment/constraint
   rows. Prove every asserted `X -> X` mapping and compactness statement.
2. Identify the nonlinear map `M:X times R^d -> X` (with the target's actual
   parameter dimension). Distinguish regularity of an individual Cao profile
   from Fréchet regularity of `M`. Check each source, Green, tag/action and
   parameter composition, including any `M R` or equivalent reconstruction,
   at the claimed second order.
3. Reconstruct Cao's positive-part interface bootstrap without a `C^1`
   shortcut: uniform `W^(2,q)` first, then `C^(1,alpha)`, then Hölder
   regularity of the right-hand side, and only then the source-space Schauder
   conclusion. Verify that Cao Lemma 3.9 is used in its actual source space and
   does not itself supply a stronger Schauder or Banach-map theorem.
4. Check the Nemytskii differentiability at the actual `p`, smooth
   zero-extension of the tag, and the `C^1` Maxwell elimination on the same
   whole-space class. A charged Green block may not silently upgrade the
   hydrodynamic map's regularity.
5. Derive the `I-compact` or equivalent Fredholm structure from the actual
   global map. Track whether compactness survives every finite row and the
   axis/exterior tail.

**Unit A verdict:** separately adjudicate the global compact-Fredholm
construction, `C^2` hydrodynamic/tag map, and `C^1` Maxwell elimination.

## Unit B — `Q A Q` soft-cell complement

6. Define the exact left and right scaled spaces, their norms and embeddings,
   the soft-cell projection `Q`, and the operator `A`. Do not identify spaces
   merely because their blow-up limits share a formula.
7. Reconstruct the complement estimate for `Q A Q`, including support,
   core-interface, axis, whole-space Green, and finite-row terms. Track every
   epsilon weight in both directions.
8. Verify that the projected inverse maps the declared right space back to
   the declared left space with uniform control, and that the projection is
   bounded in both. A scalar Lane--Emden estimate does not alone prove the
   finite-core complement theorem.
9. Separate the soft translation cell from true symmetries and bordered
   parameter directions. No inversion with both physical parameters fixed is
   demanded.

**Unit B verdict:** state whether the exact finite-core `Q A Q` complement is
established on the two-sided scaled spaces.

## Unit C — allowed-growth exhaustion and `q_a/q_t` matrices

10. Classify every allowed blow-up limit under the actual scaled-space growth
    bound, not only bounded Lane--Emden kernel elements. Prove that no constant,
    logarithmic, affine, harmonic, or source-generated alternative survives
    the matching and finite rows.
11. Track separately profile-dominated and parameter-dominated normalized
    sequences. State which limiting translation/profile/scale modes remain
    before bordering.
12. Derive the physical `q_a/q_t` matrices, their row/column conventions,
    signs, scales and determinants for both the impulse-target `B_I` system and
    the frozen-radius-target `B_R` system.
13. For `B_R`, retain fixed physical `(kappa,R)` and identify exactly which
    chemical-potential/speed or other parameters remain unknown. Do not
    substitute a fixed-speed source family for the frozen-target bordered
    system.

**Unit C verdict:** adjudicate allowed-growth exhaustion and the two physical
finite matrices separately; neither may inherit the other's determinant.

## Unit D — index-zero bordered injectivity and surjectivity

14. Assemble the exact bordered operator from Unit B and the appropriate Unit
    C matrix. Prove index zero on the declared product spaces.
15. Prove injectivity for every normalization regime, including sequences in
    which parameter coordinates dominate while the profile component tends to
    zero. A profile-normalized contradiction alone is incomplete.
16. Prove surjectivity, either directly or from index zero plus established
    injectivity, with all target finite rows. A left inverse or a nonzero
    formal determinant alone is insufficient.
17. Check uniform bounds, dependence on finite core scale, and compatibility
    with Unit A's nonlinear-map derivative. Adjudicate the full augmented
    inverse, not a fixed-`(mu,c)` subblock.

**Unit D verdict:** determine whether both bordered injectivity and
surjectivity are established for `B_I` and `B_R` at their claimed scopes.

## Unit E — finite-`L` subluminal charged IFT

18. Freeze the finite-`L` parameter regime, subluminal inequality, signed
    coupling convention, tag normalization, and exact uncharged base. Verify
    that the Maxwell elimination is `C^1` and enters the reduced fluid map at
    order `g^2` on the Unit-A spaces.
19. Apply the IFT only with the Unit-D bordered inverse matching the chosen
    physical targets. Track circulation, impulse or radius, center, speed,
    support/interface, axis, affine Coulomb tail, and charge rows.
20. State the actual neighborhood and why it remains finite-`L`; no
    `L -> infinity`, all-frequency, or subluminal-threshold passage may be
    inferred.
21. Keep the theory boundary exact. A charged branch belongs to the declared
    transported-tag `U(1)` Euler--Maxwell extension. It is not a solution of
    bare unforced Euler and supplies no action/charge selection, stability,
    analyzer, P4, or particle theorem.

**Unit E verdict:** adjudicate the finite-`L` charged branch only in its
declared extension and physical target slice.

## Unit F — conditional leading Schur jet and local path

22. Reconstruct the formal/source-derived leading triangular coefficient
    `S_lead`, including its normalization and the value-level `c`, `mu`, and
    `I` asymptotics. Decide whether any exact finite-core nonvanishing follows;
    do not presuppose a sharp finite-core theorem. Keep the `C^1`
    differentiated remainders and Lyapunov--Schmidt tangent identity open
    unless they are proved on the declared spaces.
23. Distinguish a value-level determinant sign from differentiability of a
    carrier parameter. If a local `epsilon` path is asserted, identify the
    additional `C^1` parameter dependence and nonzero derivative hypothesis.
24. Keep every local-path statement explicitly conditional unless the
    derivative/open-mapping row is proved on the same `B_R` spaces. A leading
    `delta=C epsilon[1+o(1)]` value asymptotic does not by itself give
    `delta'(epsilon)!=0`.
25. Do not require or infer quantitative geometric-`delta` coverage,
    common-domain graph/Riesz continuity, or a rational-ray IVT. Those are
    downstream active work outside final `0080` and excluded with `0083`.

**Unit F verdict:** adjudicate the conditional leading Schur jet, preserve the
blocked sharp finite-core sign/differentiated boundary, and classify the local
`epsilon` path at its exact conditional boundary.

## API, oracle, and final-artifact discipline

After activation, the importable API and focused tests will be mapped to their
actual predicates. The review will inspect whether they encode exact Green,
scaled complement, bordered, Schur, and finite-`L` algebra, and whether a
sign, normalization, rank, or parameter-dominated mutation is exposing.
Passing tests cannot prove the analytic interface, Banach-map, Fredholm, or
conditional local-path rows. No unchanged production suite will be rerun for
a count.

The reviewer will preserve the strongest supported theorem and request at
most one bounded evidence-driven correction if a finite repair is exposed.
Final `review.md` and `verdicts.yaml` give Units A--F separate route verdicts
and name the next construction for every open unit. No charged, analyzer, P4,
particle, or campaign-level conclusion is licensed.
