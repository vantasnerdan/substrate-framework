# P253/0091: independent review of final 0090 shell flux and switching theorem

## Frozen review object

This README preregisters one content-blind, fixed-boundary, independent
non-author/non-implementer review of final root-owned `0090`. The target is
adjudicated in six separate units whose evidence and verdicts do not transfer:

- **A — conserved neutral source and outgoing normalization:** net-neutral
  comoving continuity, the outgoing resolvent sign and physical prefactor,
  and correct real-phasor treatment without double counting;
- **B — shell geometry and trace:** the full shell gradient, coarea measure,
  trace domain, and the exact Gaussian curl-current oracle;
- **C — conserved switching and late field:** the
  `a J+a_prime K` conservation switch, the full-real-current definition of
  `F_T`, and the exact late free-field spectral norm;
- **D — long-time and endpoint energy:** the fixed-ramp,
  integrable-first-correlation law `T P+O(1)`, the sharp weak-switch endpoint
  energy, and the comoving world-tube relative flux;
- **E — reciprocal decay convention:** the conditional factor-two relation
  between amplitude and energy/action decay rates; and
- **F — interpretation boundary:** the explicit absence of an actual
  carrier current, KKS normalization, gate, resonance, action selection,
  Born/reset/statistics mechanism, or electron/neutrino result.

The review will preserve the strongest exact and conditional positive theorem
supported by the actual construction. It will audit the weighted outgoing
resolvent, trace, time-limit, and flux bridges directly; a passing algebraic
oracle will not stand in for those analytic steps. It will distinguish an
exact source-to-radiation theorem from application of that theorem to a still
unconstructed physical Cao/analyzer current.

Active `0088` is excluded in full. The completed `0089` review and correction
transaction are also excluded: no `0083`, `0087`, `0088`, or `0089` body, receipt,
API, test, verifier, result, or conclusion may be opened or used in this
review. No P4/P5 completion, particle model, electron/neutrino identity,
measurement rule, or stability conclusion is in scope.

## Independence and activation boundary

`particle-balance-review` authored or implemented none of `0090`. It owns only
append-only review artifacts under `attempts/0091` and will not edit `0090`,
source/API/test files, central/generated records, memory, or Git history.

Before this freeze no `0090` file body, API body, test body, verifier body,
receipt text, manifest text, or captured output was opened. Inspection was
limited to registry/Git metadata, filenames, and content-blind SHA-256 hashes.
The excluded `0088` and `0089` review/correction surfaces remained unopened.

The target may be opened only after root:

1. centrally registers `0091` with this exact reviewer and scope;
2. records the repository-schema command, stdout, stderr, and exit under
   `attempts/0091`; and
3. makes `attempts/0091/activation-schema.exit` contain exactly `0`.

After activation, the reviewer will recheck every frozen hash, open only the
final completion receipt, manifest, source audit, and precision receipt before
the derivation, and then inspect the frozen API, tests, verifier, and captured
outputs. One substantive pass may request at most one evidence-driven bounded
correction package, followed by one correction-only check. Final writes are
limited to `0091/review.md` and `0091/verdicts.yaml`, followed by a report only
to `w4:p2`.

The authority base at freeze is accepted release `v0.183.0`.
`governance/releases/current.yaml`, `governance/claims.yaml`, and the P253
proposal registry have SHA-256 values
`7b889a3c0186c0c22e698fd4ccc69e19066b8262018696cfc41e573ec2f0ffd8`,
`78dbe694359739ee15f7db2b2cb1ccde7f48b5af34f6d9cd451282732dd08983`,
and
`a59fdcce299df25b83b37ed4799f839f3ed1beb7d61ef8040c2bdd5e90187274`.
The observed repository head is
`76f0bc37f8130645685731ad89d1317b4afcee4e`. The shared worktree contains
unrelated active changes, all preserved and excluded.

This is an exact PDE, Fourier-resolvent, trace/coarea, conservation, and
large-time audit. No production numerical design or run is frozen. Exact
symbolic checks may validate shell algebra and Gaussian transforms but cannot
prove weighted limiting absorption, late-time convergence, time-integrated
flux, or applicability to an unspecified physical current.

## Frozen final 0090 inventory

The content-blind final package selectors are:

- `0090/author-completion-receipt.md` SHA-256
  `9ff5b80f0855fc0aa74a9c33ba5d7b9dd21fdfc7bc9060c338f03b370e75ee0e`;
- `0090/postactivation-energy-units-receipt.md` SHA-256
  `5737d9d5dc81b761af6853fe46cf65657bbb44931557078f1692c94ed205e67c`;
- `0090/artifact-hashes.sha256` SHA-256
  `53acca1d5d7a54e6246c686c3c8794f7be52fe90fcb31fc802bcbd97b8edc3e4`.

The complete top-level target inventory frozen content-blind is:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0090/README.md` | `073fe1a19f64ec1590a6dc2920795c9b5f865bb018519214025eb6da18587d79` |
| `0090/activation-schema.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0090/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0090/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0090/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0090/artifact-hashes.sha256` | `53acca1d5d7a54e6246c686c3c8794f7be52fe90fcb31fc802bcbd97b8edc3e4` |
| `0090/author-completion-receipt.md` | `9ff5b80f0855fc0aa74a9c33ba5d7b9dd21fdfc7bc9060c338f03b370e75ee0e` |
| `0090/derivation.md` | `20c882b24bc788ea4bb6c962d89b789efbec1e19300b0a35e76cfd5c4a11e802` |
| `0090/exact.command.txt` | `dd2bec300a7bc7504704adc88ff2022058c2256004d7bdcd66efee95d32ff7b1` |
| `0090/exact.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0090/exact.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0090/exact.stdout` | `a0b6d501ba6ab0b56b148522b2f6cabdd3d1962db799c67d8df584f5cd4a44d2` |
| `0090/focused.command.txt` | `6f07c5449f506c8c460bec79dc81ff598177b7874c5c4d5b62620b7b46da221e` |
| `0090/focused.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0090/focused.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0090/focused.stdout` | `c5f70ee4aeebea0dedd14d505d0bf885b820f25812eec4e571303a08c9fc1a3f` |
| `0090/postactivation-energy-units-receipt.md` | `5737d9d5dc81b761af6853fe46cf65657bbb44931557078f1692c94ed205e67c` |
| `0090/preenergy-activation-schema.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0090/preenergy-activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0090/preenergy-activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0090/preenergy-activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0090/repository-validation.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0090/repository-validation.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0090/repository-validation.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0090/repository-validation.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0090/result.yaml` | `f64defb0400dc348ff68ced7356630cbec4717677a307728262d47230e66bd25` |
| `0090/source-audit.md` | `dc1450eb9ec25b14797a3f04c30c1ecad7650e42b4edbdef29c03f6a2314f60b` |
| `0090/validation.md` | `869722796ad79c8d25437d0b85f1765878a539a5ad6ce32a933da47ca53e403b` |
| `0090/verify_shell_flux.py` | `b536d374ab54c9a7f54931a3644334162bc6627ea294aa0e8514faa4a74522d1` |

The frozen importable API and focused consumer test are:

| Artifact | Preactivation SHA-256 |
|---|---|
| `src/substrate_framework/euler_maxwell_radiation_flux.py` | `74881415a94e53350d2db5edb328a0d8611d3eace9b1c48fb680f0c0db14fb4d` |
| `tests/test_euler_maxwell_radiation_flux.py` | `1bc141cc20da1f5cb8c628a3f714da51e6c0b1428fbb5c33fd50304907d618d9` |

The sorted top-level 0090 filename manifest has SHA-256
`3a0f161b05b4f96d6922478240af6dda5f745d0bb254069a56c7e9427c429469`.
No differently named source module, test, verifier, or external source is
presumed.

## Source and authority audit

After activation, every primary source actually used will be checked at its
exact version and proposition. The review will freeze Fourier-transform,
phasor, Maxwell-energy, Poynting, comoving-frame, and limiting-absorption
conventions before transferring any standard identity. Independently reviewed
framework suppliers will be reused only at their pinned statements. A source
that proves stationary outgoing resolvents will not silently supply a switched
time-domain limit or a carrier-specific current.

## Unit A — neutral continuity and outgoing normalization

1. Reconstruct charge/current continuity in the laboratory and comoving
   coordinates, including every translation term and the spatial integral.
   Verify the exact decay/first-moment assumptions under which net neutrality
   follows and no boundary flux is lost.
2. Freeze the Fourier and time conventions. Derive the outgoing `+i0` or
   `-i0` prescription from retarded Maxwell evolution rather than selecting a
   sign by analogy.
3. Reconstruct all electromagnetic constants, Fourier factors, Jacobians,
   polarization projections, and time-average factors in the outgoing-power
   prefactor.
4. Represent a real current with its conjugate frequency partner. Verify that
   the phasor convention neither counts both signs twice nor drops the factor
   distinguishing complex amplitude, cycle-averaged energy, and real-field
   energy.

**Unit A verdict:** adjudicate neutrality, outgoing sign, and physical
prefactor as independent load-bearing rows.

## Unit B — shell coarea and exact Gaussian curl current

5. Derive the complete gradient of the Doppler denominator with respect to
   all three Fourier variables, not only its radial derivative. Prove it is
   nonzero on each claimed shell component in the stated subluminal,
   nonzero-frequency domain.
6. Apply coarea with the full Euclidean gradient and reconcile any radial
   parameterization Jacobian with the surface measure. Check shell components,
   orientation independence, endpoints, and all physical prefactors.
7. State the precise Sobolev/weighted domain whose Fourier transform has a
   shell trace, including the declared `<x> J in L2` first-moment hypothesis
   and its exact trace consequence. Compact smooth sources may make this
   immediate; a generic finite-energy current does not.
8. Independently transform the declared Gaussian curl current, verify exact
   divergence freedom and zero total charge/current where claimed, apply the
   transverse projection, and evaluate its shell functional without inserting
   the expected answer.

**Unit B verdict:** separate the exact shell geometry and exposing oracle from
the analytic trace theorem used by a general source.

## Unit C — conserved switching and exact late free field

9. Starting with a conserved monochromatic pair `(rho,J)`, derive a switching
   current of the form `a J+a_prime K`. Verify the sign and exact spacetime
   divergence cancellation, support/decay of `K`, net neutrality, and behavior
   at both switch endpoints. Keep the zero-mean hypothesis on `rho` distinct
   from the weighted-domain hypothesis on `P_T K`; neither substitutes for the
   other.
10. Define `F_T` from the full real switched current, including negative
    frequencies and cross terms. Do not substitute a one-sided complex phasor
    unless the normalization map is proved.
11. Solve the retarded Maxwell Cauchy problem through the switch and derive the
    late source-free field. Establish its spectral energy norm on the exact
    transverse/constraint space, including any Coulomb endpoint field and the
    zero-frequency row.
12. Determine the function spaces and weighted resolvent estimates that make
    the late-field map continuous. The algebraic Fourier formula alone is not
    the analytic theorem.

**Unit C verdict:** adjudicate exact conservation and the late free-field norm
separately from any carrier or resonance interpretation.

## Unit D — long-time growth, weak switching, and relative flux

13. For a fixed ramp and a long monochromatic plateau, derive the spectral
    autocorrelation representation of the radiated energy. State the exact
    integrable-first-correlation hypothesis and prove
    `E_T=T P+O(1)` with an error uniform in `T`.
14. Derive the weak-switch endpoint theorem and decide in what norm and under
    what spectral regularity it is sharp. Separate vanishing endpoint energy
    from vanishing plateau radiation.
15. Formulate Poynting balance on a world tube moving with the carrier. Track
    the Reynolds-transport term, field energy inside the tube, mechanical
    source work, tube normal/orientation, and the relative flux
    `S-c_g e_z u_EM`.
16. Prove the large-radius and large-time limits in their stated order, with
    sufficient spatial decay or weighted local-energy bounds. A stationary
    shell calculation does not alone justify either interchange.

**Unit D verdict:** adjudicate the exact finite-time identity, asymptotic
`T P+O(1)` theorem, endpoint energy, and comoving flux at their distinct
hypotheses.

## Unit E — amplitude versus energy/action decay

17. Freeze whether the effective variable is a complex amplitude, quadratic
    energy, or KKS-normalized action. If `A(t)` obeys a conditional damping
    law, derive the induced law for `|A|^2` and expose the reciprocal factor of
    two in lifetimes/rates.
18. Check that no amplitude equation, Markov approximation, resonance pole, or
    normalized mode action is inferred from a classical source-radiation
    identity alone.

**Unit E verdict:** preserve the factor-two conversion as an exact conditional
identity, not as an established decay theorem for the target carrier.

## Unit F — supplier and interpretation boundary

19. Identify exactly which inputs remain missing before applying the source
    theorem to a physical Cao/analyzer mode: actual current, its constraint and
    shell trace, physical KKS/action normalization, mixer response, gate time,
    and the outgoing coupled resolvent or resonance construction.
20. Keep switching/repreparation external unless an autonomous mechanism is
    derived. A finite radiated-energy formula supplies no Born rule, reset,
    exchange character, or stochastic detector.
21. State explicitly that the target does not establish action selection,
    stability, an electron or neutrino, P4/P5, or a bare-Euler finite-speed
    theorem.

**Unit F verdict:** preserve the strongest source-to-radiation theorem while
leaving every physical-carrier and particle inference at its named next
construction.

## Final decision contract

The final review will give Units A--F separate route verdicts: established as
stated, refuted with a named mechanism, or blocked with the missing
construction named. It will record source applicability, exact conventions,
API/oracle evidence roles, the strongest combined theorem, and the next
executable dependency. No unit verdict automatically decides a sibling or the
P253 parent campaign.
