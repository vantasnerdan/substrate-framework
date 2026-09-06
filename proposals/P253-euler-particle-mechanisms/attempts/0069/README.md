# P253/0069: joined independent review of 0064 and 0065

This README preregisters one joined fixed-boundary, non-author review of the
root-authored target attempts `0064` and `0065`, with a separate verdict unit
for each target. Unit A reviews `0064`'s exact compact-vorticity
multipole/pressure/interaction identities and local isotropic
scalar-to-transverse obstruction. Unit B reviews `0065`'s puncture-flux energy
bound, loop circulation, helicity, and nonlocal-multiplier constructions. The
joined transaction asks what exact charge-like Euler structures these two
routes genuinely supply; neither unit inherits the other unit's evidence or
verdict.

The strongest useful exact statement will be preserved wherever its equations,
domain, regularity, parity, and finite-energy hypotheses close. A route-scoped
obstruction will not be promoted into a global charge no-go, P5 completion, or
a particle claim.

## Independence, ownership, and activation boundary

`particle-balance-review` authored or implemented neither `0064`, `0065`,
their derivations, APIs, tests, verifiers, strategy/scope receipts, nor evidence
receipts. It owns only append-only review artifacts under `attempts/0069`.
The reviewer will not edit a target, source/API/test, central or generated
record, memory file, or Git history.

Before this freeze, no `0064` or `0065` README, derivation, result,
source-audit, validation, verifier, output, hash manifest, strategy/scope
receipt, API, or test body was opened. Inspection was limited to filenames,
content-blind SHA-256 hashes, ownership and Git metadata, accepted-authority
metadata, and the user's frozen checks. In particular, the supervisor's
`0064/sol-high-recoil-scope.md` and the `0065` pre-review
`sol-high-linear-scope-correction.md` are frozen below but have not been read.

The review opens only after the coordinator at `w4:p2`:

1. centrally registers `0069` with this exact joined scope and reviewer;
2. records the repository-schema command, stdout, stderr, and exit under
   `attempts/0069`; and
3. makes `attempts/0069/activation-schema.exit` contain exactly `0`.

After activation, the reviewer will recheck all frozen hashes, open each source
audit and the two supplied scope/strategy records first, pin the exact primary
source versions and propositions used, and only then open the remaining target
bodies, API, tests, and receipts. One substantive joined pass may request at
most one evidence-driven bounded correction package, followed by one
correction check limited to changed statements and directly affected edges.
Final writes are limited to `0069/review.md` and `0069/verdicts.yaml`, followed
by a report to `w4:p2`.

The accepted authority base is release `v0.183.0`.
`governance/releases/current.yaml` and `governance/claims.yaml` have SHA-256
values
`7b889a3c0186c0c22e698fd4ccc69e19066b8262018696cfc41e573ec2f0ffd8`
and
`78dbe694359739ee15f7db2b2cb1ccde7f48b5af34f6d9cd451282732dd08983`.
The observed repository head is
`78b12e7fc64db7ebad7631a4932d96de880b74a1`. The shared dirty worktree
contains accepted 0063 closure, active 0066--0068/0070 work, and newer gauge
APIs/tests. Those surfaces belong to other agents or excluded attempts and will
be preserved. The physics-skill preflight completed with `ok=7`, `warn=0`.

The target verifiers and focused test are exact symbolic/algebraic evidence;
the user has excluded production numerics. No small-ratio numerical
prescription binds this review.

## Frozen Unit A inventory: P253/0064

The complete current `0064` package is frozen content-blind:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0064/README.md` | `84b86e0eae440081f7db84c36e6e6db799e091a8645c1b78c2032cb8aaad4808` |
| `0064/derivation.md` | `39afd045652ca89b3e888e5394857f9d54969a4493b81eb69d781e4be0ff9fd7` |
| `0064/result.yaml` | `b205029a3e1d28e8cc646c3b22906f197ceaccc03cf444bc4c4d750100082fb6` |
| `0064/source-audit.md` | `de171b82b976fc991883608805fe4e8df45232232480e8e4d313894c55d81309` |
| `0064/validation.md` | `20c07764f38f821bd265aa76aeb031cfaa74d7d0c620f3bd5a39b7ba7f7fbff5` |
| `0064/sol-high-recoil-scope.md` | `d3382d11402fdedbd0c1510164443605582a2a6aff7c137740ea932fe63cb7e4` |
| `0064/verify_charge_multipoles.py` | `1f7b6cd9b36aa731d07de1094bff8420c1d06ec980152447fb0ab70a63e2fb55` |
| `0064/artifact-hashes.sha256` | `79aa8d6e641e004a7a190d3b7e7a509f738a09f6197e06b038609596c333d91f` |
| `0064/activation-schema.command.txt` | `0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8` |
| `0064/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0064/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0064/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0064/exact-check.command.txt` | `929838d2fe11cc7842cec61795dcb6c2033227a01286f9221007d06248337d42` |
| `0064/exact-check.stdout` | `eef3e2054a3ffa0a52f4a546a5ad7163499a185bb9e04634da757027e15fb6c6` |
| `0064/exact-check.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0064/exact-check.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0064/focused-pytest.command.txt` | `07bebd6f400c6a8f8c3a3969355ff430dcff59d7b9fe785cd0c4c05b30256203` |
| `0064/focused-pytest.stdout` | `f50ba8b7f6a3d95a43d580422e49c56f7a02f2710de3a321bd99fe4eb177642e` |
| `0064/focused-pytest.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0064/focused-pytest.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

## Frozen Unit B inventory: P253/0065

The complete current `0065` package is frozen content-blind:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0065/README.md` | `78d178714b88c5ac4e55e6a917fbb9207b8d8013abc864c478d00e11b569dfa7` |
| `0065/derivation.md` | `2b1a34608665741b41a0add3d0a941ce68d5b0df8bee8a0f14a886a1d97555ab` |
| `0065/result.yaml` | `d86bd7589121f82d67b24635764bffef906f756ecd1d8b4851d4f8d130703807` |
| `0065/source-audit.md` | `8a777460dd0c7b568a06e069711d95959bbdb1ade7118da2d4604587076e6fd8` |
| `0065/validation.md` | `f02ea9aff67da31acf3ad1faf0a8cef429e683b23021437f7dba71be35e7588e` |
| `0065/sol-high-linear-scope-correction.md` | `d025fd5c1156cc4ed0d5921fc787c3a0cbb2eaf49448b4272a3a18f999b2e757` |
| `0065/verify_topological_charge.py` | `9396c82dae6867d3c7cd98e2d0cb31a00bd5c8553f05ee3aab2953672e4daefa` |
| `0065/artifact-hashes.sha256` | `01cb76225f015a71ffd479a44f8e657433c5abc287db9eb3b49bfee15965f30b` |
| `0065/activation-schema.command.txt` | `0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8` |
| `0065/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0065/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0065/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0065/checkpoint-validation.command.txt` | `865ed1b4f098f8c9b967918f5b77f563acd4929fde0d4532dfe2110bb816fd6c` |
| `0065/checkpoint-validation.stdout` | `3f0e009abc09810e3af6bb670217ccfcd1f10532edf46787606444045551ed0f` |
| `0065/checkpoint-validation.stderr` | `bdfcc6ae2ff88fb63bef3f9aca46020788794fa092c3334a343bcc0c466a68d1` |
| `0065/checkpoint-validation.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0065/exact-check.command.txt` | `7b4a22d7873c0f5a99e363fced922e70d54020317e809cb4ca860869b4603c66` |
| `0065/exact-check.stdout` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0065/exact-check.stderr` | `1834fa5b64bb568a96a02efdd8eba66b8664a7a03c552743f968864941dc1fb3` |
| `0065/exact-check.exit` | `4355a46b19d348dc2f57c046f8ef63d4538ebb936000f3c9ee954a27460dd865` |
| `0065/exact-check.repaired.command.txt` | `7b4a22d7873c0f5a99e363fced922e70d54020317e809cb4ca860869b4603c66` |
| `0065/exact-check.repaired.stdout` | `a7f36cb8636e4899c2be8f4d83a70189c42459fb1120c0d1d8d9d3244c82f670` |
| `0065/exact-check.repaired.stderr` | `9c4918c65c1245adba03c1941e04e4953107c6e6f96d29a6720ac0e6a2a2677f` |
| `0065/exact-check.repaired.exit` | `4355a46b19d348dc2f57c046f8ef63d4538ebb936000f3c9ee954a27460dd865` |
| `0065/exact-check.final.command.txt` | `7b4a22d7873c0f5a99e363fced922e70d54020317e809cb4ca860869b4603c66` |
| `0065/exact-check.final.stdout` | `3b32a6c42bb21654148e29bac28e24faaa36d5b39d551c54fa6455be8ae761cb` |
| `0065/exact-check.final.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0065/exact-check.final.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

## Frozen importable API and focused test

The target commit names one importable API/test pair for this transaction:

| Artifact | Preactivation SHA-256 |
|---|---|
| `src/substrate_framework/euler_charge_multipole.py` | `42be8f6699584984b53293b5e6283d30d3f6879ce300136812dd7cc5e0fcc86c` |
| `tests/test_euler_charge_multipole.py` | `52638a95caa79356eb7b09b1b08d62f3d7cb29321bc09e5ce580282232cfce91` |

The untracked `euler_gauge_current.py`, `euler_gauge_emergence.py`, and their
tests are later 0067/0068 work. They are expressly excluded and will not be
opened or cited. Package placement, a check tally, and an exit-zero receipt are
provenance or regression evidence, not scientific authority.

## Fixed joined target and non-inheritance

The transaction has two independent verdict units:

- **Unit A (`0064`)** may establish exact compact-vorticity multipole,
  pressure, and initial pair-interaction identities, or a local
  scalar-to-transverse obstruction, without establishing a localized Euler
  charge or a universal no-go.
- **Unit B (`0065`)** may establish a punctured-domain energy obstruction, a
  topological circulation statement, an exact helicity identity, or a
  conditional nonlocal transverse construction without inheriting Unit A's
  compact-support, interaction, or regularity licenses.

The supervisor `0064` strategy audit and the `0065` pre-review scope receipt
are exposing instructions and claim-boundary records. They do not predetermine
either verdict and cannot substitute for the equations or cited sources.

## Unit A check 1: compact-vorticity multipoles and pressure

1. Freeze the Euler convention actually used: physical versus kinematic
   pressure, density `rho_0`, Green-function normalization, Levi-Civita and
   Fourier signs, and the definitions of velocity, vorticity, impulse, and all
   multipole tensors.
2. Reconstruct the Biot--Savart and pressure-Poisson equations directly from
   incompressible Euler. Check every index contraction, trace subtraction,
   symmetry, sign, density factor, `4*pi` coefficient, integration by parts,
   and surface term.
3. Derive which compact-vorticity moments vanish from divergence/curl
   structure rather than from informal isotropy. Include translated moments,
   parity, origin dependence, tensor symmetries, decay needed for convergence,
   and the distinction between compact vorticity and compact velocity.
4. Derive the exterior velocity and pressure asymptotics with a controlled
   remainder on the actual domain. State whether identities are exact,
   asymptotic, initial-time only, or propagated by Euler; no category may be
   upgraded by prose.

**Maximum verdict:** the exact compact-vorticity multipole and pressure theorem
with its actual regularity, support, parity, and finite-energy scope.

## Unit A check 2: pair interaction

5. Reconstruct the cross-energy, pressure source, force/centroid derivative,
   and any quadrupole or higher-moment interaction for the precise two-field
   configuration. Check displacement conventions, action--reaction signs,
   density, normalization, convergence, and whether disjointness holds only at
   the initial time.
6. Separate velocity tails, pressure tails, material acceleration, centroid
   acceleration, and a mechanical force. Determine which observable is finite
   and which conservation law licenses its evolution.
7. Test special orientations, parity reversals, exchanged carriers, and a
   vanishing leading tensor as exact sensitivity cases. A far-field power law
   does not become a Coulomb law or a point-particle equation.

**Maximum verdict:** the strongest exact or controlled-asymptotic interaction
identity at the time interval and observable actually established.

## Unit A check 3: local isotropic scalar-to-transverse obstruction

8. State the class of maps under review: local finite-order differential,
   translation/rotation invariant, linear or nonlinear, scalar input versus
   pseudoscalar input, and the precise transverse output space. Derive its
   Fourier/principal symbol and representation content rather than relying on
   a slogan about isotropy.
9. Check zero frequency, singular coefficients, parity-odd epsilon tensors,
   background vectors/tensors, multiple scalars, nonlinear gradients, and
   nonlocal projectors. Identify exactly which added structure escapes the
   obstruction.
10. Keep the conclusion local and class-scoped. It cannot exclude nonlocal
    Hodge multipliers, oriented sources, topology, boundaries, anisotropic
    backgrounds, or every Euler environment.

**Maximum verdict:** an exact no-map theorem for the explicitly quantified
local isotropic scalar construction class, with named escape routes.

## Unit B check 1: puncture flux and finite-energy lower bound

11. State the punctured domain, deleted set, number of ends, regularity, trace
    class, divergence equation, and whether the claimed flux is velocity,
    vorticity, or another current. Verify flux independence across homologous
    surfaces using the correct divergence theorem.
12. Reconstruct the spherical Cauchy--Schwarz/coarea lower bound with the exact
    surface area, radial measure, density, and energy factor. Determine at
    which endpoint it diverges and whether the conclusion concerns local or
    total finite energy.
13. Audit distributional sources at the puncture, admissible boundary data,
    harmonic Hodge rows, multiple punctures and compensating flux. Smooth
    source-free Euler on all of `R^3` cannot silently inherit a punctured-domain
    flux sector.

**Maximum verdict:** the exact flux-energy theorem on its declared punctured
domain and the corresponding finite-energy incompatibility, not a global
charge no-go.

## Unit B check 2: loop circulation

14. Define the loop, linking class, spanning surface, field regularity, and
    domain topology. Check when Stokes' theorem applies and when no spanning
    surface exists inside the domain.
15. Reconstruct circulation orientation, signs, and physical units. Separate a
    fixed-loop topological period from Kelvin circulation of a material loop;
    state the Euler smoothness, pressure/barotropic, boundary, and time-interval
    hypotheses needed for conservation.
16. Determine whether circulation labels a cohomology class, a carrier orbit,
    or only chosen initial data. Nonzero circulation is not by itself a local
    scalar charge, pair force, or particle observable.

**Maximum verdict:** the exact loop/cohomology and Kelvin statements with their
distinct domains and evolution hypotheses.

## Unit B check 3: helicity and nonlocal multiplier

17. Derive the helicity integral and any local/nonlocal density decomposition
    with boundary terms, gauge/Hodge choices, density convention, decay,
    integrability, and conservation hypotheses. Check whether compact
    vorticity suffices for each claimed integral and whether velocity tails
    alter convergence.
18. Audit parity under proper and improper rotations. Distinguish scalar from
    pseudoscalar, polar from axial vectors, and orientation-dependent from
    isotropic constructions.
19. For each nonlocal multiplier, write the exact Fourier symbol, order,
    reality/Hermitian symmetry, transverse projector, kernel, zero-frequency
    prescription, and mapping domain. Check infrared and ultraviolet behavior
    in the claimed finite-energy/Sobolev class and any singular support.
20. Separate algebraic transversality from an Euler solution, conservation
    law, local mediator, or compactly supported construction. A nonlocal symbol
    may supply a conditional field map without supplying charge selection,
    interaction dynamics, or P5 completion.

**Maximum verdict:** exact helicity identities and the strongest well-defined
nonlocal transverse map on its declared functional domain.

## Source applicability and oracle audit

21. Open the two frozen source audits first after activation. For every cited
    article, theorem, textbook identity, or repository dependency, pin its
    version/path and match its hypotheses, domain, regularity, orientation,
    pressure convention, and conclusion to the exact proposition imported.
    A source name or analogous formula does not transfer a theorem.
22. Map each symbolic verifier and API/test predicate through the full bridge:
    computed expression -> exact proposition -> object/domain licenses ->
    maximum route verdict -> effect on P253. Inspect unevaluated symbolic
    objects, copied literals, hard-coded expectations, branch/zero-mode
    assumptions, and wrong-sign or wrong-normalization sensitivity.
23. Reuse unchanged captured receipts only after their source, command, inputs,
    and exit/tally agree at the frozen boundary. No production numerical run is
    authorized or required. A focused test supports only the API algebra it
    executes.
24. Apply finite-energy, regularity, parity, and domain checks to every row,
    not only the headline formulas. Explicit hypotheses are legitimate scope;
    absent hypotheses or category changes are support gaps, not automatic
    falsifiers.

## Verdict discipline, exclusions, and deliverables

Each unit receives its own `established`, `refuted`, or `blocked` route verdict
and four-axis verification/review/compatibility/epistemic classification.
`Refuted` requires an explicit counterexample under the stated hypotheses;
missing support yields the strongest established substatement plus the one
missing construction. Any minimum correction states the exact unsupported
extension, the smallest truthful repair, and the evidence that would restore
the stronger result.

Excluded from this transaction are:

- attempts `0067` and `0068`, including their bodies, APIs, tests, and any
  post-freeze conclusions;
- P5 completion, a localized-electron or particle conclusion, and any quantum
  interpretation;
- any global charge no-go or exhaustion inference from these finite route
  checks;
- promotion, central/source/API/test edits, production numerics, commits, and
  downstream work not directly affected by the frozen statements.

After activation, the deliverables are `0069/review.md` and
`0069/verdicts.yaml`. They will preserve the strongest supported result for
each unit, distinguish exact proof from API regression and conditional routes,
name the remaining dependency without broadening it, and report final artifact
hashes and verdicts to `w4:p2`.
