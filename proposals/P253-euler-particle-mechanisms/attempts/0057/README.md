# P253/0057: joined independent review of 0052 and 0054

This README preregisters one joined, fixed-boundary, non-author review of the
worker-owned attempts `0052` and `0054`, with separate verdict components.
Unit A reviews the exact corrected resonance, physical KKS normalization, and
fixed-domain/chart claims actually earned by `0052`. Unit B reviews the claimed
`HJ2` moving-boundary whole-space Green/Leray transfer in `0054`. The joined
review then asks whether those two units really close the full graph-resolvent
and Riesz bridge.

The review is designed to preserve exact finite-core inputs even if a global
operator step remains conditional. Conversely, an exact chart identity or
local shape coefficient cannot inherit all-sector complement control merely
because it is adjacent to an exact verifier. No nonlinear Euler branch,
stability, P2 completion, particle, or quantum claim is in scope.

## Independence, ownership, and activation boundary

`particle-balance-review` did not author or implement `0052`, `0054`, their
derivations, corrections, verifiers, or receipts. It owns only append-only
review artifacts under `attempts/0057`. `particle-foundations` owns both target
units and retains responsibility for any scientific correction, source or API
change, central proposal, registry, validation, memory, and commits. The
reviewer will not edit either target, any source/API/test, any central or
generated file, memory, or Git history.

Before this freeze, no `0052` or `0054` README, derivation, result, source
audit, validation, correction, verifier, command, stdout, stderr, or exit body
was opened. Only filenames, content-blind SHA-256 hashes, central ownership
metadata, accepted authority metadata, Git metadata, and the user's specified
review questions were inspected. Substantive source and memory searches are
deferred until activation so the target argument cannot influence the frozen
criteria.

No target body may be opened until root:

1. registers `0057` centrally with this independent owner and exact joined
   scope;
2. records the repository-schema command, stdout, stderr, and exit under
   `attempts/0057`; and
3. makes `attempts/0057/activation-schema.exit` contain exactly `0`.

After activation, the reviewer will recheck every frozen hash, open the target
and declared source bodies, and perform one substantive joined pass while
keeping Unit A and Unit B verdicts separate. At most one evidence-driven
bounded correction package may follow, followed by one check limited to its
changed statements and direct dependency edges. The reviewer will then write
only `0057/review.md` and `0057/verdicts.yaml` and report to `w4:p2`.

The accepted authority base is release `v0.183.0`. At preregistration,
`governance/releases/current.yaml` and `governance/claims.yaml` have SHA-256
values
`7b889a3c0186c0c22e698fd4ccc69e19066b8262018696cfc41e573ec2f0ffd8`
and
`78dbe694359739ee15f7db2b2cb1ccde7f48b5af34f6d9cd451282732dd08983`.
The observed repository head is
`d7357eff06c30fd91cfc68de02f8edabe458005f`. The shared dirty worktree
contains unrelated active work under `0055` and `0058` and will be preserved.
The physics-skill preflight completed with `ok=7`, `warn=0`.

The discovered verifiers are symbolic/operator-form checks, not production
small-ratio numerics. If the activated bodies introduce a numerical spectral
gap, resolvent edge, soft eigenvalue, or shape remainder near a discretization
floor, the small-ratio skill must be loaded before that evidence is
adjudicated. Otherwise no numerical floor prescription binds.

## Frozen Unit A inventory: P253/0052

The complete current `0052` package is frozen content-blind:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0052/README.md` | `219334fcd5bbf40f69d691dc45596a82926b5f25e63836dc3a9eee5990ba6bd1` |
| `0052/derivation.md` | `d353a4ffb8ea1ce65b419f8152b3292add53671a3d5e94195685722ae95f1d5b` |
| `0052/result.yaml` | `9bab0def7992099eace135783d77e4e2a84f96f55792089d0f51177912afc8a8` |
| `0052/source-audit.md` | `1dd942371cf4591ad991dfe9bff8ab67b9a64895ffc84f06442adc04dfd558c2` |
| `0052/validation.md` | `e7a9ee362443ceca42dbb631977ad87218fe6a3353ff9722079007f7810de509` |
| `0052/cross-delta-volume-correction.md` | `339f1ed35a6a2213aa161fa32ba0968c6e22429ff03dd05d1e993f16200e068f` |
| `0052/verify_fixed_domain.py` | `1caf2ac69fdd619b9ebcf01b98ea61d9e97328680f870e52f2c54f12be4d2884` |
| `0052/activation-schema.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0052/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0052/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0052/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0052/fixed-domain.command.txt` | `99d1d6a1b3976fa83b93de87c0b20417027c369da6835e305cbce99844e5856a` |
| `0052/fixed-domain.stdout` | `ea2b14df27d533e1b4810f1e80e9a9fa7e7abf79404955461c51ffda7345c182` |
| `0052/fixed-domain.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0052/fixed-domain.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

## Frozen Unit B inventory: P253/0054

The complete current `0054` package is frozen content-blind:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0054/README.md` | `a13c902a9798e48897aa787f3eecdaf331843a78d0ade783690d8ca476e96a8c` |
| `0054/derivation.md` | `62e3d45cd8c45099c51b47467339a0824b3fdc1ddedcc371645258723575d760` |
| `0054/result.yaml` | `f060d63e92cdf69d414a0a795dcbd35d4e5e17da6ae0557de6d5a7e4d28eebe6` |
| `0054/source-audit.md` | `4ba9323c8f2e16f3418d1809d1245b49e9d688f924f7a50d6b6b952e6c508f6f` |
| `0054/validation.md` | `e5ee164a75b4e75a13a48d2056c9a1f15d867c3a96ae7902881522a9ec88e829` |
| `0054/verify_hj2.py` | `0f04bad95ee0954050490a372f16bea651c0e0f92838d9d6dfa36f26aabae83c` |
| `0054/activation-schema.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0054/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0054/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0054/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0054/hj2.command.txt` | `6863c9a3244091d67bb87b33224ed38ec7d9e5907ea8fc55fe9353a0641d01c8` |
| `0054/hj2.stdout` | `e1250273994c4ce3baa04a9b3308424f0fcd10e651a4d2115efef3e367160e8b` |
| `0054/hj2.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0054/hj2.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

No new canonical Hanzawa/Piola, HJ2, or graph-Riesz API/test pair was found by
filename discovery. The only nearby canonical supplier is the previously
reviewed twisted-carrier pair; it will be treated only at its established
scope if either target actually declares it:

| Possible supplier | Preactivation SHA-256 |
|---|---|
| `src/substrate_framework/euler_twisted_carrier.py` | `f2a2c58235a0f2c8ee358123b383b8e919d22ad095bc0f70b67142c7fdbd5214` |
| `tests/test_euler_twisted_carrier.py` | `f7ae11d6a2bfc292675740e43357d9140257f5a206a582fafb903d03042956cc` |

Package placement and zero exit establish provenance, not the continuum
operator theorem. After activation, each verifier predicate will be mapped to
the exact proposition it proves and the operator steps it cannot prove.

## Fixed joined target and non-inheritance

This is one coherent operator-transfer review with two independently
adjudicated components:

- **Unit A (`0052`)** may establish corrected resonance/KKS/chart or
  fixed-domain finite-core statements without proving `HJ2`, a full
  complement, or a Riesz family.
- **Unit B (`0054`)** may establish a genuine whole-space first/second shape
  derivative or conditional transfer without inheriting `0052`'s finite-core
  normalization, DA quotient, all-sector resolvent, or contour bounds.

The joined result is earned only if the exact outputs of Unit A match the
inputs and norms of Unit B and every remaining complement/contour step is
closed. A local source expansion, formal commutator count, or finite symbolic
matrix cannot silently substitute for a common-domain whole-space operator
estimate.

The target is a fixed theorem/operator claim, not a mechanism selection, and
has no empirical comparator. Exact source versions and cached content hashes
will be recorded from each target's source audit after activation, before the
source theorem is used.

## Check 1: corrected resonance, KKS, and cross-`delta` chart

1. Reconstruct the corrected resonance scale, mode index/remainder, physical
   KKS sign and magnitude, density/Fourier/Jacobian factors, and the exact
   finite-core parameter hierarchy used by `0052`.
2. Identify the reference and deformed domains at two unequal parameter or
   volume values. Verify that the Hanzawa map is orientation preserving,
   globally one-to-one, collar supported as claimed, and has the declared
   determinant rather than being silently volume preserving.
3. Check the Piola transform with its determinant and cofactor factors for
   velocity, vorticity, pressure reconstruction, normal traces, divergence,
   and the energy/KKS pairings. An unequal-volume map cannot be replaced by a
   volume-preserving formula.
4. Construct the common dynamically accessible class: push forward compact
   generators, constraints, circulation/impulse/center slices, stabilizers,
   and the quotient. Verify that changing `delta` does not identify inequivalent
   KKS leaves or discard a harmonic/closed row.
5. State both graph domains and prove two-sided energy/graph norm equivalence
   with uniform constants in the claimed parameter range.

**Maximum Unit A verdict:** the exact corrected resonance/KKS quantities and a
genuine unequal-volume fixed-domain DA equivalence at the stated finite-core
scope. This check alone cannot supply the Green/Leray jet or Riesz projection.

## Check 2: source augmented-map surjectivity and branch identification

6. Pin the actual theorem, proposition, function spaces, symmetry class,
   parameters, gauges, and smallness hypotheses used from every source.
   In particular, treat Cao equation (3.36) as an auxiliary asymptotic
   parameter relation unless the source explicitly proves it is the exact
   steady constraint used by the target.
7. Write the exact augmented nonlinear map `F`, its domain, codomain, Fredholm
   index, derivative, kernel/cokernel, side constraints, parameter rows, and
   parameter Schur complement. Verify actual surjectivity or a bounded right
   inverse in precisely the spaces required by `0054`.
8. Distinguish existence, local uniqueness modulo symmetries, and uniqueness
   within a restricted source ansatz. Determine whether the theorem identifies
   the same centered, fixed-circulation polyhomogeneous branch across `delta`,
   including its logarithmic terms, or merely one solution at each parameter.
9. Check that source uniqueness controls the selected Hanzawa/Piola pullback
   and its first and second parameter derivatives. A unique physical profile
   does not automatically identify a differentiable graph-domain branch.
   Derive the exact restriction relating the source variables `epsilon`, `q`,
   and `tau` before invoking an implicit-function or uniqueness theorem.

**Maximum verdict:** source-licensed local branch and augmented inverse with
the quantifiers needed for the claimed polyhomogeneous transfer.

## Check 3: `p>=6` interface traces and Hardy equivalence

10. Track the exact order at which the compact profile and vorticity vanish at
    the free boundary for `p>=6`. Compute every coefficient or reciprocal
    weight used by the first and second derivatives.
11. Verify interior and exterior trace spaces, jumps, normal derivatives,
    transmission conditions, axis regularity, and collar compatibility after
    two shape derivatives. Include the traces of `A_delta eta`, not only those
    of `eta` or the base field. No interface distribution may be hidden by a
    flat extension.
12. Prove both directions of every weighted Hardy estimate and norm
    equivalence, with constants uniform over the allowed deformations. Prove
    that the transformed range is the common DA range in both directions and
    include the composite homogeneous/zero-moment Hodge bounds used to recover
    velocity and pressure. Check endpoints and whether derivative loss
    requires more than `p>=6`.
13. Match those spaces to the actual Euler graph domain and its adjoint, not
    only to scalar profile regularity.

**Maximum verdict:** the precise `p>=6` trace/Hardy license for the claimed
two-jet and common graph norm.

## Check 4: whole-space Biot--Savart/Leray shape derivatives

14. Differentiate the pulled-back whole-space Biot--Savart/Green kernel,
    including determinant, cofactor, evaluation-point, source-point, cutoff,
    gauge, and far-field terms. State the operator spaces and remainder norm.
15. Derive the first and second Leray derivatives either from the exact
    Biot--Savart curl representation or from the full pressure Poisson problem.
    Verify sign, tensor contraction, boundary/interface contribution, and
    commutators with the Piola transform.
16. Retain the velocity and pressure tails outside compact vorticity support.
    A local meridional source operator or compact-core coefficient does not
    prove a whole-space shape derivative.
17. Prove differentiability in operator/graph norm, not only pointwise kernel
    differentiation, and give the uniform small-`delta` bounds used downstream.
    A matched scalar kernel integral does not by itself establish an
    `H^s -> H^(s+1)` symbol remainder; prove the differentiated multiplier or
    pseudodifferential estimate with the stated uniformity.

**Maximum Unit B verdict:** exact first and second whole-space Green/Leray
shape derivatives on one common graph domain with the stated remainder.

## Check 5: front/side matching and global remainder

18. Identify every front, side, collar, axis, interface, and far-exterior
    region, its local coordinates, and the overlap domains. Verify that the
    expansions agree to the claimed order in each overlap.
19. Check cutoff commutators and partition-of-unity derivatives, including
    logarithmic scale factors and the order in which the core radius, collar,
    and exterior scale are taken.
20. Prove the assembled remainder in the declared global energy/graph norm.
    Pointwise front and side errors with incompatible weights do not imply one
    global `o(delta^2)` estimate.

**Maximum verdict:** a uniform matched whole-space remainder with no uncovered
region or double-counted logarithmic term.

## Check 6: three commutators and two-index summability

21. Identify the exact generator, basis, weights, and two mode indices in the
    claimed matrix coefficients. Write the boundary and interface terms for
    each of the three integrations by parts or commutators.
22. Verify that the operator and basis vectors remain in the required domains
    after every commutator. Track derivative loss, `m`, radial-index, `delta`,
    collar, and exterior constants.
23. Derive the actual off-diagonal decay. If the target's equation (48) uses
    one generator `Lambda`, three commutators with that single `Lambda` control
    only the corresponding combined spectral difference; they do not by
    themselves yield a product of independent decay factors in `m` and `n`.
    Require separate index generators, mixed commutators, or one directly
    proved two-dimensional summable difference bound, with uniform eigenbasis
    constants. One-index decay or an asserted product bound is not
    automatically summable over both indices.
24. Check the diagonal and near-diagonal bands separately and prove the exact
    Schur/Hilbert--Schmidt/trace-class conclusion used later.

**Maximum verdict:** the declared uniform two-index summability and the exact
operator ideal it earns, not merely formal high-mode decay.

## Check 7: all-sector Grushin/complement closure

25. Define the finite block, injection/projection, complement, Grushin
    operator, graph norm, spectral set, and symmetry quotient. Check physical
    KKS biorthogonality and neutral/stabilizer rows.
26. Cover every azimuthal and radial sector, including `m=0`, high `|m|`,
    critical layers, the axis, collar/interface traces, and whole-space
    exterior pressure. State uniform constants and the allowed dependence on
    `delta` and contour size.
27. Prove kernel and cokernel control, adjoint-domain compatibility, Fredholm
    index, and invertibility of the complement block. A finite list of modes or
    sequence-level Schur estimate cannot close the full complement.
28. Verify the effective finite matrix and its reconstruction error only after
    the complement inverse is earned. For the nonnormal all-sector operator,
    reconstruct the actual `X^* -> D^*` adjoint resolvent estimate; black-box
    spectral separation or a distance-to-spectrum bound is not a resolvent
    norm estimate.

**Maximum verdict:** a full-sector graph-domain Grushin inverse and complement
closure on the common DA quotient.

## Check 8: contour-relative resolvent and Riesz norm

29. State the reference and deformed operators on the same space and domain,
    the contour, its distance to each spectrum, and all parameter quantifiers.
30. Derive the relative resolvent identity with the perturbation in the actual
    graph-to-energy norm. Check whether the estimate is absolute, relative to
    the reference resolvent, or weighted by the spectral gap.
31. Integrate the correct contour bound to obtain the Riesz projection norm,
    including contour length, gap powers, `delta` dependence, and the rank
    stability argument. A pointwise finite matrix residual is insufficient.
32. Distinguish bounded projection convergence, graph-norm convergence,
    similarity of invariant subspaces, and an approximate spectral window.
    Record exactly which one is proved.

**Maximum joined verdict:** the claimed contour-relative norm-resolvent and
Riesz theorem with explicit common-domain constants. It does not imply a
nonlinear branch or all-time nonlinear stability.

## Exact finite-core atoms and conditional bridges

The review will separately preserve any independently supported:

- corrected resonance location and remainder;
- physical KKS sign, density, volume, and real/complex normalization;
- local Cao/profile cells and polyhomogeneous coefficients;
- exact Hanzawa/Piola shape identities for a declared finite deformation;
- exact toroidal spectral-distance identities and leading thin-core scaling;
- first or second kernel derivative at a local or weak operator scope;
- front/side asymptotic coefficient;
- one-index or sequence-level commutator/Schur bound; and
- conditional finite Grushin or Riesz implication.

None will be discarded merely because a later global bridge is missing. Each
will be labeled exact, source-conditional, formal, or regression-supported at
its actual domain. The two unit verdicts will not inherit evidence from one
another without an explicit matching map.

## Source, oracle, and correction protocol

After activation, every primary source named by the two source audits will be
pinned by exact version and content hash before use. The review will inspect
the actual source theorem and proof ingredients for domain, symmetry,
parameter, regularity, uniqueness, and operator applicability. A source name
or nearby elliptic theorem is not a substitute for the required whole-space
Euler map.

The attempt-local verifiers will be audited only at their actual predicates.
The review will inventory assertion sites, runtime tally, mutations or
counterexamples, and exact quantities, but will not treat symbolic finite
algebra as proof of an infinite-dimensional graph theorem. Existing receipts
will be reused while code and inputs remain frozen; no unchanged oracle will
be rerun merely for a pass count.

If one concrete defect is found, the review will state the strongest supported
result, exact unsupported extension, minimum correction, and evidence that
would restore the stronger result. One bounded target correction package may
then receive one changed-statement/dependency check. No second substantive
pass will be performed.

## Verdict map and exclusions

The final record will contain separate Unit A and Unit B verdicts plus one
joined operator-transfer verdict. Each of the eight checks above receives one
route-scoped status: `established`, `refuted` with its mechanism, or `blocked`
with its missing construction. Conditional implications are recorded as such
and cannot be promoted by a passing test tally.

The frozen review excludes nonlinear Euler branch existence, convergence of a
formal branch, nonlinear or orbital stability, P2/P4/LP4 completion or no-go,
all particle or quantum interpretations, action or detector selection,
registry promotion, active `0055`/`0058`, source or central edits, and commits.
The parent campaign remains active regardless of the joined review outcome.
