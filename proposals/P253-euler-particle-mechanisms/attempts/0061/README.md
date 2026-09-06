# P253/0061: joined independent review of 0055 and 0060

This README preregisters one joined fixed-boundary, non-author review of the
root-owned attempts `0055` and `0060`, with separate verdict units. Unit A
reviews Euler similarity/action scaling, the compact-swirl pressure multipole,
the classical hyperbolic extension, and lattice convergence in `0055`. Unit B
reviews the exact disjoint stationary composite, proper-octahedral orbit,
vorticity moments, corrected Hodge multipole order, DA product structure, and
lattice first-moment threshold in `0060`.

The review rewards the strongest exact classical theorem in each unit. A
topological orbit may survive while a continuous action rescales; a stationary
octahedral orbit may exist without providing a positive invariant mode. These
distinctions will be adjudicated independently rather than collapsed into one
negative or particle-facing verdict.

## Independence, ownership, and activation boundary

`particle-balance-review` authored or implemented neither `0055`, `0060`,
their formulas, APIs, tests, verifiers, repairs, nor receipts. It owns only the
append-only review artifacts under `attempts/0061`. Root retains responsibility
for target corrections, source/API/test changes, central records, memory, and
commits. The reviewer will write no target, source, central, generated, memory,
or Git-history file.

Before this freeze, no `0055` or `0060` README, construction, result,
source-audit, validation, verifier, output, correction, API, or test body was
opened. Inspection was limited to filenames, content-blind SHA-256 hashes,
ownership and Git metadata, accepted-authority metadata, and the user's frozen
questions.

The source-first activation sequence is fixed as follows:

1. root centrally registers `0061` with this joined scope and independent
   reviewer;
2. root records the exact repository-schema command, stdout, stderr, and exit
   under `attempts/0061`, with exit exactly `0`;
3. after activation, the reviewer first opens the two frozen source audits and
   pins every cited primary source version and content hash;
4. only then does the reviewer open the construction/result/validation,
   verifier, API, test, and receipt bodies.

This ordering freezes source applicability before the target arguments can
influence it. One substantive joined pass may request at most one bounded
correction package, followed by one correction check limited to changed
statements and direct dependency edges. Final writes are limited to
`0061/review.md` and `0061/verdicts.yaml`, followed by a report to `w4:p2`.

The accepted authority base is release `v0.183.0`. At preregistration,
`governance/releases/current.yaml` and `governance/claims.yaml` have SHA-256
values
`7b889a3c0186c0c22e698fd4ccc69e19066b8262018696cfc41e573ec2f0ffd8`
and
`78dbe694359739ee15f7db2b2cb1ccde7f48b5af34f6d9cd451282732dd08983`.
The observed repository head is
`bd4209612ebde3869bdb23a822a18baa4c91bffd`. Unrelated untracked `0057`
checkpoint-validation receipts are present and will be preserved. The physics
skill preflight completed with `ok=7`, `warn=0`.

The discovered checks are exact symbolic/API predicates, not production
small-ratio numerics. No small-ratio numerical prescription binds at this
freeze. If activated bodies unexpectedly use a soft numerical spectrum,
near-floor force, or near-floor energy difference as evidence, that evidence
will be typed separately before it is adjudicated.

## Frozen Unit A inventory: P253/0055

The complete current `0055` package is frozen content-blind:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0055/README.md` | `b39546f5f3f34bccf1bae2aa6b3993886ff37a6d9085e0db6585a5233018daf1` |
| `0055/construction.md` | `372af8de9c86538155f34e7cff597e5bc523b42370caad00eacb75f177f49c3c` |
| `0055/result.yaml` | `ac67b78a4cba97b9bab9fe40da1f3832365cd0765a556130dac92312df0ea81e` |
| `0055/source-audit.md` | `45d69883d8fa7b7e52901a0f760746ed5d56ff83411ac41d749d27523440600d` |
| `0055/validation.md` | `ae9c041f4f6684b30cdddb34b7c62cac9fc9b429e5eb925a86ca79d96ff69989` |
| `0055/verify_scale_causality.py` | `51c0d4ae72c49d4eba525f9d55ac1a535d998dac7d34878aaabf0d976212f5b0` |
| `0055/activation-schema.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0055/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0055/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0055/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0055/exact-check.command.txt` | `08929e42ee9d70dfd37b8dd0840e1c8e37d76c063600fb5ddc8402b7cc119a68` |
| `0055/exact-check.stdout` | `c0b917e062f32d5e9bbaf92204c83b9c3ca4ae730c566266bfa33403d187ef4f` |
| `0055/exact-check.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0055/exact-check.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0055/extension-ledger.command.txt` | `08929e42ee9d70dfd37b8dd0840e1c8e37d76c063600fb5ddc8402b7cc119a68` |
| `0055/extension-ledger.stdout` | `ea56f577a3df62732bcb84c6e363a91d3f890b61851084313290a218916ccbe2` |
| `0055/extension-ledger.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0055/extension-ledger.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0055/focused-pytest.command.txt` | `38b521d1230ab7f164d0d586f0139bb64b1568969fbaa530e8f0e1db7b4fb2d3` |
| `0055/focused-pytest.stdout` | `e6eae79f2ecf62b9d87a811f2f0ee6d9af22657a9582b085e2318955b72f2e3c` |
| `0055/focused-pytest.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0055/focused-pytest.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0055/checkpoint-validation.command.txt` | `ed7bf1f44573eb825083eb58251f6d54dcd7133ba93abad4d213ba18f55d245b` |
| `0055/checkpoint-validation.stdout` | `cb1af9f059e33644ec7e9e61d161d9f6236a5aca551e667daa613c3898458993` |
| `0055/checkpoint-validation.stderr` | `bdfcc6ae2ff88fb63bef3f9aca46020788794fa092c3334a343bcc0c466a68d1` |
| `0055/checkpoint-validation.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

The direct importable pair is frozen separately:

| Artifact | Preactivation SHA-256 |
|---|---|
| `src/substrate_framework/euler_scale_causality.py` | `b15f59edfc3d951d031d4927e50e210e321d6470ccb419c21158c6ec02e280c8` |
| `tests/test_euler_scale_causality.py` | `afa4a16fbaa9f6da505efd126519e220170e985b2813ec08a43e1ff9369597f8` |

## Frozen Unit B inventory: P253/0060

The complete current `0060` package, including the repaired multipole and test
receipts, is frozen content-blind:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0060/README.md` | `74a72c03f156f931888bb7d6b60b13098b7e703a51b2e53aa4a7948f75413475` |
| `0060/construction.md` | `1e5d2c769daf81ea688d1a029fa2f2c46ed5f98b5ad57559b10d1a1062cda372` |
| `0060/result.yaml` | `46f8a0da69e30a01bc46cee0d946f0fe3ad3411deba00eefd84d41a57bf459f7` |
| `0060/source-audit.md` | `d7154821f6953cf0b7e661dd78c90abd4f84426a0d762bfcac71b267d15244f4` |
| `0060/validation.md` | `06d3487954e6af9139b5f8a27423170f413821757630f9cd5a2bc02c4feb84e8` |
| `0060/verify_neutral_cell.py` | `895cbd35955209b8b56bab3560b59809affec613190ff26b8ca3cf4c49369c0f` |
| `0060/multipole-order-correction.md` | `9b3589b0062f502d4ad91c60a0296e7330309e1a21441ec288fcdbb1dd7f181d` |
| `0060/activation-schema.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0060/activation-schema.stdout` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |
| `0060/activation-schema.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0060/activation-schema.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0060/exact-check.command.txt` | `bf4b1fdb53ebc627f10b8797d255d6b9aa2a2f4f55767342825c6047f86f1233` |
| `0060/exact-check.stdout` | `eeb32162de01dbd947deb12e83c177b7e1a2764e174e9dee0aa7db0e2168bfbe` |
| `0060/exact-check.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0060/exact-check.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0060/multipole-order-correction.command.txt` | `bf4b1fdb53ebc627f10b8797d255d6b9aa2a2f4f55767342825c6047f86f1233` |
| `0060/multipole-order-correction.stdout` | `f3a908db95d4fd54c845b0f58223111d596199df74b0f2f568a173f0379e0b07` |
| `0060/multipole-order-correction.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0060/multipole-order-correction.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0060/focused-pytest.command.txt` | `f680d7558c8cd5f6550eafc7430e1feb65a941a2707015ea12dfd9c1b17c59dc` |
| `0060/focused-pytest.stdout` | `b4a5d3759e611cd9dde20a8804ee3e75fd0949886be9c311379bcfc57e88e536` |
| `0060/focused-pytest.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0060/focused-pytest.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0060/focused-pytest-repaired.command.txt` | `f680d7558c8cd5f6550eafc7430e1feb65a941a2707015ea12dfd9c1b17c59dc` |
| `0060/focused-pytest-repaired.stdout` | `50b027a51efd2571ad9c171d98a5ccef20afaedbdcce9b2ed22bed18d22b292e` |
| `0060/focused-pytest-repaired.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0060/focused-pytest-repaired.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

The repaired direct importable pair is frozen separately:

| Artifact | Preactivation SHA-256 |
|---|---|
| `src/substrate_framework/euler_neutral_cell.py` | `b666df746bf7bb96a04fa7abfa3e08b230622d2c8b89cc6bf5d59fb843f4fe63` |
| `tests/test_euler_neutral_cell.py` | `a3b59c32618fa2e5847a3df0aa7bd81ccec258fed0cf232fbe66c24c880da56a` |

Package placement and a zero exit establish provenance, not the continuum
theorem. After activation, each assertion and API output will be mapped to the
mathematical proposition it actually checks.

## Frozen source inventory protocol

The target source inventories themselves are frozen at hashes
`45d69883d8fa7b7e52901a0f760746ed5d56ff83411ac41d749d27523440600d`
for `0055/source-audit.md` and
`d7154821f6953cf0b7e661dd78c90abd4f84426a0d762bfcac71b267d15244f4`
for `0060/source-audit.md`. No attempt-specific external PDF was discoverable
by filename alone. After activation these two files will be opened first; each
cited paper, accepted supplier, local theorem, and canonical module will be
pinned by exact version/path and SHA-256 before either construction body is
opened. A source name, symmetry label, or familiar multipole formula will not
substitute for its actual hypotheses.

## Fixed joined target and non-inheritance

The transaction has two independent verdict units:

- **Unit A (`0055`)** may establish exact similarity weights, a pressure-tail
  multipole, a classical hyperbolic extension, or convergence thresholds
  without selecting an action quantum, a finite-speed cone, or a particle.
- **Unit B (`0060`)** may establish an exact stationary composite and its
  proper-octahedral orbit or far-field cancellation without constructing a
  positive invariant normal mode, Bloch band, or lattice dynamics.

The joined review asks whether these classical scale and neutral-cell atoms
are mutually consistent. Neither unit inherits the other's dynamical or
topological conclusion. A compact stationary product does not make the 0055
hyperbolic coordinate physical, and a convergent lattice sum does not turn a
single-cell asymptotic into an autonomous propagation law.

## Unit A check 1: exact Euler similarity weights

1. State the exact active/passive similarity convention, time scaling,
   velocity, vorticity, pressure, density, support, and integration-measure
   transformations in three spatial dimensions.
2. Recompute the weights of circulation, kinetic energy, helicity, axial
   angular momentum/impulse, and the angular KKS action directly from their
   defining integrals and forms. Distinguish a KKS period from the numerical
   action carried by one orbit.
3. Check every physical dimension and every factor of density and `2*pi`.
   Track whether the dilation changes the Euler solution, the background, the
   coadjoint leaf, or only coordinates.
4. Verify composition and inverse laws for positive scale and time parameters,
   including sign/orientation restrictions and the natural smooth finite-energy
   class.

**Maximum verdict:** exact classical similarity and observable weights in the
declared Euler class.

## Unit A check 2: topology versus continuous action scaling

5. Identify the exact configuration/orbit space and homotopy or exchange class
   claimed invariant under dilation. Prove that positive scaling gives a
   continuous path or homeomorphism on that space without crossing forbidden
   configurations.
6. Separately compute how the physical KKS/action integral varies. Topological
   invariance alone cannot quantize a continuously rescalable action, and a
   continuously varying action does not erase a genuine discrete topology.
7. Check whether circulation or another conserved quantity has been fixed. If
   fixing it removes a scale freedom, state the remaining similarity family
   and whether it still changes the action.
8. Preserve any exact no-cone counterexample at its actual quantifiers. Euler
   similarity and continuous action scaling can refute universal scale or cone
   selection from topology and conservation alone; they do not refute every
   Euler environment with a boundary, background, constitutive scale, or other
   explicitly scale-breaking structure.

**Maximum verdict:** an exact sufficiency/selection boundary—topology may be
preserved while Euler scaling leaves the action unselected.

## Unit A check 3: compact-swirl pressure quadrupole

9. Start from the Euler pressure Poisson equation with the declared density and
   sign convention. Expand the Cartesian index contraction for the explicit
   compact swirl rather than importing a radial slogan.
10. Establish smoothness, compact support or decay, zero lower moments, and
   convergence of every source moment. Include all integration-by-parts and
   boundary terms at the axis, cutoff edge, and infinity.
11. Compute the first nonzero pressure multipole tensor, its trace-free
    quadrupole, angular dependence, sign, density normalization, and Green
    function coefficient. Check a concrete observation direction and symmetry
    transformation as exposing cases.
12. Differentiate the exterior pressure to the acceleration field. Verify that
    a nonzero pressure quadrupole gives `p=O(r^-3)` and acceleration
    `-grad(p)/rho=O(r^-4)`, with the advertised coefficient only where the
    exact tensor contraction is nonzero. Distinguish exact leading coefficient,
    big-O bound, and possible nodal directions.

**Maximum verdict:** the exact initial or stationary pressure quadrupole and
its convergent `r^-4` acceleration tail at the stated field scope, not a
Coulomb law or strict causal cone.

## Unit A check 4: classical hyperbolic extension and lattice thresholds

13. Write the extension's symplectic/KKS form, Hamiltonian, equations, and
    flow. Verify that the curvature/boost scale `K` and propagation parameter
    `c` are genuinely independent inputs and cannot be eliminated or equated by
    a hidden normalization.
14. Check whether the extension is derived from Euler, merely compatible with
    the classical observable algebra, or explicitly postulated. Preserve a
    useful classical construction without upgrading it to a physical cone.
15. Re-derive lattice absolute convergence in three dimensions. For a tail
    `|R|^-p`, verify the base sum threshold `p>3`; for its first spatial moment,
    verify `p>4`. Distinguish absolute, conditional, symmetric, and
    shape-dependent summation.

**Maximum verdict:** exact classical extension algebra and exact lattice
thresholds at their declared conditional scope.

## Unit B check 1: exact disjoint stationary composite

16. State the regularity, finite-energy class, supports, and ambient domain of
    each compact-velocity stationary Euler cell. Verify divergence-free Euler,
    pressure normalization, and smooth flatness at the velocity-support edge.
17. Prove that exterior pressure is constant on the relevant connected
    component and flat to the order required for gluing. Check whether compact
    velocity, compact vorticity, and compact pressure perturbation are all
    actually established rather than interchanged.
18. Audit the target's pressure-moment convention explicitly. If `p` denotes
    physical pressure, the isotropic moment identity must read
    `M_ik=-(delta_ik/rho_0) integral p`; the formula without `rho_0^-1` is
    correct only for kinematic pressure `p/rho_0`. Verify that exterior
    constancy and cancellation of all exterior pressure multipoles survive the
    chosen convention.
19. For pairwise disjoint translated/rotated supports, substitute the complete
    sum into Euler and show every cross-advection and pressure term vanishes
    pointwise, including collars between supports. Record whether the composite
    is exact stationary Euler or only initial data.

**Maximum verdict:** an exact smooth finite-energy stationary composite on the
declared domain, with its genuine exterior-pressure property.

## Unit B check 2: proper-octahedral rotation orbit

20. Enumerate the orientation-preserving octahedral group, its action on
    physical space and vector fields, the chosen carrier centers/orientations,
    and the stabilizer. Verify orbit size after quotienting duplicate
    configurations.
21. Distinguish the six `+/-` coordinate-axis sites from a free 24-site orbit.
    Each axis site has a `C4` stabilizer. For a generic tangent seed, the four
    stabilizer rotations are different tangents supported on the same physical
    component, and a nonzero seed can have zero projection onto the trivial
    representation. Disjointness of the six sites therefore does not prove
    that the full 24-element group average is nonzero.
22. Require one of three concrete constructions for a nonzero invariant
    tangent: first exhibit a nonzero `C4`-invariant projection and average over
    the six cosets; use 24 genuinely disjoint sites in a free proper-octahedral
    orbit; or prove directly that the overlapping stabilizer average is
    nonzero and lies in the tangent space of the chosen product leaf.
23. Check that every group image is the same Euler construction under an actual
    physical rotation and that parity/improper reflections are not silently
    included in the proper group.
24. Distinguish a finite physical rotation orbit from an internal state space,
    invariant linear mode, representation selected by dynamics, or exchange
    character.

**Maximum verdict:** the exact proper-octahedral physical orbit and stabilizer
at configuration-space scope.

## Unit B check 3: vorticity moments and corrected Hodge tail

25. Derive the compact-vorticity moment identities from the actual compact
    velocity and curl relation. Check which zeroth, first, second, and higher
    moments vanish automatically, which vanish by octahedral symmetry, and
    which depend on centering or orientation choices.
26. Include the translation of every component in the global moments. Verify
    `integral omega=0`, define the first and second vorticity moments with all
    spatial and vector indices visible, and derive their symmetry and trace
    relations before applying the group average.
27. Prove the proper-octahedral identities on those exact translated tensors.
    In particular, a rank-three averaged tensor vanishes only after using the
    symmetry in its first two spatial indices; isotropy prose alone does not
    rule out the proper-octahedral pseudoscalar tensor.
28. Expand the whole-space Biot--Savart/Hodge kernel with Cartesian indices and
    tensor symmetries. Include trace, antisymmetric, and divergence constraints
    before declaring the first surviving multipole.
29. Verify the repaired order: the claimed moment cancellations must make the
    exterior velocity `O(r^-5)` and its gradient `O(r^-6)`. Check whether this
    is an upper bound or a nonzero exact leading term, and test an allowed
    observation direction for accidental cancellation.
30. Audit the correction receipt, verifier, API, and repaired tests specifically
    for sensitivity to the old off-by-one multipole order. A literal expected
    exponent or copied tensor symmetry is regression evidence unless it
    evaluates the defining moments and kernel expansion.

**Maximum verdict:** the exact moment-cancellation theorem and corrected Hodge
tail order in the stated norm/asymptotic sense.

## Unit B check 4: DA product and lattice first moment

31. Construct dynamically accessible tangents from compact generators on each
    cell. Verify support, stabilizer, harmonic rows, pressure/Hodge tails, and
    whether cross KKS or Hessian terms vanish on the direct product.
32. For the same-field composite, specify the displacement on every component,
    preservation of its support/collar, and its circulation and centering rows.
    Only after these componentwise data are fixed does the assembled DA tangent
    define a tangent of the direct-product orbit.
33. Distinguish an exact product of stationary leaves from an invariant
    positive mode or a coupled Hamiltonian band. Exact disjointness may remove
    coupling rather than generate it.
34. Apply the corrected `r^-5` velocity or the actual claimed interaction tail
    to a three-dimensional lattice. Verify absolute convergence of the field
    and the `p>4` first-moment threshold, including uniformity over cell
    orientation and lattice summation shape.
35. Separate this pointwise lattice fact from a Bloch `C^1` operator theorem.
    The latter additionally requires one uniformly graph/domain-bounded
    coupling `J(R)` between normalized cell modes, with pressure/Leray
    derivatives and mode normalization obeying the same `r^-5` estimate.
    Pointwise decay of the velocity alone does not construct that operator.

**Maximum verdict:** a genuine DA/direct-product support statement and exact
lattice first-moment convergence boundary.

## Source, API, oracle, and correction protocol

After activation, source theorem hypotheses will be matched to the exact
domain, regularity, support, topology, symmetry, and parameter quantifiers used
by each target. Repository code is algebraic support only until its predicates
are bridged to those continuum claims.

For each verifier and test pair, the review will inventory assertion sites,
runtime result, status zero, wrong-sign/order/normalization sensitivity, and
the exact proposition each check can establish. The 0060 original and repaired
receipts will be distinguished; the repaired receipt is authoritative for the
corrected API boundary. Existing receipts will be reused while their frozen
inputs remain unchanged.

If one concrete defect appears, the review will state the strongest supported
result, exact unsupported extension, minimum repair, and evidence that restores
the stronger result. One bounded target correction package may receive one
changed-statement/dependency check. The correction protocol does not reopen
unaffected exact atoms.

## Verdict map and exclusions

The final record will give separate Unit A and Unit B verdicts plus a joined
compatibility verdict. Every route receives exactly one status: `established`,
`refuted` with its mechanism, or `blocked` with its next construction.
Conditional extension algebra and convergence implications remain conditional
and cannot inherit physical realization from adjacency.

The strongest supported output may include exact Euler similarity weights,
continuous-action/topology separation, a pressure quadrupole and acceleration
tail, exact classical extension algebra, exact convergence thresholds, an
exact stationary octahedral composite, or corrected Hodge multipoles. None of
those statements by itself supplies:

- an invariant positive Euler mode or finite-dimensional autonomous doublet;
- a Bloch band or propagating lattice dynamics;
- a strict finite-speed or Lorentz cone;
- numerical action selection, `hbar`, or a quantum representation;
- a particle, electron, neutrino, Born/reset, or exchange conclusion;
- P4 or parent-campaign completion; or
- evidence from active work outside `0055` and `0060`.

No nonlinear stability theorem, particle/quantum interpretation, registry
promotion, central edit, source edit, API/test edit, memory edit, or commit is
part of this review. The parent campaign remains active whatever the joined
verdict.
