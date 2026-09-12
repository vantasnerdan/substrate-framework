# P253/0101: independent final review of P253/0100

## Frozen review object

This README preregisters one content-blind, fixed-boundary, independent
non-author/non-implementer review of final root-owned `0100`. The review has
nine separate one-verdict units:

- **A0 — unforced Ertel material invariant/current:** derive the unforced
  material invariant, local continuity current, true-scalar parity, and exact
  smoothness, divergence and boundary/decay domain;
- **A1 — forced Ertel flux:** derive the forced source/flux law with its exact
  sign, curl/gradient contractions and boundary domain;
- **B0 — closed-line real constant lock:** adjudicate the kinematic
  refutation of a nonzero real constant coefficient on a closed line;
- **B1 — exact forced lock condition:** independently derive and domain-check
  `chi D_t lambda=(curl f) dot grad chi`, keeping it distinct from the
  unforced condition `D_t lambda=0`;
- **B2 — force-selected lock:** determine whether any force actually selects
  a nontrivial admissible `lambda`, preserving that construction as open when
  only the compatibility equation is supplied;
- **C0 — Cao regular density/convective current:** check the claimed zero
  regular density and convective part on the Cao carrier;
- **C1 — charged Cao magnetization:** audit the magnetization superpotential
  as a possibly nonzero identically conserved flux, distinct from C0;
- **D1 — punctured-domain azimuthal phase:** adjudicate only the equilibrium
  integral and axisymmetric no-swirl scope, with
  `D_t theta=u_theta/r` as the dynamical boundary; and
- **D2 — moving defect continuation:** keep a moving advected defect and
  nonaxisymmetric physical construction separate unless the target supplies
  its domain, dynamics and same-carrier map.

Each unit receives exactly one route-scoped verdict: `established as stated`,
`refuted` with its mechanism, or `blocked` with its missing construction. No
unit inherits a verdict from a sibling, and this finite unit list decides
neither P6 nor the parent campaign.

P253/0095 is excluded in full. Its body, module, tests, verifier, receipts,
result, source audit and conclusions will not be opened or used. Every parent,
P2/P4/P5/P6, particle, Born/reset, exchange/statistics, electron and neutrino
conclusion is outside this review.

## Independence and activation boundary

`particle-balance-review` authored or implemented none of `0100`. It owns
only append-only review artifacts under `attempts/0101` and will not edit
`0100`, source/API/test files, central/generated records, memory, or Git
history.

Before this freeze no `0100` README body, derivation, source audit, result,
validation, verifier, manifest, receipt, command, output, API, or test body
was opened. Inspection was limited to central registry metadata, filenames
needed to identify the target inventory, and content-blind SHA-256 hashes.
The excluded attempt was not opened. An earlier unactivated `0101` draft used
the pre-forced-lock inventory; it was never reported for activation and is
superseded by this README.

The target may be opened only after the coordinator:

1. centrally registers `0101` with this exact reviewer and scope;
2. records the repository-schema command, stdout, stderr, and exit under
   `attempts/0101`;
3. makes `attempts/0101/activation-schema.exit` contain exactly `0`; and
4. acknowledges activation against this exact README hash.

After activation the reviewer will first recheck the frozen hashes, then open
the final manifest, completion and continuation receipts before the claim
bodies. One substantive pass may request at most one evidence-driven bounded
correction package, followed by one correction-only check. Final writes are
limited to `0101/review.md` and `0101/verdicts.yaml`, followed by a report only
to `w4:p2`.

The authority base at freeze is accepted release `v0.183.0`. The frozen
SHA-256 values of `governance/releases/current.yaml`,
`governance/claims.yaml`, and the P253 proposal registry are
`7b889a3c0186c0c22e698fd4ccc69e19066b8262018696cfc41e573ec2f0ffd8`,
`78dbe694359739ee15f7db2b2cb1ccde7f48b5af34f6d9cd451282732dd08983`,
and
`4c98a7a4c3dd296f7a5f8588fc9f4255a2dfedae7abed1d317bdb4563a6ddea4`.
The observed repository head is
`3e93e41a236b42b4f10db9911cbbd79677590026`. Unrelated shared-worktree
changes are preserved and excluded.

This is an exact Euler/material-invariant, forced-current, topology,
axisymmetric carrier, superpotential, punctured-domain and evidence-role
audit. No production numerical design or run is frozen. Exact scripts may
expose vector identities, signs and integral cancellations but cannot by
tally prove global defect dynamics, force selection, nonaxisymmetric
continuation, a physical weak current, or a particle mechanism.

## Frozen final P253/0100 inventory

The final content-blind selectors are:

| Artifact | SHA-256 |
|---|---|
| `0100/README.md` | `fa0a44de2aff51b1cccab740fe35dc2688d4fb3673c2552f5c411e158066962b` |
| `0100/derivation.md` | `3998a2f55bc3f819170b0ba0a51baf5728594e64c2fe13d2ead1bbc1fb38369b` |
| `0100/result.yaml` | `f1edc4ca9148cc930a44e312ade60ed03ea380a09e786d93ef394546466e7697` |
| `0100/source-audit.md` | `01c32ef1eed59c9fcbe7907caa5d1ae122790e92dcc6ae6012de5059c80d4ccf` |
| `0100/validation.md` | `2fa59c5788a07c1902770d804208d520d77883a5f12c734d1a8510510370e4d0` |
| `0100/verify_ertel_current.py` | `acb960106b55886e52ea5b6314af4edbe2616cc8b8e66390ed42e52acdded4f9` |
| `0100/author-completion-receipt.md` | `53fde48bb1495daee2a4d055f504cd501a1d7d742a2b80a0ccb1f5e81a825397` |
| `0100/failure-derived-continuation-receipt.md` | `c32d1298554a36074fb38fdde67269d67dffb0fd22c8906784bd2aef9029a3ed` |
| `0100/artifact-hashes.sha256` | `dff6ddccd1166019720f6a545344d246195075f9cf462833ea68886c6af8d54f` |
| `src/substrate_framework/euler_ertel_current.py` | `7f1f8233416d2bb817818fdabc72a27480864e86740269cd133483f7497c2ec3` |
| `tests/test_euler_ertel_current.py` | `f2d200423d5f36cdd4cb83299ca712f90b027f598a8a85d6e9a2a65f89ed1bb4` |

The final receipt set is pinned through the manifest:

| Final receipt | SHA-256 |
|---|---|
| `0100/exact-v5.command.txt` | `aadc3d807c45e48274609967ce7c70814b16cff773e7152c2e88be70da09bf93` |
| `0100/exact-v5.elapsed.txt` | `269aa5bbe6bb07f6a2e5dd0f959ba2eace7a546423324991c57df1fb1a34bfa2` |
| `0100/exact-v5.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0100/exact-v5.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0100/exact-v5.stdout.txt` | `f774ef4f01b073bb4bbef97001461523fc72408f5ca345997a750e3f43f0303a` |
| `0100/focused-v5.command.txt` | `11f70c0cc6052cb9f1f1a648874d1e3be865409bc55834f30a1f3b8710d3bda0` |
| `0100/focused-v5.elapsed.txt` | `83dcdc58f86b7badd27b427c487c06d27dacc93c88e9164d4d0b7d2d39c03c9e` |
| `0100/focused-v5.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0100/focused-v5.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0100/focused-v5.stdout.txt` | `1ce03978dbcbfe5f9195563e532d113b23517586b892d27d23ba9a460466a9ab` |
| `0100/repository-v5.command.txt` | `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f` |
| `0100/repository-v5.elapsed.txt` | `274449b2fe7ff8b940768d71728840444a4e0d3a1639c7fc01bee928744e4123` |
| `0100/repository-v5.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0100/repository-v5.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `0100/repository-v5.stdout.txt` | `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f` |

The SHA-256 of the sorted top-level `0100` filename list is
`13af0a6fd5f9017fec64a3c21b992c7374599e39fcc92df0f6deac936eb19efb`.
The pre-forced-lock manifest `07074a74...` and completion receipt
`68869ab7...`, together with earlier v1--v4 receipts, are historical
chronology rather than final evidence. No differently named target module or
test is presumed.

### Mixed-base manifest convention

The final manifest is treated as a mixed-base selector: bare attempt artifact
names resolve under `attempts/0100`, while `src/` and `tests/` entries resolve
from the repository root. After activation each entry must match at its own
declared base. Historical receipts establish chronology, not scientific
authority; the final v5 receipts govern the completed package.

## Frozen review method

### Unit A0 — unforced Ertel material invariant and current

1. Starting from smooth incompressible Euler and a true transported scalar,
   derive the vorticity and scalar-gradient material equations and recompute
   `D_t(omega dot grad chi)` without assuming the cancellation.
2. Convert material invariance to the exact local continuity current using
   incompressibility. State the smoothness, boundary-flux and decay hypotheses
   for local and integrated conservation separately.
3. Check parity under `O(3)`: `chi` is a true scalar, `grad chi` is polar,
   `omega` is axial, and the contraction must use the resulting exact
   transformation law consistently.

### Unit A1 — forced Ertel source and flux

4. Derive the forced vorticity equation from the actual force convention and
   contract with `grad chi`. Verify every curl, divergence and
   integration-by-parts sign.
5. Distinguish a local source from a flux divergence and state the boundary
   conditions under which the integrated quantity changes by the claimed
   term.

### Units B0, B1 and B2 — lock compatibility and selection

6. Integrate the proposed real constant lock around a closed
   material/vorticity line with orientation, period and regularity explicit.
   Verify the kinematic mechanism refuting a nonzero constant coefficient.
7. In the unforced system independently derive `D_t lambda=0`. Do not extend
   the constant-lock refutation to an admissible variable advected
   coefficient; audit its loop monodromy and zero-mean constraint.
8. In the forced system independently derive
   `chi D_t lambda=(curl f) dot grad chi`, including the zero set of `chi`,
   regularity, closed-line solvability and boundary domain. Do not conflate
   this condition with unforced material advection.
9. Determine separately whether the target constructs a force and coefficient
   satisfying the forced condition on the claimed carrier. If not, preserve
   force-selected lock as the named next construction. Neither compatibility
   law supplies compact character, quantization, chirality, spin or flavor.

### Unit C0 — Cao regular density and convective current

10. Evaluate the Ertel density and convective current on the actual regular
    axisymmetric no-swirl Cao field and tag class. Check cylindrical
    components, axis behavior, parity, support, and whether each vanishing is
    pointwise or only integrated.
11. Keep this carrier-specific vanishing distinct from a universal
    Euler-current no-go and from the magnetization term in C1.

### Unit C1 — charged Cao magnetization superpotential

12. Reconstruct the antisymmetric magnetization tensor or curl current and
    verify its divergence identity, gauge/domain regularity, compact support
    or decay, and physical factors.
13. Determine whether the flux is identically conserved, dynamically
    transported, or sourced. Its possible nonzero value must not be erased by
    C0's zero density/convective part; conversely a superpotential does not by
    itself yield nonzero net charge, weak coupling, radiation, or a particle
    observable.

### Unit D1 — punctured equilibrium phase

14. Specify the punctured domain, branch of azimuthal angle, axis exclusion,
    loop orientation and admissible tag/field regularity. Check the equilibrium
    integral independently.
15. Enforce `D_t theta=u_theta/r`. The equilibrium conservation is claimed
    only for axisymmetric no-swirl flow, where this term vanishes; it is not an
    advected-phase theorem for general Euler evolution.

### Unit D2 — moving advected defect

16. Audit whether the target constructs a moving puncture/defect, material
    transport of its branch data, finite-energy core, nonaxisymmetric domain,
    and same-carrier physical current. If not, name that object as the next
    construction without weakening D1.

### Evidence and agreement

17. Inspect the final v5 verifier and focused module/test receipts separately.
    Record which signs and identities are independently derived, which are API
    regressions, and the maximum verdict each can support.
18. Reconcile the failure-derived continuation, historical v1--v4 receipts,
    final completion receipt, mixed-base manifest, source audit, validation,
    result and derivation for stale or crossed claims.
19. Do not rerun unchanged evidence merely for counts. Only a predicate defect
    may trigger the single bounded correction.

## Final decision contract

The final review will publish separate verdicts for A0, A1, B0, B1, B2, C0,
C1, D1 and D2, followed by the strongest combined theorem and its next
construction. Every route receives exactly one verdict. Explicit hypotheses,
carrier-specific zeros, a potentially nonzero identically conserved
magnetization flux and honest parent exclusions are not defects. A bounded
correction will be requested only from concrete contrary evidence and will
preserve every stronger true statement that survives.
