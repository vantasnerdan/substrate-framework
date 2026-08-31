---
description: Independent constructive review of C-M5C-001 through C-M5C-004
author: codex-independent-p249-reviewer
created: '2026-08-31'
updated: '2026-08-31'
tags:
- substrate-framework
- claim-review
- P249
category: decisions
confidence: high
status: active
---

# Independent review of the P249 O6 claim transaction

This review retains the exact compact-clock construction and positive-action
algebra, but requests changes on all four submitted claims because the full
six-component tensor contains an omitted charged shear doublet.  That omission
changes both the Noether inertia and the lightest charged continuum edge, so it
also removes the strict hylomorphy premise used for existence and stability.

## Claim and Positive Role

The frozen transaction proposes four additive claims for the canonical
auxiliary-frame completion in `claim-draft.yaml`:

- `C-M5C-001` asserts a faithful diagonal compact action on an arbitrary
  symmetric spatial tensor `S` and a complex scalar, exterior fixation,
  Noether charge `Q_clock=omega I`, and fixed-charge Legendre reduction with
  `I=int(|psi|^2+4||B||^2)`.
- `C-M5C-002` asserts global positivity, the complete six-channel exterior
  pencil, and the charged generator-frequency edge `19/4`.
- `C-M5C-003` uses a `39/16` split plateau and Benci--Fortunato Theorem 18 to
  assert a translation-compact full-field fixed-charge minimizer with an
  intrinsic open-space width and positive-mass tails.
- `C-M5C-004` promotes that minimizing set to full fixed-charge Hessian
  nonnegativity, Lyapunov orbital stability, and the stated radiation scope.

These are meaningful positive roles and are reviewed at their submitted
full-field scope, not reduced to the earlier two-amplitude ansatz.

## Frozen Transaction

The review boundary is repository `HEAD`
`386763ef2b18434ad29892e1ac85b499e6501d46` plus the following uncommitted
content hashes:

- claim draft: `c7bcdcc849b37d496627f4df32efde9df86d31bf5132c0c50ab7a3d48142e1c9`
- attempt-0006 manifest: `6585fbbb6a71eb2ca840fe7940bc51c880b4162af01d9cb403b3ede110de74f1`
- attempt-0006 derivation: `3b6d6cab386956c47d41e354658632508dda8366eb752f193198a24497950429`
- attempt-0006 verifier: `46acb1cd7576b3ddb8821ea06614579648219233f71c1ed7493c06eb46d106f8`
- theorem bridge: `54f90297cee5d2c4dc81013eca0f3509cdb0b5aef52f4fb185dda339d19b84b6`
- dependency audit: `68c630ae469438a787e959c58ae91c9091c2fbd7f512058b2ba04196475ab4ce`
- canonical module: `e9c45e8ee91d8887c81c9cfd9cd6f230209ad4f89f0df42a30348bfb304b2a5c`
- canonical tests: `e5ae7a8660145af003cbe66753066693f27a86b3284e12bc59473b69eb86b55a`
- companion checks: `7f8bd8673417f80c57d6133ba65506db261159719b72af6b1fda98a17f9e22eb`

The changed implementation is the additive module
`src/substrate_framework/m5_exterior_clock.py`; no accepted claim or existing
canonical symbol is changed yet.  The declared external dependency is
Benci--Fortunato, arXiv:1212.3236, Theorem 18.  Direct affected consumers are
the new canonical test, companion exact checks, four claim records, and their
eventual registry/release/docs/memory records.  Existing adjacent P239/P240
claims and modules are outside the scientific review boundary.

## Strongest Supported Positive Statement

The submitted artifacts exactly support a useful full-field core: the
potential is nonnegative with the unique aligned exterior; the diagonal
`S1` action is smooth, faithful through `psi`, and fixes that exterior; the
six-channel Hessian and kinetic metric printed in the draft are correct; the
split plateau has the stated local potential and inertia values when the
axis--tangent shear is zero; and the complete fixed-charge differentiation
formula is correct for the *correct* inertia functional.

For arbitrary symmetric `S`, however, the generator also rotates
`(u,v)=(S_12,S_13)` with weight one.  The exact tensor contribution to the
inertia is

`Tr((A S-S A)^2)/2 = u^2+v^2+4p^2+4q^2`,

not only `4||B||^2=4p^2+4q^2`.  With the submitted axis-lock coefficient, the
same weight-one shear doublet has generalized mass square `1`.  The true
charged generator-frequency edge of the submitted action is therefore `1`,
not `19/4`, and the `39/16` plateau lies `23/16` *above* that edge.  The
strongest supported result consequently stops before theorem-backed existence,
localization, or stability.

## Evidence Map

The evidence groups contribute at the following exact scopes.

| Evidence | Proposition established | Role | Bridge to claim | Limit |
| --- | --- | --- | --- | --- |
| Attempt 0006 positivity and exterior-Hessian algebra | Global M5 lower bound, unique aligned exterior after positive locks, Hessian `diag(5,7,19,1,1,19)`, kinetic metric `diag(1/2,1,1,1,1,1)` | exact proof | Supports the positivity and pencil core of C-M5C-002 | Does not classify the six channels under the compact generator |
| `m5_exterior_clock.py` plus canonical tests | Exact basis, equivariance of `Q(psi)`, potentials, plateau values, Legendre derivatives | regression/exact subclaims | Supports implementation of parts of C-M5C-001 and C-M5C-002 | `clock_inertia_density` omits the weight-one shear pair and the tests do not vary it |
| Attempt 0006 and companion exact checks | Submitted plateau ratio `39/16`, B-unsplit lower bound `5`, source `-12 Q(psi)`, Derrick formula | exact proof of those algebraic predicates | Supplies candidate binding data for C-M5C-003 | Compares against the wrong charged edge, so it does not establish strict hylomorphy |
| arXiv:1212.3236, Theorem 18 | Conditional existence and stability of hylomorphic solitons when `E,C` on an FT dynamical phase space satisfy EC-0 through EC-3 and strict hylomorphy | external exact theorem | Could bridge C-M5C-003 and C-M5C-004 after all hypotheses are discharged | The submitted hylomorphy inequality is false for this action, and the bridge records only an `H1` configuration space rather than the full canonical phase space |
| Adjacent canonical tests | Existing P239 action and kinetic-axis consumers still pass | regression | Confirms the additive module did not break those consumers | Does not validate the new claim semantics |

## Oracle Audit

The submitted exact verifier reports
`ALL 40 CHECKS PASS [P249-O2-O5-0006]`; the companion reports
`ALL 22 CHECKS PASS [P249-O6-CANONICAL]`; and the targeted test command reports
`29 passed`.  Repository validation also passes with 248 accepted claims.
These are valid execution receipts but false green for the headline
transaction because none audits the generator weight of `(u,v)`.

An independent exact generator calculation gives

```text
delta S = A S-S A
        = [[0,-v,u],[-v,-2q,2p],[u,2p,2q]]
Tr(delta S^2)/2 = u^2+v^2+4p^2+4q^2.
```

Combining this with the verifier's own generalized masses
`(10,7,19,1,1,19)` in `(a,t,p,u,v,q)` order proves that `(u,v)` is a
charge-one mass-square-`1` physical doublet.  Hence
`min(6,1,19/4)=1`, while `39/16-1=23/16>0`.  This is a direct exact
counterexample to the submitted inertia and edge predicates, not a numerical
resolution issue.

The existing mutations delete the M5 quartic, selected kinetic entries, or the
phase lock.  They do not mutate the group representation or test every charged
channel, which is the concrete reason the tally remains green.

The external source was checked at its primary statement.  Theorem 18 assumes
`C1` energy and charge functionals on an FT dynamical system and obtains stable
fixed-charge minimizing sets after EC-0 through EC-3 and strict hylomorphy.
Its concrete field-theory setup includes canonical momentum/velocity variables
in the phase space.  P249 instead names only
`H1(R^3;Sym(3)+C)` and replaces the energy/charge ratio by a reduced profile
`J/K` comparison without explicitly defining the full phase-space `E,C`, the
canonical flow, or their EC properties.  A correct Legendre reduction can
potentially supply that bridge, but the current artifact does not.

## Findings

Each finding below is inside the frozen transaction and directly affects a
submitted statement or its dependency closure.

| Finding | Direct evidence | Blocking in boundary / minimum correction / follow-up | What would resolve or overturn it |
| --- | --- | --- | --- |
| F1: C-M5C-001 omits a physical contribution to its Noether inertia | Exact generator calculation gives `I_S=u^2+v^2+4p^2+4q^2` for arbitrary symmetric `S` | **Blocking.** Preserve the compact action, exterior fixation, and Noether/Legendre construction, but correct `I`, the canonical API, tests, and all dependent formulas to include the axis--tangent shear doublet | An exact derivation showing that `(u,v)` is constrained out of the declared physical field space, or the corrected full inertia with equivariance and Legendre checks |
| F2: C-M5C-002's charged edge `19/4` is false for the submitted action | `(u,v)` has weight one and generalized mass square `1`; the draft's own `K,G` matrices establish the mass | **Blocking.** For the current action the edge must be `1`.  If the action is repaired to lift this channel, recompute the complete Hessian, weights, inertia, edge, positivity, and witness in one corrected transaction | A corrected action/constraint whose complete generator-resolved pencil has no lighter charged channel, verified by a channel-by-channel weight test |
| F3: C-M5C-003 does not satisfy the cited theorem's strict hylomorphy premise | `39/16 > 1` by `23/16`; the theorem bridge compares only to the omitted-channel edge `19/4` | **Blocking.** Retain the exact plateau, unsplit lower bound, B source, fixed-charge calculus, and Derrick identity, but do not infer a minimizer.  Produce a strict witness below the *complete* charged edge or a different exact existence proof | A corrected full action and witness with `J/K` strictly below every charged vanishing channel, plus a complete concentration/compactness or theorem bridge |
| F4: the external theorem is not yet typed on the declared full canonical phase space | The primary theorem assumes dynamical-phase-space `E,C`; `theorem-bridge.md` supplies only the configuration `H1` space and assertions of splitting/coercivity | **Blocking for C-M5C-003 and C-M5C-004.** Define the canonical phase space including momenta, conserved energy and charge, well-posed flow, and show that reduction to `E_Q` preserves the theorem's minimizer and stability conclusions | A sourced or self-contained phase-space construction discharging EC-0 through EC-3, strict hylomorphy, flow invariance/well-posedness, and the relative-equilibrium reconstruction |
| F5: C-M5C-004 misclassifies the mass-square-`1` shear doublet as neutral and inherits an unproved minimizing set | The compact action rotates `(u,v)` with weight one; C-M5C-004 calls mass squares `10,7,1,1` neutral and uses C-M5C-003 | **Blocking.** Recompute active harmonics and radiation inequalities from the corrected representation, then apply the stability theorem only after C-M5C-003 is repaired | A complete corrected charged/neutral decomposition with all charged harmonics below their thresholds and a valid full-phase-space minimizing-set stability bridge |

No prior independent review has narrowed these claims, so there is no recursive
scope-shrink history to escalate.  The findings arise from new contrary exact
evidence and name a finite repair path.

## Compatibility and Consumers

The compact action, positive kinetic convention, mostly-plus Lorentz typing,
and aligned exterior are compatible with the stated conditional auxiliary-frame
extension.  No M7 threshold or empirical comparator is imported.  The module is
additive and adjacent existing tests pass.

The defect is internal to the new transaction: the implementation and both new
verifier surfaces encode an incomplete generator representation.  Therefore
the new canonical test and companion exact checker are affected consumers that
currently pass the wrong proposition.  Registry, release, generated docs, and
accepted memory must not consume C-M5C-001 through C-M5C-004 until correction.

## Four-Axis Decision

The decisions are claim-by-claim and do not reject the candidate concept.

### C-M5C-001 — REQUEST_CHANGES (high confidence)

- Verification: `symbolic_verified` for compact action/exterior fixation;
  contradicted as submitted for the full inertia formula
- Review: `audited`, not accepted
- Compatibility: `compatible_extension` in its retained conditional scope
- Epistemic: `proposed`
- Relationship: additive conditional claim
- Strongest accepted or proposed statement: the full field has a faithful
  exterior-trivial diagonal `S1` and a Noether/Legendre map with inertia that
  must also contain the weight-one shear norm

### C-M5C-002 — REQUEST_CHANGES (high confidence)

- Verification: `symbolic_verified` for positivity and the displayed `K,G`;
  contradicted as submitted for charged-channel classification and edge
- Review: `audited`, not accepted
- Compatibility: `compatible_extension`
- Epistemic: `proposed`
- Relationship: depends on corrected C-M5C-001
- Strongest accepted or proposed statement: the submitted action is positive
  with generalized masses `(10,7,19,1,1,19)`, and its complete charged edge is
  `1` unless the action or physical constraints are changed

### C-M5C-003 — REQUEST_CHANGES (high confidence)

- Verification: `unverified` for existence/localization; exact subclaims only
- Review: `audited`, not accepted
- Compatibility: `unassessed` pending a valid existence bridge
- Epistemic: `proposed`
- Relationship: depends on corrected C-M5C-001 and C-M5C-002
- Strongest accepted or proposed statement: exact positive functional,
  plateau data, unsplit bound, nonzero B source, Legendre calculus, and Derrick
  identity; no full-field minimizer follows from the current evidence

### C-M5C-004 — REQUEST_CHANGES (high confidence)

- Verification: `unverified` for orbital stability/radiation; conditional
  variational principle only
- Review: `audited`, not accepted
- Compatibility: `unassessed`
- Epistemic: `proposed`
- Relationship: depends on a repaired C-M5C-003 minimizing set
- Strongest accepted or proposed statement: if a full canonical
  fixed-charge global minimizing set is established with well-posed conserved
  dynamics, its constrained Hessian is nonnegative and the standard
  energy-charge argument is available; those premises are not yet established
  here

## Promotion Transaction

No claim in this frozen transaction is promotion-eligible.  A bounded repair
must update the action/inertia representation, generator-resolved pencil,
binding witness, theorem bridge, canonical module/tests, companion checker, and
the four draft statements together.  Registry, campaign adjudication, release,
generated docs, and accepted memory remain unchanged until the correction
check passes.

## Correction Check

The correction check is pending and must be limited to F1--F5 and their direct
dependency/consumer edges.  It should verify: the exact full generator norm;
all channel weights and generalized thresholds; a strict binding inequality
against the complete edge; the full canonical phase-space theorem hypotheses;
and the resulting charged/neutral stability and radiation statements.  It is
not a second substantive review of unrelated M5 claims.

## Result and Frontier

The review preserves a strong algebraic candidate: a positive canonical
full-spatial M5 completion with a genuine compact exterior-trivial action and
an exact split-core coupling.  The decisive missing construction is a
generator-complete binding and phase-space existence bridge.  Correct the
weight-one shear sector before any promotion or numerical stability work; the
parent positive objective remains open on that finite analytic repair.

## Cross-References

The reviewed boundary is issue #189, proposal P249, attempts 0001--0006
(especially 0006), `claim-draft.yaml`, `evidence/theorem-bridge.md`,
`evidence/dependency-audit.yaml`, `src/substrate_framework/m5_exterior_clock.py`,
`tests/test_m5_exterior_clock.py`, and `companion/exact_checks.py`.  The external
theorem source is <https://arxiv.org/abs/1212.3236>.

## Correction Check — completed

This is the single bounded correction check authorized for F1--F5.  It does
not reopen the original substantive review or inspect unrelated claims.  The
repaired boundary is pinned by these hashes:

- claim draft: `af10c01d55ef6c7846f8db7d14a97e555c44a77394437eabc28794319eef73c8`
- canonical module: `c484acb1ada8e746663123f148ee56f172c32d0e6f4b55fa02becd8eaf518777`
- canonical tests: `4db9e44d658dad02928b19b72c197cba33310d106bbd57e1c75181a8905aba6b`
- companion checks: `ce095386917d3cc31ce5c3d41c1399513c60fa32dbcfd66cff35677771af3e0f`
- theorem bridge: `e7325af2cd8fc0e7d5ffd917fec7bbcbdf27bb2b7e7e68065391961b45273354`
- attempt-0008 manifest: `6c7e20c8ea2d2223c552364b113b7dab7cfab5e28e1a35346ec2252b267d0f10`
- attempt-0008 derivation: `fbef866510166ba06d014f38cb17bc1ebd961d23b9a1f729ffcd1e71daf5a5a7`

The direct correction results are:

| Original finding | Correction evidence | Result |
| --- | --- | --- |
| F1, omitted shear inertia | `clock_inertia_density` now computes `Tr([A,S]^T[A,S])/2`; independently and in tests it gives `u^2+v^2+4p^2+4q^2`.  The canonical phase-space APIs give `C=omega I` and `E_vel=omega^2 I/2` exactly. | Resolved |
| F2, incomplete charged edge | The unit axis lock gives generalized masses `(10,10,22,4,4,22)`.  The complete weights are scalar one, shear one, and clock tensor two, hence `min(6,4,22/4)=4`. | Resolved |
| F3, failed hylomorphy | The repaired plateau has `V=45/64`, `I=1/2`, and squared optimized ratio `45/16`, strictly below `4` by `19/16`.  The full `B=0` ratio is at least `4`, so a minimizing state below the edge necessarily has a nonzero clock split. | Resolved |
| F4, untyped theorem phase space | The bridge now defines `X=H1 x L2`, the complete positive Hamiltonian, bilinear Noether charge, canonical velocities, EC-0 through EC-3, the translation-vanishing edge, and reconstruction of the relative equilibrium by equivariance and uniqueness of the energy-subcritical canonical flow. | Resolved |
| F5, wrong radiation representation | The repaired draft classifies the shear doublet as weight one, leaves only the two mass-square-10 channels neutral, and uses the corrected inequalities `6-omega^2>2`, `4-omega^2>0`, and `22-4omega^2>6`. | Resolved |

For the strict frequency step in F5, the primary theorem's constructed state is
a minimizer of `J_delta=Lambda+delta E`.  At positive charge its first
variation gives

`E'=[Lambda/(1+delta C)] C'`.

The canonical phase-space reduction identifies this Lagrange multiplier with
`omega`.  Since the repaired hylomorphy threshold is `Lambda_0=2`, this gives
`0<omega<2`, and hence all three strict tail inequalities in the corrected
claim.  This closes the objective bridge rather than inferring the frequency
from the plateau value alone.

The corrected oracle receipts are:

- `ALL 28 CHECKS PASS [P249-O6-CANONICAL]`;
- `31 passed` for `test_m5_exterior_clock.py`,
  `test_m5_covariant_action.py`, and `test_m5_kinetic_axis.py`;
- repository validation: `WORKFLOW VALID: 248 claims, 248 accepted, 12
  proposals, 4 skills; MIGRATION QUEUE: 218 units, 0 pending, 0 partial`.

The correction-check decisions supersede the earlier request-changes verdicts
at exactly this repaired boundary:

### C-M5C-001 — ACCEPT (high confidence)

The exact compact action, exterior fixation, full Noether inertia, regular
Legendre map, and `dE/dQ=omega` are symbolically verified on the declared
six-component tensor plus scalar field.  Final axes are
`symbolic_verified / accepted / compatible_extension / active`.

### C-M5C-002 — ACCEPT (high confidence)

Global positivity, the unique aligned exterior, complete positive kinetic
pencil, generator weights, generalized masses, and charged edge `4` are
symbolically verified without an imported threshold.  Final axes are
`symbolic_verified / accepted / compatible_extension / active`.

### C-M5C-003 — ACCEPT (high confidence)

The corrected strict hylomorphy witness and canonical phase-space theorem
bridge establish a nonempty translation-compact family of full-field
fixed-charge relative equilibria with necessary nonzero M5 core split, finite
energy/inertia, intrinsic open-space localization, positive-mass tails, and
the stated Derrick identity.  Final axes are
`symbolic_verified / accepted / compatible_extension / active`.

### C-M5C-004 — ACCEPT (high confidence)

The full minimizing set has nonnegative constrained second variation and the
energy-charge Lyapunov argument applies on the declared well-posed canonical
flow.  The repaired generator-resolved channel decomposition supports the
stated subthreshold radiation scope.  Final axes are
`symbolic_verified / accepted / compatible_extension / active`.

All F1--F5 correction debts are closed.  These decisions authorize the
claim-promotion transaction; registry insertion, release pinning, generated
docs, memory synchronization, and final validation remain operational O6 work
rather than further scientific review.

The companion hash above was reconciled operationally after two transient
extra assertions were removed to retain the reviewed 28-check surface.  The
current checker still exercises the same F1--F5 propositions and reports
`ALL 28 CHECKS PASS [P249-O6-CANONICAL]`; this hash correction is record
reconciliation, not a new substantive review pass.
