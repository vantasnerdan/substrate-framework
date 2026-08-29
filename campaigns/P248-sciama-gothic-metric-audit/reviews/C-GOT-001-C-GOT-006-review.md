# Independent review of C-GOT-001 through C-GOT-006

## Claim and Positive Role

The frozen transaction proposes six exact claims. C-GOT-001 preserves the
paper's valid determinant-slaved metric as a constrained theorem. C-GOT-002
adds the missing lapse and proves a complete ADM-domain bijection. C-GOT-003
states the exact material/harmonic compatibility law. C-GOT-004 fixes the weak
Newtonian energy-sign ledger. C-GOT-005 supplies a conformally separated shear
variable. C-GOT-006 composes the complete map with accepted field and worldline
claims without promoting a microscopic ether interpretation.

## Frozen Transaction

The review boundary is release `v0.169.0` at
`443f1bd6edf0c1997f42b05c3b6719142e7a11c1`, `claim-draft.yaml`, the additive
`optical_gothic.py` API and tests, P248 companion verifiers and captured
outputs, `OpticalGothic.lean` and its axiom audit, source adjudication records,
and the impact replay. The only accepted propositions composed are
`C-STG-001`, `C-WLN-001`, and `C-WLN-002`; their acceptance is not reopened.

## Strongest Supported Positive Statement

The strongest supported result is all six claims after the bounded C-GOT-006
correction. In chart coordinates `q=(N,v,gamma)` the nonzero exact Jacobian
transfers metric Euler equations to the optical variables. `C-STG-001` supplies
the Einstein-scalar equations, `C-WLN-001` supplies the massive square-root and
Hamiltonian sector, and `C-WLN-002` supplies the massless null and affine
geodesic sector. No microscopic or empirical interpretation is imported.

## Evidence Map

The evidence roles remain deliberately separate.

| Evidence | Proposition established | Role | Limit |
| --- | --- | --- | --- |
| `exact_metric_checks.py` and `optical_gothic.py` | block map, determinant, inverse, Jacobian, round trip, variational chain | exact proof | no material realization |
| `exact_compatibility_checks.py` | compatibility identity, sign conflict, shear split | exact proof | time gauge only; no dynamics |
| `OpticalGothic.lean` | finite compatibility, sign, shear, nonzero-factor, and composition logic | corroborating subclaim | does not encode the full block matrix calculus |
| `numerical_checks.py` | SPD regressions, printed SI intervals, forced-wave convergence | regression/applicability | binary64 and declared meshes/tolerances |
| literature and dependency audits | exact scope of historical and accepted imports | provenance | not a proof of the new algebra |

## Oracle Audit

The reviewer reused the unchanged receipts: 31 symbolic runtime passes, 18
numerical passes, eight SPD round trips below `2e-13`, a second-order forced-wave
ladder with a rejected sign mutation, 62 affected tests, and a successful Lean
build. The original and corrected `commonMetricComposition` theorem both had an
empty axiom footprint. Check-call, mutation, runtime, and assertion inventories
are kept distinct in `evidence/verifier-audit.yaml`.

## Findings

The substantive review made individual decisions and found two blockers only
in the original C-GOT-006 wording.

| Claim | Initial decision | Confidence | Finding |
| --- | --- | ---: | --- |
| C-GOT-001 | ACCEPT | 0.99 | Exact constrained map, image constraint, and excluded general-metric witness close. |
| C-GOT-002 | ACCEPT | 0.99 | Exact complete inverse and nonzero Jacobian close; numerics are regression only. |
| C-GOT-003 | ACCEPT | 0.99 | Exact identity and counterexample close, with Lean corroboration. |
| C-GOT-004 | ACCEPT | 0.98 | The weak energy sign contradiction and separate stress coefficient close. |
| C-GOT-005 | ACCEPT | 0.99 | Additive counterexample and log-unimodular replacement close. |
| C-GOT-006 | REQUEST_CHANGES | 0.99 | The statement named `V` while quoting the `v=V/c` Jacobian, and it credited C-WLN-001 with a timelike-geodesic result absent from that accepted statement. |

The minimum correction changed the chart to `q=(N,v,gamma)`, described
`v=V/c`, replaced the unsupported massive timelike-geodesic phrase by the
accepted massive square-root/Hamiltonian sector, and removed the assumed
massive-geodesic bridge from the Lean composition. An exact V-coordinate
Jacobian with its `c^-3` factor and an accepted massive coordinate
Euler-Lagrange theorem would restore the two stronger formulations if desired.

## Compatibility and Consumers

All claims are compatible extensions. They preserve mostly-plus signature,
`x0=c*t`, exact inputs, SPD spatial slices, positive lapse, declared flow
orientation, and the accepted dependency scopes. No existing canonical symbol
was changed. The additive module, formal umbrella, accepted-dependency tests,
linearized-Einstein tests, and Gordon-scalar consumers form the replay boundary;
the targeted replay passed 62 tests.

## Four-Axis Decision

The final individual decision is identical for all six claims: verification
`symbolic_verified`, review `accepted`, compatibility
`compatible_extension`, and epistemic `active`. For C-GOT-003 and C-GOT-005,
Lean is a formal corroborating subclaim rather than the claim-level
verification label. C-GOT-006 additionally has `category: synthesized` and
`layer: core`.

## Promotion Transaction

Promotion adds the canonical module/tests and formal file, registers all six
claims with individual evidence and assumptions, pins the resulting accepted
set in the next release, moves P248 to immutable campaigns, renders generated
docs and accepted memory, and runs the full validation boundary selected by
the impact record.

## Correction Check

The bounded correction check inspected only the two repaired C-GOT-006
statement/dependency edges, the corrected Lean composition, the three accepted
dependency scopes, and `lean-correction-output.txt`. It returned `ACCEPT` at
confidence `0.99`: both blockers are closed, the formal build passes, and
`commonMetricComposition` remains axiom-free. No second substantive pass was
performed.

## Result and Frontier

All six claims are accepted for promotion. The paper's microscopic substrate
action, universal material coupling, spin-two emergence, value of `G`,
preferred foliation, and empirical equivalence remain explicit exclusions,
not hidden debt. A future accepted massive coordinate Euler-Lagrange theorem
could strengthen C-GOT-006's worldline clause without changing this result.

## Cross-References

The review links `claim-draft.yaml`, `evidence/claim-results.yaml`,
`evidence/verifier-audit.yaml`, `evidence/impact-analysis.yaml`, attempts 0002
and 0003, `src/substrate_framework/optical_gothic.py`,
`tests/test_optical_gothic.py`, `formal/SubstrateFramework/OpticalGothic.lean`,
and accepted dependencies C-STG-001/C-WLN-001/C-WLN-002.
