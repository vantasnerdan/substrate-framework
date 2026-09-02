---
description: Shell, wall, and bubble sectors of the exterior-degenerate M5 clock (issue #195)
author: Main
created: '2026-09-02T07:10:00+00:00'
updated: '2026-09-02T07:10:00+00:00'
tags:
- substrate-framework
- campaign-proposal
- m5
- wall-tension
- phase-decoupling
- fixed-charge-bag
category: proposals
confidence: exploratory
status: active
---

## Question and Positive Deliverable
This campaign constructs, on the accepted exterior-degenerate clock completion
C-M5C-001..004, the strongest validated account of finite degenerate
interfaces and the charge carriers they enclose: a stationary shell through
the degenerate locus with exact tension sigma_0, exact phase decoupling of two
clock regions across such a shell, a fixed-charge shell-bubble with
volume-vs-surface radius selection and dE/dQ = omega closure, and the
strongest honest comparison of the P249 exterior-degenerate picture against
the Discussion #186 isofrequency-wall picture. Positive completion is the
S1-S5 object of issue #195 with reusable implementation, individual claim
review, promotion, release, generated documentation, synchronized memory, and
an empty in-scope debt ledger. A failed route or finite-rung obstruction is
not completion.

## Obligation Graph and Closure Map
The campaign freezes the following dependency graph; a passing node activates
the next unsatisfied parent obligation, while a failed route leaves that node
active and changes method, representation, or candidate in the same run.

| Node | Positive intent | Requires | Pass licenses | Does not license | Maximum verdict | Failure scope | Unlocks | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| S1 | Stationary shell through the degenerate locus: exact slice reduction, EL system, first integral, exact static-wall nonexistence, exact Maxwell frequency, wall existence, tension functional and bounds | C-M5C scopes | L-S1-WALL | Decoupling, bag, stability | SHELL_EXISTS | Tested slice/profile/frequency branch | S2, S3 | active |
| S2 | Exact phase decoupling: orbit-fixed locus, zero mismatch energy, no volume term in (omega1-omega2)^2 | L-S1-WALL | L-S2-DECOUPLING | Bag, force law | DECOUPLING_ESTABLISHED | Tested two-region construction | S3, S4 | active |
| S3 | Fixed-charge bag: thin-wall selection R = 2 sigma/p, exact Legendre closure, limits | L-S1-WALL | L-S3-BAG | Stability, two-clock force | BAG_RELATIVE_EQUILIBRIUM_EXISTS | Thin-wall family/charge window | S4 | active |
| S4 | One wall mechanism, two arrangements; P249 and wall-picture limits compared | L-S2-DECOUPLING, L-S3-BAG | L-S4-COMPARISON | M5.32 equivalence | COMPARISON_ESTABLISHED | Comparison scope | S5 | active |
| S5 | APIs/tests, individual reviews, consumer replay, registry/release/docs/memory, empty debt | L-S4-COMPARISON | L-S5-PROMOTION | Excluded phenomenology | CAMPAIGN_SUCCESS | C-M5W transaction | Terminal success | active |

## Base Release and Provenance
The accepted base is release `v0.171.0` at
`substrate-framework@8b74d3ab5fdd307ba564291ff7bfa702db5d3330`, which equals
the hash frozen in issue #195. The campaign branch is
`research/P250-shell-bubble-clock`. Proposal/source hashes are frozen in
`proposals/P250-shell-bubble-clock/proposal.yaml`. The PR #190 defect class
(literal-returning oracle functions compared against literals) was verified
present at the frozen baseline in `m5_exterior_clock.py:300-311`; every
P250 oracle is instead derived symbolically from the canonical potential and
carries a channel-local mutation check that must fail.

## Source Inventory and Access Gate
All load-bearing sources are accessible in hand or via open URLs; every
scientific fact reused from memory is verified against these sources.

| Source | Access status | Extracted claims |
| --- | --- | --- |
| `governance/releases/current.yaml`, `governance/claims.yaml` | in-hand repository source | v0.171.0 boundary, exact C-M5C-001..004 statements |
| `src/substrate_framework/m5_exterior_clock.py`, `tests/test_m5_exterior_clock.py` | in-hand repository source | canonical action primitives, clock basis, lock strengths, witness data |
| `campaigns/P249-exterior-degenerate-clock/` | in-hand repository source | immutable route history, review conventions, companion check style |
| Issue #195 | open URL, body plus owner comment read at freeze | objective S1-S5, candidate universe, exclusions, executor protocol, no-comment instruction |

## Invariants, Conventions, and Allowed Imports
The campaign preserves the C-M5C-001 action as the only source of potentials,
spectra, and thresholds; the auxiliary frame stays a declared constrained
quotient; the exterior vacuum stays fixed pointwise; wall energetics come from
wall-localized observables, never from large-subtraction; area versus volume
discrimination is decided by exact locality, never by a box average; and every
claimed gate carries a mutation check that must fail. Allowed inputs are
C-M5C-001..004 in their exact accepted scopes plus standard linear algebra,
variational calculus, Noether theory, ODE monotonicity arguments, Sturm root
isolation, and declared SymPy methods. The M7 0.786 threshold, M5.32-derived
numbers in either direction, fitted frequencies, imposed tapers or masks, and
any other model's spectra are excluded.

## Frozen Candidate Universe
The in-scope route families come from the issue's objective, the accepted
invariants, and the permitted imports; the boundary expands append-only when
failure or history teaches a new concept and never shrinks without user
approval.

| Route family | Why in scope | Known concepts | Coverage strategy |
| --- | --- | --- | --- |
| Static degenerate-shell profile | Issue candidate seed; tests whether pure statics can hold the shell | mechanical co-conservation obstruction expected | derive exact obstruction or construct |
| Rotating-frame Maxwell wall | The fixed-frequency stationary kink at degenerate depths | exact slice functional, Maxwell frequency | exact characterization + Sturm isolation |
| Orbit-fixed-locus decoupling | The issue's S2 mechanism via the vanishing clock tangent | orbit-fixed locus, free phase slip | exact symmetry theorem |
| Fixed-charge thin-wall bag | The issue's S3 volume-surface selection | reduced (R, omega) functional | exact stationarity + Legendre closure |
| Geometry ladder | Planar through spherical arrangements | area scaling, curvature corrections | exact leading laws, curvature as frontier |
| Representation variants | Director/bundle charts through the locus | n ~ -n identification | justify orientation-jump bookkeeping |

## Candidate Preregistration
The mechanism is not being selected against an external rival; the fixed
theorem target is the wall sector of the declared action, so one complete
construction route is registered per node, with the static-shell route
explicitly registered as the competing concept expected to fall to an exact
obstruction, which is decisive evidence for the rotating-frame construction
rather than a manufactured rival.

| Candidate | Description | Assumptions | Parameters | Framework-fit prediction | Decisive test |
| --- | --- | --- | --- | --- | --- |
| WALL_STATIC_PROFILE | omega = 0 finite-energy shell | none beyond canonical action | none | falls to exact obstruction | exact co-conservation identity |
| WALL_ROTATING_MAXWELL | kink at the Maxwell frequency | rotating frame, fixed omega | omega (eliminated) | exists with exact tension functional | exact EL + first integral + existence |
| BAG_THIN_WALL | thin-wall fixed-Q bag | thin-wall window near Maxwell | Q via omega | radius selected exactly at leading order | volume-surface stationarity |
| SLIP_LOCUS_DECOUPLING | two ticking regions, one shell | orbit-fixed locus | none | zero mismatch exactly | energy bookkeeping identity |
| BUNDLE_DIRECTOR_CHART | director chart through locus | same action, chart change | none | same wall objects | chart-invariance of tension |

## Selection Criteria and Blinding
Structural criteria ordered: canonical-action-only derivation, exact analytic
tractability, assumption economy, correct limits to accepted C-M5C sectors,
downstream reach. No numerical agreement ever selects a construction. The
accepted witness ratio 45/16 is a provenance cross-check only; the exact
Maxwell frequency is derived before any comparison with it.

## Proposed Claim Delta
The campaign proposes claims C-M5W-001 (exact wall-slice reduction, EL system,
first integral, and real-psi/shear invariance), C-M5W-002 (static-shell
nonexistence with exact mechanism), C-M5W-003 (orbit-fixed locus and exact
phase decoupling with zero mismatch coefficient), C-M5W-004 (Maxwell
frequency, wall existence, tension functional with bounds), and C-M5W-005
(thin-wall bag selection, Legendre closure, and limits), with dependencies
001 -> 002/004, 003 requiring 001, 005 requiring 004, and the S4 comparison
riding on 003+005 inside the claim statements. Identifiers P250 and C-M5W-001..
C-M5W-005 were collision-searched against governance, proposals, campaigns,
memory, and src at freeze and are free. Consumers are the generated docs, the
claim graph, downstream two-clock frontier work (frontier only), and the
release manifest.

## Analytic Specification and Numerical Licenses
Every node's typed chain is frozen in `proposals/P250-shell-bubble-clock/proposal.yaml`.
The load-bearing exact objects are: the clock-slice functional with kinetic
metric diag(1/4,1/2,1/2,1/2) in (m,c,b,f); the mechanical first integral
sum K u'^2 - V_omega = const; the orbit-fixed locus (psi = 0, tangent
isotropic, shear free) at which the generator norm Tr[A,S]^2/2 + |psi|^2
vanishes identically; and the Maxwell condition min V_omega = 0. The only
candidates for a named numerical remainder are the tension value sigma_0 and
the exact algebraic isolation of the Maxwell frequency; production numerics
stays blocked until the exact polynomial elimination and Sturm isolation
route is tried and recorded, and any eventual tension evaluation enters
through the geodesic functional with an itemized error budget, never through
large-subtraction.

## Implementation and Oracle Plan
The implementation order is analytic. The canonical module
`src/substrate_framework/m5_wall_clock.py` will expose pure SymPy APIs for
the slice potential, EL system, first integral, locus orbit-invariance, slip
bookkeeping, Maxwell system with exact elimination and Sturm isolation, the
tension functional with exact bounds, and the thin-wall bag laws; every
public check carries a mutation that must fail (lock strength, W
coefficients, lock stiffness, inertia normalization), directly closing the
#190 defect class. Campaign verifiers run with PYTHONPATH=src; reusable
logic is imported from the package, never copied. Attempt artifacts are
captured under the proposal attempts directory on first execution.

## Attempts and Continuation
Attempts are append-only with active obligation, route verdict, failure
scope, licenses earned and missing, and the next materially different route.
A failed route is a branch point: method repair, representation change, or a
failure-generated candidate is executed in the same run.

| Attempt | Obligation | Route | Verdict and failure scope | Licenses earned/missing | Method repair | Representation change | Alternative concept | Routes remaining |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |

## Debt Ledger
The ledger is empty at freeze. Declared hypotheses (thin-wall window, planar
exactness, spherical leading order), honest exclusions (stability, thick-wall
regime, two-clock force), and the open frontier are not debt.

## Review and Promotion Plan
One claim-level review per proposed C-M5W claim with the frozen transaction;
package extraction into `src/substrate_framework/` with tests; registry and
release update; `scripts/render_docs.py`; accepted-memory synchronization;
one impact-selected validation receipt reused by author, reviewer, and
merger. Reserve `refuted` for exact contradictions; evidence-attachment
roles are declared per artifact.

## Scientific Exhaustion Certificate
This section stays active and incomplete unless positive success passes
first. To claim exhaustion the full candidate universe above would need every
route verdict, an independent adversarial candidate-generation artifact, and
coverage or no-go arguments for the open classes (thick-wall bags,
non-spherical shells, dynamical shells).

| Candidate class | Routes and equivalence partition | Method repair | Representation change | Alternative concepts | Coverage evidence | Terminal verdict |
| --- | --- | --- | --- | --- | --- | --- |

- Independent candidate-generation artifact:
- Infinite-class coverage/no-go artifacts:
- Routes remaining:
- Exhaustion review:

## Done Gate
The campaign PR gate opens only on complete positive success or the scientific
exhaustion certificate. Until then the campaign branch banks checkpoint
commits and the executor continues on the active obligation.
