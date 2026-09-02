---
description: Shell, wall, and bubble sectors of the exterior-degenerate M5 clock (issue #195)
author: Main
created: '2026-09-02T07:10:00+00:00'
updated: '2026-09-02T08:55:00+00:00'
tags:
- substrate-framework
- campaign-proposal
- m5
- wall-tension
- phase-decoupling
- fixed-charge-bag
category: proposals
confidence: established
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
| S1 | Stationary shell through the degenerate locus: exact slice reduction, EL system, first integral, exact static-wall nonexistence, exact Maxwell frequency, wall existence, tension functional | C-M5C scopes | L-S1-WALL | Decoupling, bag, stability | SHELL_EXISTS | Tested slice/profile/frequency branch | S2, S3 | established by attempt 0001 |
| S2 | Exact phase decoupling: orbit-fixed locus, zero mismatch energy, no volume term in (omega1-omega2)^2 | L-S1-WALL | L-S2-DECOUPLING | Bag, force law | DECOUPLING_ESTABLISHED | Tested two-region construction | S3, S4 | established by attempt 0001 |
| S3 | Fixed-charge bag: thin-wall selection R = 2 sigma/p, exact Legendre closure, limits | L-S1-WALL | L-S3-BAG | Stability, two-clock force | BAG_RELATIVE_EQUILIBRIUM_EXISTS | Thin-wall family/charge window | S4 | established by attempt 0001 |
| S4 | One wall mechanism, two arrangements; P249 and wall-picture limits compared | L-S2-DECOUPLING, L-S3-BAG | L-S4-COMPARISON | M5.32 equivalence | COMPARISON_ESTABLISHED | Comparison scope | S5 | established by attempt 0001 |
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
variational calculus, Noether theory, ODE monotonicity arguments, exact
interval-arithmetic certification, and declared SymPy methods. The M7 0.786
threshold, M5.32-derived numbers in either direction, fitted frequencies,
imposed tapers or masks, and any other model's spectra are excluded.

## Frozen Candidate Universe
The in-scope route families come from the issue's objective, the accepted
invariants, and the permitted imports; the boundary expands append-only when
failure or history teaches a new concept and never shrinks without user
approval.

| Route family | Why in scope | Known concepts | Coverage strategy |
| --- | --- | --- | --- |
| Static degenerate-shell profile | Issue candidate seed; tests whether pure statics can hold the shell | mechanical co-conservation and Derrick obstructions | exact obstruction derived (attempt 0001) |
| Rotating-frame Maxwell wall | The fixed-frequency stationary kink at degenerate depths | exact slice functional, Maxwell system | exact characterization + interval certification |
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
Maxwell frequency was derived before any comparison with it.

## Proposed Claim Delta
The campaign proposes claims C-M5W-001 (exact wall-slice reduction, EL system,
first integral, and real-psi/shear invariance), C-M5W-002 (static-shell
nonexistence with exact mechanism), C-M5W-003 (orbit-fixed locus and exact
phase decoupling with zero mismatch coefficient), C-M5W-004 (Maxwell
frequency, wall existence, tension functional), and C-M5W-005 (thin-wall bag
selection, Legendre closure, limits, and comparison), with dependencies
001 -> 002/004, 003 requiring 001+004, 005 requiring 001+003+004. Identifiers
P250 and C-M5W-001..C-M5W-005 were collision-searched against governance,
proposals, campaigns, memory, and src at freeze and are free. Consumers are
the generated docs, the claim graph, downstream two-clock frontier work
(frontier only), and the release manifest.

## Analytic Specification and Numerical Licenses
Every node's typed chain is frozen in `proposals/P250-shell-bubble-clock/proposal.yaml`.
The load-bearing exact objects are: the clock-slice functional with kinetic
metric diag(1/4,1/2,1/2,1/2) in (m,c,b,f); the mechanical first integral
T - V_omega = const; the orbit-fixed locus (psi = 0, tangent isotropic,
shear free) at which the generator norm Tr[A,S]^2/2 + |psi|^2 vanishes
identically; and the Maxwell condition min V_omega = 0. The named numerical
remainders were: exact isolation of the Maxwell frequency (closed by
Krawczyk interval certification after the elimination and PSLQ routes were
recorded and superseded/refuted) and the tension value sigma_0 (carried as
an exact variational functional; no claim needs its value). No production
numerics was opened; interval certification entered only with a fully frozen
design (rational box, 40-dps directed rounding, standard Krawczyk formula).

## Implementation and Oracle Plan
The implementation order is analytic. The canonical module
`src/substrate_framework/m5_wall_clock.py` exposes pure SymPy APIs for the
slice potential, EL residuals, first integral, locus orbit-invariance, slip
bound, Maxwell system, exact resultants, tension integrand, and thin-wall
bag laws; the test suite carries channel-local mutation checks (lock
strength, scalar potential) that must fail, directly closing the #190
defect class. Campaign verifiers run with PYTHONPATH=src; reusable logic is
imported from the package, never copied. Attempt artifacts are captured
under the proposal attempts directory on first execution.

## Attempts and Continuation
Attempts are append-only with active obligation, route verdict, failure
scope, licenses earned and missing, and the next materially different route.
A failed route is a branch point: method repair, representation change, or a
failure-generated candidate is executed in the same run.

| Attempt | Obligation | Route | Verdict and failure scope | Licenses earned/missing | Method repair | Representation change | Alternative concept | Routes remaining |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0001 | S1 | WALL_STATIC_PROFILE | refuted for the objective's purpose: no omega = 0 interior state (exact unique-zero decomposition); radial static shell killed by exact Derrick; planar static loop left as frontier | static-state nonexistence earned | none needed | rotating frame | WALL_ROTATING_MAXWELL | - |
| 0001 | S1 | WALL_ROTATING_MAXWELL | established: exact slice, EL, first integral, Maxwell system, Krawczyk-unique solution in rational box, Gershgorin-positive Hessian, exact witnesses 5/3 and 45/16 | Maxwell local-isolation + nondegeneracy earned; global SOS named as closing route | lex-Groebner and degree-108 Sturm superseded as computationally prohibitive; PSLQ shortcut refuted at 90 dps | m = 0 deep branch (failure-generated) | - | geometry ladder, dynamical shells |
| 0001 | S2 | SLIP_LOCUS_DECOUPLING | established exactly: orbit-fixed locus, all density layers orbit-invariant, zero-cost phase slip, mismatch coefficient exactly zero, no volume term | all S2 licenses earned | - | - | - | curvature corrections |
| 0001 | S3 | BAG_THIN_WALL | established exactly in the reduced family: R = 2 sigma/p, dE/dQ = omega envelope, iota_int envelope, wall-to-infinity limit | all S3 thin-wall licenses earned | fixed-omega vs fixed-Q ensemble bookkeeping repaired en route | - | - | thick-wall regime |
| 0001 | S4 | comparison | established: one orbit-fixed-locus wall mechanism, two arrangements, exact limits | all S4 licenses earned | - | - | - | - |
| 0002 | G1 (S1 value layer) | DIRECT_LONG_BOX_COLLOCATION | blocked as production route: Newton falls into multi-kink trains from every direct long-box guess (basin structure, documented probes) | landscape license earned: V_w >= 0 on the corridor, pass height ~0.309 | L-continuation from short boxes | warm-started Dirichlet continuation | string/NEB paths (string: shelf-stalled, kept as L=2 guess and order-check; NEB: diverged, recorded) | h-ladder, BC cross-check, budget items |
| 0002 | G1 (S1 value layer) | L_CONTINUATION_KINK | established (numeric, SIGMA0_SLICE_VALUE): sigma_0 = 0.72929841787(21) at omega*^2; profile constructed, monotone, T = V_w pointwise; routes agree 1.6e-11; BC agreement 1.1e-13; budget total 2.0e-7 PASS | tension-value and profile-existence licenses earned; claim candidate C-M5W-006 named | - | - | - | G2 bag, G3 SOS, G4 stability |
| 0003 | G2 (S3 value layer) | BAG_DELTA_CONTINUATION | established (numeric, BAG_EXISTS slice scope): R(w) family over delta = w^2 - w_*^2 in [1e-3, 7e-3], R: 1217.4 -> 171.1; chi = Rp/(2 sigma_0) -> 1 (0.99956 at delta=0.001, linear beyond-thin-wall deviation ~ -1.75 delta); envelope dE/dQ = omega verified to 2e-4; virial/monotone identities hold; direct singular-term full-domain solves stall (residual metric col_res/h spiral + trivial basin) — domain-reduced clamped windows with amplitude-matched distances (m*dist = 9.2) are the working formulation | bag-family, thin-wall-law, and envelope-value licenses earned; thick-wall continuation (delta > 7e-3) named frontier | Lommel S-formulation replaced by domain reduction | - | - | G3 SOS, G4 stability |
## Debt Ledger
The ledger is empty. Declared hypotheses (thin-wall window, planar
exactness, spherical leading order, local scope of the Maxwell
certification), honest exclusions (stability, thick-wall regime, two-clock
force), and the open frontier are not debt.

## Review and Promotion Plan
One claim-level review covering the frozen C-M5W-001..005 transaction by a
distinct reviewer agent; package extraction already in
`src/substrate_framework/m5_wall_clock.py` with tests; registry and release
update (v0.172.0); `scripts/render_docs.py`; accepted-memory synchronization;
one impact-selected validation receipt reused by author, reviewer, and
merger. Reserve `refuted` for exact contradictions; evidence-attachment
roles are declared per artifact.

## Scientific Exhaustion Certificate
This section stays active and incomplete; positive success is the campaign's
terminal path. The open frontier classes (thick-wall bags, non-spherical
shells, dynamical shells, global SOS certificate) are continuations, not
exhaustion evidence.

| Candidate class | Routes and equivalence partition | Method repair | Representation change | Alternative concepts | Coverage evidence | Terminal verdict |
| --- | --- | --- | --- | --- | --- | --- |

- Independent candidate-generation artifact:
- Infinite-class coverage/no-go artifacts:
- Routes remaining:
- Exhaustion review:

## Continuation State (owner-directed extension, 2026-09-02)

The owner accepted the promoted claims and directed that the concrete-object
gaps now be cleared on the SAME branch and PR #196: numerics are authorized
because the concept is sorted, and they must stay concept-first (analytic
receipt before any solver). The full technical brief with the exact frozen
state, gap list, priorities, and paid-for gotchas lives in
`proposals/P250-shell-bubble-clock/continuation.md` — read it first on
resume. Gap status: G1 DONE (attempt 0002), G2 DONE (attempt 0003:
R(w) family over delta in [1e-3, 7e-3], chi -> 1 thin-wall law verified,
envelope dE/dQ = omega at 2e-4, bag profiles exported; thick-wall
continuation delta > 7e-3 named frontier). Remaining: G3 global SOS
certificate closing omega_c^2 = omega_*^2 (corridor nonnegative, interior
min 0.015 at m=0.1, pass height 0.309 — supportive); G4 bag stability.
Claim candidates: C-M5W-006 (numeric tension + profile, G1),
C-M5W-007 (bag family + law + envelope, G2) — formalize after owner
checkpoint.
Attempt records continue under proposals/P250-shell-bubble-clock/attempts/
(0002 = G1). The terminal PR stays open (Fixes #195) and is extended by
this work per owner authority; distinct merger still required; no
issue/discussion comments until PR approval.

## Done Gate
The campaign PR gate opened on positive success (PR #196, Fixes #195) and is
extended by owner direction with the G1-G4 continuation. The branch banks
each gap as a checkpoint; the PR is updated in place; distinct merger only.
