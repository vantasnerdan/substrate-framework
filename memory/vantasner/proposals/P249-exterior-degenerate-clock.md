---
description: Construct a licensed exterior-degenerate SO(2) M5 clock with intrinsic fixed-J width
author: Main
created: '2026-08-31T10:55:00+02:00'
updated: '2026-08-31T13:19:02+02:00'
tags:
- substrate-framework
- campaign-proposal
- m5
- fixed-j-clock
category: proposals
confidence: exploratory
status: active
---

## Question and Positive Deliverable
This campaign asks whether one explicitly declared local Lorentz-covariant M5
action can carry a smooth globally normalized compact SO(2) clock action fixing
the exterior vacuum pointwise and, under the complete fixed-J equations, a
finite-energy stationary open-space clock with finite positive inertia,
nonzero frequency, dynamically selected box-independent core width, and a
licensed stability statement. Positive completion is the O1-O6 object in issue
#189, including reusable implementation, individual claim review, promotion,
release, generated documentation, synchronized memory, and an empty in-scope
debt ledger. A failed route or finite-rung obstruction is not completion.

## Obligation Graph and Closure Map
The campaign freezes the following dependency graph; a passing node activates
the next unsatisfied parent obligation, while a failed route leaves that node
active and changes method, representation, or candidate.

| Node | Positive intent | Requires | Pass licenses | Does not license | Maximum verdict | Failure scope | Unlocks | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| O1 | Global normalized compact SO(2), pointwise exterior-vacuum fixation, exact invariance, Noether charge, regular Legendre map | Declared action/field space/charts | L-O1-CLOCK | Mode, solution, localization, stability | ESTABLISHED_SYMMETRY | Tested action/topology/generator | O2, O3 | established by attempt 0001 |
| O2 | Full action-specific principal symbol, kinetic rank/sign, positive physical subspace, essential/radiation thresholds | L-O1-CLOCK | L-O2-ASYMPTOTIC | Nonlinear clock or stability | ASYMPTOTIC_DYNAMICS_ESTABLISHED | Tested action/linear representation | O5 | established by attempt 0002 |
| O3 | Complete fixed-J relative equilibrium with field, split, generator, inertia varied together and dE/dJ=omega | L-O1-CLOCK | L-O3-RELATIVE-EQUILIBRIUM | Isolation or stability | RELATIVE_EQUILIBRIUM_EXISTS | Tested branch/representation | O4 | active |
| O4 | Finite limiting energy/inertia, nonzero omega, localized cumulative J, intrinsic width, natural tail | L-O3-RELATIVE-EQUILIBRIUM | L-O4-ISOLATION | Full stability or two-body law | ISOLATED_CLOCK_EXISTS | Tested branch/function space | O5 | pending |
| O5 | Complete constrained spectrum and no radiative growth in the declared class | L-O2-ASYMPTOTIC, L-O4-ISOLATION | L-O5-STABILITY | Particle/gravity/two-body interpretation | STABLE_ISOLATED_CLOCK | Represented physical tangent/dynamics | O6 | pending |
| O6 | Canonical APIs/tests, individual reviews, consumer replay, accepted claims, release, rendering, memory, no debt | L-O1 through L-O5 | L-O6-PROMOTION | Excluded phenomenology | CAMPAIGN_SUCCESS | Frozen C-M5C transaction | Terminal success | pending |

## Base Release and Provenance
The accepted base is release `v0.170.0` at
`substrate-framework@443f1bd6edf0c1997f42b05c3b6719142e7a11c1`. The campaign
branch starts from process-contract commit
`4071e5d469462b8c77612103d589a764fcc7d8b8`. The authority boundary is the
accepted release and registry, not the later campaign branch. Proposal/source
hashes are frozen in `proposals/P249-exterior-degenerate-clock/proposal.yaml`.

## Source Inventory and Access Gate
All load-bearing sources required to start O1 are accessible; every scientific
fact reused from memory will be verified against these sources.

| Source | Access status | Extracted claims |
| --- | --- | --- |
| `governance/releases/current.yaml`, `governance/claims.yaml` | in-hand repository source | v0.170.0 boundary and exact C-M5S scopes |
| Canonical M5 action, kinetic-axis, stationary-field, fluctuation-spectrum, and numerical modules | in-hand repository source | Starting APIs named by issue #189; bodies remain unopened until schema freeze passes |
| P239, P240, P244, P247 proposal/evidence surfaces | in-hand repository source | Conditional construction, fixed-J mechanism, spectrum machinery, route mechanisms to reproduce before reuse |
| Discussion #186 | open GitHub URL, opening post plus current comments/replies read via API | Next independent step, prior-route scope, M5.32 provenance caveat |
| Issue #189 | open GitHub URL | Positive objective, O1-O6 graph, candidate universe, gate pedigree, exclusions |

## Invariants, Conventions, and Allowed Imports
The campaign preserves local Lorentz/index validity, a fixed normalized compact
period, pointwise exterior-vacuum fixation, bounded Hamiltonian, regular
Legendre map, and exact or explicitly challenged static 3x3 recovery. It may
import accepted C-M5S-001 through C-M5S-020 only in their exact scopes and may
reuse P239/P240/P244/P247 machinery only after source reproduction/audit. The
M7 value 0.786, fitted clock frequency, imposed taper/wall/mask localization,
and static-profile post-processing are explicitly excluded inputs.

## Frozen Candidate Universe
The original and append-only-expandable candidate universe comes directly from
issue #189 and may not shrink without user approval.

| Route family | Why in scope | Known concepts | Coverage strategy |
| --- | --- | --- | --- |
| Exterior-degenerate spectral clock | Primary requested mechanism | Full-field SO(2) conjugation fixing equal-eigenvalue vacuum | Exact centralizer/orbit/invariance analysis, then complete fixed-J calculus |
| Bundle-aware ring clock | Director topology may obstruct a global oriented generator | n~-n charts, transition functions, holonomy | Explicit descent and global group-action proof |
| Auxiliary/gauged amplitude clock | A vacuum-vanishing amplitude/connection can make phase invisible outside | One compact phase plus amplitude or compensator | Type all fields/gauge signs and compare assumption cost |
| Non-spherical/core-selected continuation | P247 leaves fixed-axis and core-selected branches open | Fixed-axis rotor, localized defect core without imposed support | Full angular representation and dynamical width equation |
| M5.32 action-specific dynamics | The certified action requires its own quadratic/asymptotic analysis | Direct pencil, nonlinear mode, minimal covariant kinetic completion | Reproduce action, derive stabilizer/rank; continue if direct quadratic route degenerates |

## Candidate Preregistration
Candidate selection is structural; numerical closeness cannot select a clock.

| Candidate | Description | Assumptions | Parameters | Framework-fit prediction | Decisive test |
| --- | --- | --- | --- | --- | --- |
| A_EXTERIOR_DEGENERATE | Global spectral SO(2) with equal exterior eigenvalues and defect-lifted core split | Smooth full-field spectral action/descent exists | None beyond declared action coefficients | Best economy and direct fit to issue mechanism | Exact global action, invariance, Noether and regular-inertia proof |
| B_BUNDLE_RING | Compact action on disclination-ring bundle | Nontrivial chart/holonomy data allowed | No fitted clock constants | Natural topology if global orientation fails | Explicit transition-compatible action and conserved charge |
| C_AUXILIARY_GAUGED | Vacuum-vanishing amplitude or compensating connection | Additional typed field/gauge structure | Every new coefficient declared | Higher assumption cost but possible regular clock | Positive Hamiltonian, no ghost, regular quotient and Legendre map |
| D_NONSYMMETRIC_CORE | Fixed-axis or nonspherical core-selected rotor | Symmetry reduction widened | No imposed support scale | Can lift P247 width quasi-modulus | Complete angular fixed-J equations select width |
| E_M532_DYNAMICS | Certified M5.32 action or minimal preserving kinetic completion | Pinned external action reproduction | No borrowed M7 threshold | Direct provenance but possibly degenerate quadratic sector | Exact action-specific stabilizer, kinetic rank and principal symbol |

## Selection Criteria and Blinding
The frozen order is exact Lorentz/index validity; global symmetry/topology and
fixed period; Hamiltonian boundedness and regular Legendre map; exact or
explicitly challenged static recovery; natural exterior and dynamically
selected width; assumption/field/parameter economy; correct limits and accepted
sector fit; analytic tractability; numerical conditioning; downstream reach.
Historical values exposed by issue #189 are provenance only. Published clock
frequencies, de Broglie normalization, and M7's 0.786 remain outside derivation
and selection until equations, conventions, exact variations, observables,
tolerances, and budgets freeze.

The requester additionally directs the campaign to stay within algebra and
calculus for as long as possible. O1 is exact-only. O2 begins with exact
second-variation and polynomial/pencil analysis. O3's complete fixed-J first
and second variations must be derived symbolically before a numerical solver is
eligible. Numerical work is reserved for existence, localization, or spectral
questions that remain analytically unresolved after those gates.

## Proposed Claim Delta
The collision search across registry, proposals, campaigns, and memory found no
P249 or C-M5C identifiers on 2026-08-31. P249 provisionally reserves
C-M5C-001 for compact symmetry/Noether structure, C-M5C-002 for action-specific
asymptotic dynamics, C-M5C-003 for the isolated fixed-J relative equilibrium and
intrinsic width, and C-M5C-004 for the declared physical stability class. These
identifiers remain reserved even if a route fails. Dependencies follow the
O1-O5 graph. Consumers are the canonical M5 action/clock/spectrum APIs, their
tests, downstream isolated-clock/two-clock campaigns, generated claims/releases,
and durable memory.

## Mathematical and Numerical Licenses
Each node's complete object-to-verdict chain is machine-frozen in the proposal
manifest. Production numerics are forbidden until the upstream license resolves.
O1 consumes no numerical evidence for its load-bearing conclusion. O2 must
differentiate the full unreduced action before background restrictions and type
the physical positive-kinetic subspace. O3 must vary the generator, inertia,
constraints, and rank-one Legendre term. O4 requires unsupported open-space
tails and crossed domain/mesh designs. O5 binds the small-ratio prescriptions:
full operator, zero-mode gauge, lambda_min/lambda_2, eigenpair residual,
observed-order extrapolation, independent boundary/discretization treatment,
non-soft monitor, deterministic evaluator-floor measurement, and epsilon-jitter
sign test.

## Implementation and Oracle Plan
The implementation order is analytic. Attempt 0001 reconstructs the exact
P239/P240 action and reduced clock map, types fields and measures, computes the
global exterior-degenerate centralizer action, proves or disproves exact action
invariance, derives the Noether current and Legendre form, and records every
chart/topology obstruction. SymPy is regression evidence for exact matrix
identities; hand-derived indexed calculus remains visible in the evidence. A
numerical sample cannot repair a failed global group action.

After O1, O2 expands the complete second variation and derives the asymptotic
kinetic/principal-symbol polynomial exactly, including rank, sign, constraints,
zero modes, and channel thresholds. After that, O3 derives the full fixed-J
Routhian and its first and second variations before any ansatz restriction.
Only equations irreducible by exact or qualitative calculus become numerical
BVP/PDE/continuation targets. Those later solvers must use canonical numerics,
record precision/thread/library/hardware conditions, use compensated reductions,
capture stdout on first execution, cross h and L independently, and carry the
itemized small-ratio budget. No numerical rerun counts as independent evidence
for an exact identity.

Reusable exact transformations and variational operators will live under
`src/substrate_framework/`; campaign scripts remain thin. Compatibility
preflight rejects direct/imported/dynamic `np.trapz` and eager nested fallbacks.
Targeted tests, mutation of signs/normalizations/generators, exact limits, and
impact-bounded consumer replay are frozen before promotion.

## Attempts and Continuation
Attempts are append-only under
`proposals/P249-exterior-degenerate-clock/attempts/`. Each ends with a
route-scoped established, refuted-with-mechanism, or blocked-with-missing-
construction verdict plus computed predicate, implied proposition, consumed and
earned licenses, parent effect, and what it cannot decide.

| Attempt | Obligation | Route | Verdict and failure scope | Licenses earned/missing | Method repair | Representation change | Alternative concept | Routes remaining |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0001 | O1 | A_EXTERIOR_DEGENERATE/commuting-radial-split | Refuted for fixed-width minimum by exact amplitude-width collapse; exact global SO(2) survives | L-O1-CLOCK earned | Expanded factored coefficient before extraction | Noncommuting spatial core required | A_EXTERIOR_DEGENERATE/noncommuting-core | noncommuting A plus B-E |
| 0002 | O2 | A_EXTERIOR_DEGENERATE/exact-asymptotic-pencil | Established as stated; clock split is nonlinear/singular at quadratic order | L-O2-ASYMPTOTIC earned | None | Exact 10-channel rational pencil | Nonlinear fixed-J or positive kinetic completion | nonlinear A, completion C, B, D, E |
| 0003 | O3 | A_EXTERIOR_DEGENERATE/commuting-split-plus-axis-lock | Refuted by incompatible exact scale/amplitude virials | none; O3 active | Corrected collapse to tau=4+2d inside simple branch | Noncommuting curvature or volume inertia | canonical-amplitude Q-ball | noncommuting A, auxiliary C, B, D, E |

Attempt 0001 established the strongest exact positive result available at O1.
For `A_23=-1`, `A_32=1`, the completed P239 action admits the smooth global
action `M -> exp(qA) M exp(qA)^T`, with period `2*pi`. The exterior
`diag(-g,1,0,0)` is fixed pointwise, the full density is invariant, localization
of the group parameter gives the conserved current
`j^mu=<Pi^mu,AM+MA^T>`, and a smooth Gaussian core split has finite positive
`C=6*pi^(3/2)` and nondegenerate `J=2C omega`.

The attempt also closed one O3 representation exactly. For the commuting split
`diag(-g,1,d,-d)`, `C=32 int d^2|grad d|^2` and
`V=3d^2+4d^4`. With `u=d^2`, the fixed-J energy is
`int(3u+4u^2)+J^2/(32 int|grad u|^2)`. The family
`u=a f(x/R)`, `R=a^-1`, sends every term to zero while every finite nonzero-J
profile remains positive. The route collapses instead of selecting a width.
This is representation-scoped and generates the noncommuting-core route; it is
not a no-clock result.

Attempt 0002 established the complete action-specific exterior dynamics. In
the channel order `(00,01,02,03,11,12,13,22,23,33)`, the exact stiffness is
`diag(1,0,0,0,5/2,0,0,3/2,3,3/2)` and the projector-current kinetic metric is
`kappa*diag(0,1/9,1/16,1/16,0,0,0,0,0,0)`. The three positive kinetic
branches obey `omega^2=k^2`, so the essential/radiation threshold is zero.
The `(12,13)` director shears are the two-dimensional common inert kernel. The
clock split `(22-33,23)` is an SO(2) weight-two doublet with isotropic stiffness
`3(d^2+b^2)` and zero quadratic kinetic metric. It is therefore an intrinsically
nonlinear fixed-J candidate, not a linear mode with an importable mass gap.

Attempt 0003 established the complete reusable fixed-J calculus:
`delta E_J=delta E0-omega^2 delta C` and
`delta^2 E_J=delta^2 E0-omega^2 delta^2 C+2omega^2(delta C)^2/C`, including
the mandatory positive rank-one term. Its spatial virial retains every
`E4/R`, `E2 R`, `E0 R^3`, `C1 R`, and `C3 R^3` class. For the commuting split,
the scale and amplitude stationarity equations eliminate to
`-10 U2-8 U4=0`, so neither the base potential nor Candidate K's extra positive
quadratic axis lock can yield a nonzero stationary core. The corrected
`tau=4+2d`, `R=a^-2` collapse stays strictly within the simple timelike branch;
this supersedes only attempt 0001's collapse-path admissibility detail, not its
O1 symmetry result.

## Debt Ledger
There is no in-boundary debt at freeze. Open scientific routes are frontier,
not debt. Any hidden fitted coefficient, undeclared chart restriction, action
term omitted before variation, unsupported promotion promise, or broken changed
consumer becomes debt and must be discharged before O6.

## Review and Promotion Plan
One independent substantive review will be performed per proposed or changed
C-M5C claim at the frozen transaction boundary, with one bounded correction
check if needed. Evidence attachments receive exact-proof, corroborating,
regression, applicability, or provenance roles and do not generate duplicate
claim reviews. Accepted reusable logic moves into canonical package modules and
tests, the registry and new pinned release are updated, generated docs are
rendered rather than edited, accepted memory is synchronized, affected
consumers replay, and one impact-selected validation receipt is recorded.

## Scientific Exhaustion Certificate
This section remains incomplete because the positive objective and all route
families are active. Exhaustion would require a distinct non-author adversarial
candidate-generation pass from the frozen objective/invariants/source inventory,
reconciliation against every historical, external, preregistered, and failure-
generated route, equivalence partitions, coverage/no-go arguments for open
classes, complete continuation ladders, and no plausible route remaining.

| Candidate class | Routes and equivalence partition | Method repair | Representation change | Alternative concepts | Coverage evidence | Terminal verdict |
| --- | --- | --- | --- | --- | --- | --- |
| All issue-#189 route families | active | active | active | active | absent | nonterminal |

- Independent candidate-generation artifact: absent while campaign is active.
- Infinite-class coverage/no-go artifacts: absent while campaign is active.
- Routes remaining: A_EXTERIOR_DEGENERATE, B_BUNDLE_RING, C_AUXILIARY_GAUGED, D_NONSYMMETRIC_CORE, E_M532_DYNAMICS.
- Exhaustion review: not applicable while plausible routes remain.

## Done Gate
The campaign ends positively only when O1-O6 pass and the promoted statement
has no in-boundary debt. Otherwise it continues until the full scientific-
exhaustion certificate passes. No rung, utility, subclaim, clean checkpoint, or
diminishing-return point opens a PR. The terminal PR links exactly issue #189;
it uses `Fixes #189` only on positive success and `Advances #189` on certified
exhaustion.
