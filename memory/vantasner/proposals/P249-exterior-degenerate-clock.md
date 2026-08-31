---
description: Construct a licensed exterior-degenerate SO(2) M5 clock with intrinsic fixed-J width
author: Main
created: '2026-08-31T10:55:00+02:00'
updated: '2026-08-31T16:02:00+02:00'
tags:
- substrate-framework
- campaign-proposal
- m5
- fixed-j-clock
category: proposals
confidence: established
status: archived
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
| O3 | Complete fixed-J relative equilibrium with field, split, generator, inertia varied together and dE/dJ=omega | L-O1-CLOCK | L-O3-RELATIVE-EQUILIBRIUM | Isolation or stability | RELATIVE_EQUILIBRIUM_EXISTS | Tested branch/representation | O4 | established by full-field attempt 0006 |
| O4 | Finite limiting energy/inertia, nonzero omega, localized cumulative J, intrinsic width, natural tail | L-O3-RELATIVE-EQUILIBRIUM | L-O4-ISOLATION | Full stability or two-body law | ISOLATED_CLOCK_EXISTS | Tested branch/function space | O5 | established by full-field attempt 0006 |
| O5 | Complete constrained spectrum and no radiative growth in the declared class | L-O2-ASYMPTOTIC, L-O4-ISOLATION | L-O5-STABILITY | Particle/gravity/two-body interpretation | STABLE_ISOLATED_CLOCK | Represented physical tangent/dynamics | O6 | established by full-field attempt 0006 |
| O6 | Canonical APIs/tests, individual reviews, consumer replay, accepted claims, release, rendering, memory, no debt | L-O1 through L-O5 | L-O6-PROMOTION | Excluded phenomenology | CAMPAIGN_SUCCESS | Frozen C-M5C transaction | Terminal success | active |

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
| Benci--Fortunato, arXiv:0903.3508 | open primary mathematical source; theorem and hypotheses read in PDF | Exact nonlinear Klein--Gordon hylomorphic-soliton existence/stability gate used only for the auxiliary scalar seed |
| Benci--Fortunato, arXiv:1212.3236 | open primary mathematical source; abstract theorem and stability proof read in PDF | Translation-compact full-field hylomorphic solitons under EC-0 through EC-3 and strict binding |

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

Before attempt 0004, the auxiliary seed is further frozen as a complex Lorentz
scalar with standard positive kinetic term and amplitude potential
`W(f)=3 f^2-6 f^4+4 f^6`. Its global phase has period `2*pi`, its vacuum is
`f=0`, and its fixed-charge existence, localization, and stability conditions
are tested exactly. This is deliberately only a seed: it cannot satisfy O3 as
a decoupled scalar. The same declared action must still provide a covariant
coupling that dynamically lifts the M5 tangent degeneracy, and the complete
cross-sector first and second variations must close.

Before attempt 0005, the same-action bridge is frozen. A charge-one complex
Lorentz scalar and the weight-two tangent-traceless M5 tensor transform under
one diagonal `S1`. A constrained covariantly constant oriented clock-plane
spurion `(A,H)` defines `K=(AH-HA)/2` and the equivariant tensor
`Q(psi)=Re(psi^2)H+Im(psi^2)K`. The reproduced P239 density is completed by a
positive canonical kinetic term for the clock tensor and the positive lock
`6||B-Q(psi)||^2`, together with the attempt-0004 scalar action. The uniaxial
`B=psi=0` sector is unchanged; the clock-biaxial sector is explicitly a
challenged/completed sector. No frequency or comparator selects the coefficient.

The first source-level variation of attempt 0005 caught a preregistration
defect before candidate evaluation: adding canonical `B` kinetics alone would
also retain P239's clock-plane curvature inertia
`64 int b^2|grad b|^2`, so `I=int(f^2+4b^2)` would be incomplete. The recorded
repair replaces only the `A`-directed orthogonal projection of the curvature
two-form by the positive canonical `B` kinetic term and retains the orthogonal
curvature sector. This leaves static energies unchanged and makes the reduced
Legendre map exactly the frozen two-field one; it is an explicit minimal
kinetic replacement, not a silent omission.

The first coupled exact pass then exposed a theorem-bridge gap: the
attempt-0004 sextic lies at the three-dimensional `H^1` endpoint, while the
primary source's scalar truncation argument does not automatically establish
compactness for the coupled phase-locked functional. Before any verdict, the
coupled route alone is repaired to the subcritical `C^2` radial potential
`W(f)=3f^2-4f^3+2f^4`, `f=|psi|`. It has the exact positive form
`f^2[2(f-1)^2+1]`; all interactions are strictly subcritical. The earlier
passing tally remains preserved but cannot earn a license.

The subsequent O5 audit revoked attempt 0005's provisional O3/O4 promotion.
At the reduced background `S=diag(1,b,-b)`, the exact first derivatives in the
previously frozen director and tangent-trace directions are `8b^2` and
`-6b^2`. Thus the reduced solution is not a critical point of the full M5
field space. This is a route refutation, not a scalar or phase-lock refutation.
Attempt 0006 changes representation before any spectrum: all six spatial
symmetric M5 components become canonical physical variables, the exact M5
static potential is retained, and only typed auxiliary-frame/timelike
constraints are quotiented.

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
| 0004 | O3 | C_AUXILIARY_GAUGED/canonical-amplitude-Q-ball-seed | Blocked at the missing coupled M5 bridge; scalar seed established exactly | none; O3 active | Primary-source theorem audit replaces scalar shooting | Same-action covariant M5 core coupling | covariant core-degeneracy bridge | coupled C, noncommuting A, B, D, E |
| 0005 | O3/O4 | C_AUXILIARY_GAUGED/phase-locked-M5-core-bridge | Reduced binding established, later refuted as full critical point by omitted `8b^2,-6b^2` variations | O3/O4 revoked | Full-field variation audit | All six spatial components must vary | full canonical spatial completion | full E route plus B/D alternatives |
| 0006 | O2-O5 | E_M532_DYNAMICS/full-spatial-canonical-phase-locked-completion | Established in the complete declared physical tensor space | L-O1 through L-O5 earned | Full translation-compact theorem and exact spectrum signs | Canonical all-spatial M5 tensor | O6 promotion | canonicalization/review/release |

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

Attempt 0004 established the strongest analytic auxiliary-amplitude seed
without running a solver. For the complex Lorentz scalar potential
`W(f)=3f^2-6f^4+4f^6`, exact square completion gives positivity, the asymptotic
mass is `m^2=6`, and `inf 2W/f^2=3/2`. The Benci--Fortunato nonlinear
Klein--Gordon hypotheses hold exactly, including a negative nonlinear
remainder at `f^2=3/4` and the alternative growth condition, so a nonempty open
charge set of orbitally stable hylomorphic standing waves exists in the scalar
phase space. `Omega=2` is an exact interior frequency with tail
`exp(-sqrt(2)r)/r`. The complete fixed-charge Hessian contains
`Omega^2(delta N)^2/N`. This does not earn O3: the scalar is not yet coupled to
an M5 core, and the missing construction is now precisely a bounded covariant
same-action core-degeneracy bridge with complete cross-sector variations.

Attempt 0005 closes that bridge analytically. A constrained covariant clock
plane defines one diagonal `S1` acting with weights one and two on `psi` and the
M5 tangent tensor `B`. The selected action replaces only the old clock-plane
curvature square by positive canonical `B` kinetics and adds the positive lock
`6||B-Q(psi)||^2`. After repairing the coupled scalar potential to the
subcritical `W(f)=3f^2-4f^3+2f^4`, the exact reduced potential is positive and
all compact interactions have degree below six. Its split plateau has
`J/K=37/16`, strictly below the action-derived charged edge `9/2`, whereas the
entire unsplit sector has ratio at least `5`; more directly, the `b` equation
has source `-12f^2`, so no nontrivial solution can remain unsplit. The abstract
variational theorem therefore gives an open charge set of fixed-charge
minimizers on `R^3`, with finite volume inertia, `dE/dQ=omega`, positive natural
tails, and exact Derrick closure. O1-O4 are established for this selected
action. O5 remains active on the complete orthogonal constrained tangent and
neutral massless radiation channels.

That last attempt-0005 conclusion is superseded by the append-only attempt-0006
correction: the reduced core was not a full critical point. Attempt 0006 varies
all six spatial tensor components with a positive canonical Cartan metric and
retains the exact globally nonnegative M5 static potential. Its aligned
exterior pencil has generalized mass squares `(10,7,19,1,1,19)` plus scalar
mass square `6`; the diagonal charged edge is `19/4`. The full split plateau
has `J/K=39/16`, while every unsplit state is at least `5`. The abstract
translation-compact hylomorphic theorem applies because the local subcritical
energy satisfies EC-0 through EC-3 and the strict binding condition. It yields
a stable full-space constrained minimizing set with finite inertia, intrinsic
open-space scale, positive natural tails, nonnegative complete fixed-charge
Hessian, and subthreshold active harmonics. O1-O5 are now established for the
selected auxiliary-frame quotient; O6 is active.

Independent attempt-0007 review revoked that provisional O2--O5 conclusion.
For a full symmetric `S`, the generator norm is
`Tr([A,S]^T[A,S])/2=u^2+v^2+4p^2+4q^2`; attempt 0006 omitted the weight-one
shear pair `(u,v)`. At its provisional axis strength their mass square was one,
so `39/16` lay above the true edge and the passing checks were false green.

Attempt 0008 performs the finite exact repair without a numerical proxy. The
axis-lock coefficient is the unit value, the full Noether inertia includes the
shears, and the six tensor mass squares become `(10,10,22,4,4,22)`. Including
the scalar, the complete charged edge is `min(6,4,22/4)=4`. The split plateau
now has potential `45/64`, inertia `1/2`, and squared hylenic ratio `45/16`,
strictly below the edge by `19/16`; every `B=0` state remains at ratio at least
four. The theorem bridge is now on canonical `H1 x L2` field/velocity phase
space with the full Hamiltonian and bilinear Noether charge. Exact correction
checks report 28 identities and 31 targeted tests passing. The independent
bounded correction review accepted all four claims and O6 then completed.

## Debt Ledger
There is no in-boundary debt at freeze. Open scientific routes are frontier,
not debt. Any hidden fitted coefficient, undeclared chart restriction, action
term omitted before variation, unsupported promotion promise, or broken changed
consumer becomes debt and must be discharged before O6.

## Review and Promotion Record
One independent substantive review covered each C-M5C claim at the frozen
transaction boundary. It requested changes on all four claims, and the single
bounded correction check accepted all four after F1--F5 closed. Reusable logic
is canonical, the registry and release are pinned, generated docs and memory
are synchronized, affected consumers replayed, and the impact-selected receipt
is recorded in the immutable campaign.

## Scientific Exhaustion Certificate
Not applicable: the positive campaign success contract passed. Unselected
candidate families remain historical alternatives rather than obligations
because one full-scope candidate established every O1--O6 conjunct.

## Done Gate
The campaign ends positively only when O1-O6 pass and the promoted statement
has no in-boundary debt. Otherwise it continues until the full scientific-
exhaustion certificate passes. No rung, utility, subclaim, clean checkpoint, or
diminishing-return point opens a PR. The terminal PR links exactly issue #189;
it uses `Fixes #189` only on positive success and `Advances #189` on certified
exhaustion.

## Promotion Outcome

P249 completed positively in accepted release `v0.171.0`. The independent
review first caught the omitted charged shear doublet and rejected all four
draft claims; the bounded exact repair and correction check then accepted
`C-M5C-001` through `C-M5C-004`. The canonical result uses the complete
generator, exact charged edge four, strict split binding margin `19/16`, and
the phase-space hylomorphic minimizer theorem. The impact-selected final replay
passed 28 campaign checks, 31 direct and adjacent tests, every fixed repository
check, registry/release/generated documentation closure, and all 1,028 durable
memory files. No production numerical solver, discrete Hessian, or time
evolution was required.
