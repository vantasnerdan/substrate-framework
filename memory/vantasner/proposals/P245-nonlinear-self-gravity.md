---
description: Construct and classify the nonlinear self-gravitating confined clock beyond the weak-field wall
author: Main
created: '2026-08-26T19:30:00+00:00'
updated: '2026-08-26T20:25:00+00:00'
tags:
- substrate-framework
- campaign-proposal
- nonlinear-gravity
- confined-clock
- small-ratio-numerics
category: proposals
confidence: established
status: archived
---

## Question and Positive Deliverable

This campaign asks whether the accepted confined-clock sector at the purely induced `B=0`, `xi=0` coupling has a regular horizonless nonlinear spherical completion or necessarily develops a trapped surface. The positive deliverable is an importable areal-gauge Einstein constraint for the typed clock source, a verified Misner–Sharp compactness and ADM-mass construction, and a nonlinear fate theorem at the accepted coupling: either a self-consistent regular branch or a refinement-stable trapped-surface crossing with its mechanism and branch boundary located. A failed nonlinear solver does not complete the campaign.

## Base Release and Provenance

The accepted base is release `v0.166.0` at framework commit `3e3a6133030b3041515d07cbeef74eeb3624b99b`; work starts on `campaign/p245-nonlinear-self-gravity`. Canonical issue #174 existed before implementation and PR submission. Accepted inputs read at source are C-M5S-003/004/007/010/011, C-GRV-002, C-STG-002 and qualified C-PDE-013; package inputs are `m5_covariant_action.py`, `m5_kinetic_axis.py`, `m5_stationary_fields.py`, `spherical_einstein_scalar.py`, `spherical_einstein_scalar_bvp.py`, `linearized_einstein.py`, and `numerics.py`. The source profile is the committed P240 R=12 order-16 family-S root only through C-M5S-004's accepted provenance role. P244/PR #173 is an unmerged separate transaction and is neither edited nor treated as accepted authority.

## Source Inventory and Access Gate

All load-bearing sources are in hand inside the pinned repository, so no external-access blocker exists. C-M5S-004 supplies `G_total=46.80699908016004`, classical source mass `54.70900884959007`, radius `12`, the spherical density construction and the declared exclusion of the global fixed-J term. C-M5S-007 supplies the exact induced-versus-bare dichotomy. C-M5S-011 supplies a separately typed in-model zero-point shift `72.58859645998888`; its de-boxed consumer attribution remains open and will be reported as a second source ledger rather than silently added. C-STG-002 and C-PDE-013 supply accepted areal-gauge conventions and numerical precedent, not a field map from M5 to sine-Gordon.

| Source | Access status | Extracted claims |
| --- | --- | --- |
| `governance/claims.yaml` C-M5S-003/004/007/010/011 | in hand | Coupling, source, wall, baseline dichotomy and radiative mass ledger |
| `m5_covariant_action.py`, `m5_kinetic_axis.py` | in hand | Covariant action primitives and exact fixed-J scaling |
| `m5_stationary_fields.py`, P240/P243 evidence | in hand | Biaxial profile and radial density/inertia construction |
| C-STG-002, C-PDE-013 and spherical Einstein modules | in hand | Areal-gauge equations, origin/exterior conditions and BVP precedent |

## Invariants, Conventions, and Allowed Imports

The campaign preserves the accepted mostly-plus sign, areal radius, positive-density mass orientation, fixed-J topology, boundary spectrum and action coefficients. The nonlinear radial metric function is `f(r)=1-2*G*m(r)/r`, with `m'(r)=4*pi*r^2*rho(r)` when `m` is the matter mass in C-M5S-004 units. The classical source convention and the optional C-M5S-011 radiative ledger remain separate. No equation of state, scalar proxy, or ADM mass is inferred from a shared symbol or numerical resemblance. Permitted imports are the accepted claims and pinned canonical modules named above; a new action/stress bridge must be derived and reviewed.

Small-number discipline is binding even though the headline compactness is large: nondimensionalize on `x=r/R` and `M`, accumulate the cumulative mass without subtracting large energies, pin BLAS threads, report the evaluator noise floor, refine profile evaluation and radial quadrature independently, bracket roots rather than infer them from a minimum sample, and require signal-to-error margins. Near a continuation fold, use pseudo-arclength rather than natural-parameter stepping; monitor the smallest Jacobian singular direction and constraint/virial residuals.

## Candidate Preregistration

Three mechanisms are registered before any new nonlinear profile or horizon value is evaluated. Candidate A is the frozen-source Einstein constraint: compute the exact necessary Misner–Sharp gate for the accepted C-M5S-004 density. Candidate B is the full covariant M5 continuation: derive curved-space stress and field equations, then co-solve profile and metric from `G=0` to the accepted coupling with Newton/pseudo-arclength continuation. Candidate C is a homothetic fixed-J branch: use the exact Derrick scale ledger to decide whether dilation alone can produce a horizonless stationary point.

| Candidate | Description | Assumptions | Parameters | Framework-fit prediction | Decisive test |
| --- | --- | --- | --- | --- | --- |
| A | Frozen-source Misner–Sharp constraint | Accepted spherical density held fixed | No fitted parameters | Exact necessary nonlinear gate, not a full matter solution | Refinement-stable sign-bracketed `f=0` or positive full-domain margin |
| B | Covariant M5 profile-plus-metric continuation | Spherical averaged biaxial fixed-J ansatz | Accepted action coefficients and continuation coordinate | Highest-fidelity clock completion | Reach accepted coupling with regular metric and residual closure, or locate a genuine branch endpoint |
| C | Homothetic fixed-J nonlinear branch | Accepted shape retained under scale change | Scale only | Isolates dilation escape without claiming arbitrary deformation | Stationary scale plus full-domain compactness classification |

## Selection Criteria and Blinding

Structural fit outranks numerical convenience. The ordered criteria are: exact accepted-action/stress/sign compatibility; no hidden equation of state or borrowed scalar map; positive nonlinear reach or a trapped-surface theorem; correct `G -> 0` Poisson and exterior Schwarzschild limits; parameter economy; clear separation of frozen, homothetic and co-solved scope; and residual/refinement/independent-route robustness.

The accepted C-M5S numbers were public and known before P245, so full comparator blinding is impossible and this exception is explicit. The new compactness profile, horizon radius, homothetic minimum, critical coupling and continuation values remain unopened until both contracts validate. A horizonless result requires `f(r)` to exceed one hundred times the composed absolute `f` error everywhere. Trapped-surface evidence requires a sign-bracketed `f=0` root stable under profile and radial refinement and at least ten local mesh spacings from either boundary. Failure of either band changes the method; thresholds do not move after data.

## Proposed Claim Delta

Provisional `C-M5S-012` is collision-free across the registry, campaigns and durable memory as of 2026-08-26. Its strongest intended scope is the nonlinear spherical fate of the accepted clock source and any derived critical-coupling or homothetic branch boundary. Dependencies are C-M5S-003/004/007 and, only for a separately reported source ledger, C-M5S-010/011; C-STG-002/C-PDE-013 are convention and numerical precedents rather than matter-model dependencies unless an explicit typed bridge is constructed. Direct consumers are C-M5S-005's applicability warning, C-M5S-007's nonlinear escape branch, issue #170's strong-coupling frontier, and future clock-pair dynamics.

## Implementation and Oracle Plan

Reusable equations and solvers belong under `src/substrate_framework/`; attempt scripts only load the accepted root and call them. The exact oracle derives the areal-gauge mass constraint, weak-field limit and exterior matching with SymPy or direct algebra. The numerical oracle uses SciPy cumulative integration/root bracketing and, for co-solved equations, adaptive collocation plus an independently written shooting or finite-volume route. Every run records float64 precision, mesh/domain, solver and tolerances, boundary data, status, residual norm, minimum `f`, compactness maximum, ADM/Misner–Sharp mass, error budget and runtime.

Refinements independently vary radial quadrature, profile evaluator resolution, domain and solver tolerance. Counterexamples flip the Einstein sign, remove source mass, rescale `G`, corrupt boundary matching and perturb one load-bearing profile coefficient; relevant gates must fail. The `G -> 0` solution must reduce to the accepted Poisson sign/monopole normalization. Candidate B starts from a controlled weak-coupling core and uses pseudo-arclength at a fold. GitNexus impact analysis precedes canonical-symbol edits and impact-bounded replay follows them.

## Attempts and Continuation

Attempts are immutable under `campaigns/P245-nonlinear-self-gravity/attempts/`. Attempt 0001 preserved an oracle-only abort after seven scientific gates passed: equivalent reference-scale formulas differed by `1.194e-12`, just beyond an unregistered `1e-12` absolute check. Attempt 0002 repaired only that roundoff oracle to `16*eps*max(1,maxC)` and passed 11/11 gates. The accepted frozen source has first `f=0` at `r=0.0934704`, maximum compactness `2Gm/r=893.81197`, minimum `f=-892.81197`, and matched vacuum exterior Schwarzschild radius `2GM=5121.5305`; the critical coupling for that profile is `G=0.0523678`, versus accepted `G=46.807`. Candidate C's convex-in-log-scale homothetic objective has global minimum maximum compactness `881.50538` at `R=13.02676`, with critical `G=0.0530989`, so dilation misses horizonlessness by a factor above 881. Attempt 0003 independently re-evaluated analytic Chebyshev and angular derivatives in NumPy, integrated component masses with PCHIP plus DOP853, and bracketed the horizon with Brent; all seven gates passed, mass agreed within `6.9e-7` relative, first-horizon radius within `8.0e-4`, maximum compactness within `2.1e-6`, and homothetic minimum within `1.9e-6`. Independent reviewer `ReviewCMM5S012` accepted the claim at confidence `0.97`; C-M5S-012 and release v0.167.0 are the promoted state.

## Debt Ledger

The frozen C-M5S-012 transaction has no hidden equation of state or scalar proxy: the spherical Hamiltonian constraint is necessary and closes before lapse or pressure data can restore `f>0`. It states only the accepted certified source and its exact curvature-plus-potential homothetic dilations. The global fixed-J term remains excluded under C-M5S-004's declared convention; because its physical energy is nonnegative, this omission makes the trapped conclusion conservative but supplies no local stress claim. Arbitrary covariant M5 profile relaxation would change the certified source and remains a separate future object rather than debt inside C-M5S-012. The C-M5S-011 de-boxed attribution likewise stays frontier; adding its positive mass shift can only increase compactness but is not used.

## Review and Promotion Plan

One independent claim review will freeze C-M5S-012's exact statement, implementation, evidence roles, dependencies and affected consumers. Exact identities, numeric trapped-surface/application evidence and continuation results receive distinct evidence roles without multiplying claims. A correction check will inspect only requested changes and directly affected edges. Promotion extracts reusable APIs/tests, updates the registry and release, renders generated docs, synchronizes accepted memory, validates the affected graph and records one content-addressed receipt. The author opens but does not self-merge the PR.

## Done Gate

P245 completes its user-authorized trapped-surface branch when the importable nonlinear construction, primary and independent oracles, affected-consumer replay, individual C-M5S-012 review, registry/release promotion, generated state, memory and validation receipt all agree. The positive object is the nonlinear fate theorem for the accepted certified source: a necessary Einstein constraint constructs and locates the trapped surface, while the homothetic theorem rules out pure dilation as an escape. No solver failure, unrefined sign, arbitrary-profile claim, localized fixed-J stress, black-hole interior solution or de-boxed quantum mass is smuggled into that statement.
