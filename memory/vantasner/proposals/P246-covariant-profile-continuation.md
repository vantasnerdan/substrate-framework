---
description: Construct a self-consistent covariantly relaxed nonlinear Einstein-M5 clock
author: Main
created: '2026-08-26T18:25:56+02:00'
updated: '2026-08-26T18:25:56+02:00'
tags:
- substrate-framework
- campaign-proposal
- nonlinear-gravity
- m5
category: proposals
confidence: exploratory
status: active
---

## Question and Positive Deliverable
P246 asks whether the accepted confined-clock sector at purely induced `B=0`, `xi=0` forms a self-consistent nonlinear spherical gravitating object when its M5 profiles and stress vary covariantly with the metric. The positive deliverable is a typed reduced Einstein-M5 action, regular boundary data, and a co-solved branch at the accepted coupling: either a regular horizonless clock or a genuine horizon-bearing M5 configuration. A further frozen-source no-go or solver termination is attempt evidence, not completion.

## Base Release and Provenance
The accepted base is v0.166.0 at `substrate-framework@3e3a6133030b3041515d07cbeef74eeb3624b99b`. P246 is developed on a branch stacked above draft PR #175 for provenance, but draft C-M5S-012 is only proposal evidence until merged. Accepted inputs are C-M5S-003/004/007/010/011, C-GRV-002, C-STG-001/002, and canonical M5, spherical-Einstein, and numerics modules at the pinned base.

## Source Inventory and Access Gate
All load-bearing sources are in hand in the repository. `governance/claims.yaml` supplies the accepted coupling and clock ledgers; `src/substrate_framework/m5_covariant_action.py` and `m5_stationary_fields.py` supply the action and stationary ansatz primitives; C-STG-001/002 and the spherical Einstein modules supply gravitational conventions. P245 supplies branch-local evidence that the frozen source and pure dilation are trapped, while explicitly excluding arbitrary profile relaxation.

## Invariants, Conventions, and Allowed Imports
The campaign preserves mostly-plus signature, the accepted Einstein-Hilbert normalization, areal radius outside horizons, the purely induced `B=0`, `xi=0` endpoint, the M5 coefficients and topology, and the fixed-J clock sector. Stress must be derived by metric variation, the fixed-J reduction must be varied before any localization, and no scalar field map or equation of state may be assumed. Regular solutions must satisfy the matter equations, Einstein constraints, and on-shell stress conservation.

## Candidate Preregistration
Five mechanisms are frozen before new branch outcomes are inspected. Candidate A derives and continues the full covariant spherical M5 profile and metric. Candidate B enlarges the profile space so kinetic, potential, radial, and tangential stresses can redistribute at fixed charge. Candidate C uses coupling and `xi` continuation to locate folds but must return to the target endpoint. Candidate D constructs a baseline-dominated seed and tests connectivity back to `B=0`. Candidate E constructs a horizon-regular M5 object if the regular-origin family terminates at a horizon.

## Selection Criteria and Blinding
Selection prioritizes exact action compatibility, preservation of the target endpoint and topology, absence of fitted constitutive data, finite invariant stress, full residual closure, correct weak-coupling and vacuum limits, parameter economy, and independent-route agreement. Existing accepted and P245 compactness values are unavoidably known. New branch endpoints, folds, masses, central compactness values, and horizon radii stay blinded until the manifest and this memory validate.

## Proposed Claim Delta
P246 reserves C-M5S-013 for the complete fixed-J stress and symmetry classification and C-M5S-014 for the positive horizon-penetrating initial-data construction. C-M5S-013 depends on the accepted clock source and spectrum scopes; C-M5S-014 depends on the classical source convention plus the exact nonlinear constraint. Neither claim says that the stationary axisymmetric one-clock field equations have been solved.

## Implementation and Oracle Plan
The implementation will expose pure APIs for the stationary axisymmetric action, complete metric-derived stress including momentum, fixed-J Routhian, regular-axis/origin and regular-horizon data, residual assembly, observables, and continuation. Exact SymPy checks will cover signs, metric variation, conservation identities, zero-inertia limits, and weak-coupling limits. Numerical work will use nondimensional variables, sparse Jacobians, Newton correction with pseudo-arclength near folds, explicit residual norms, mesh/domain/tolerance continuation, solver-status gates, load-bearing sign and coupling mutations, and an independent route.

## Attempts and Continuation
Attempts are append-only under `proposals/P246-covariant-profile-continuation/attempts/`. Technical failures change method; profile-space obstructions enlarge the ansatz; folds trigger pseudo-arclength; a regular-origin horizon triggers the horizon-regular candidate. Reduced coupling and baseline branches are continuation tools rather than substitutions for the target. The campaign continues until it constructs a positive nonlinear object in the declared scope.
Attempt 0003 established the required representation change before a wrong spherical solver was built. The split-free uniaxial field commutes exactly with its clock generator and has zero inertia, while the accepted split-active root has 2.3785% angular energy variation and 7.1874% integrated tangential-pressure anisotropy. Candidate A therefore continues with axisymmetric gravity and angular field modes; P245's spherical average remains only a compactness diagnostic.
Attempts 0004-0009 resolved all preregistered structural comparisons. The full stress has a 0.5095% azimuthal momentum fraction and 1.6273% radial-polar shear, selecting stationary axisymmetry. Four regular even split modes reduce energy by only 3.12e-7 relative, below the frozen nontriviality gate. Counterrotation cancels momentum but preserves energy and compactness exactly. Horizon avoidance requires either xi about -148.802 or a baseline 892.812 times the induced inverse coupling. The positive surviving construction is the Painleve-Gullstrand initial-data geometry for the accepted classical source: inner marginal radius 0.0935458851 and matched-vacuum outer marginal radius 5121.5304868, with exact Hamiltonian/momentum constraint closure and a nonsingular metric determinant at both surfaces.

## Debt Ledger
The promoted C-M5S-013/C-M5S-014 boundary has no hidden debt: the action/stress map, source convention, exact constraints, numeric refinement, exclusions, tests and impact inventory are explicit. The stationary axisymmetric fixed-J Einstein-M5 matter-plus-metric solve remains active campaign frontier rather than a promise inside either claim. P245's frozen-source result and unrelated radiative attribution remain separately scoped.

## Review and Promotion Plan
The frozen promotion transaction includes C-M5S-013 and C-M5S-014, the new stress/relaxation and horizon-penetrating APIs, exact and numerical evidence, affected M5 and gravity consumers, and one independent claim review per claim in one bounded pass plus one correction check if needed. Promotion requires the registry, release, generated documentation, durable memory, and one impact-selected validation receipt to agree.

## Done Gate
P246 does not close merely because its two independently useful claims promote. The positive horizon-penetrating classical-source object satisfies the black-hole-like candidate, but the broader one-clock objective remains open until a stationary axisymmetric fixed-J Einstein-M5 branch reaches the accepted purely induced coupling with full matter and metric residuals. A no-go, fold, disconnected baseline branch, or the scoped initial-data object cannot be relabeled as that stronger result.
