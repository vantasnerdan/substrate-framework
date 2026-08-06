---
description: Construct general relativity and its coupling from sine-Gordon microscopic dynamics
author: axis
created: '2026-08-06T18:22:18Z'
updated: '2026-08-06T18:55:44Z'
tags:
- substrate-framework
- campaign-proposal
category: proposals
confidence: exploratory
status: active
---

## Question and Positive Deliverable
This campaign asks whether accepted sine-Gordon microscopic dynamics can support one coherent framework-native construction of low-energy general relativity.

The positive deliverable is an end-to-end dependency chain that starts from accepted sine-Gordon primitives, constructs the required collective 3+1 structure, derives the Einstein dynamics and gravitational coupling from upstream parameters, and produces one static weak-field prediction plus one radiative prediction through reusable package APIs and tests. A no-go, conditional substitution into known GR formulas, or reduced-scope report does not complete the campaign.

## Base Release and Provenance
The accepted authority boundary is release `v0.159.0`, whose manifest pins source baseline `substrate@6d1f4e0`; the clean framework checkout began at `substrate-framework@e38ca619` on `main`.

The campaign question is issue 1. Closed PR 2 is attempt evidence only: its maintainer comments identify cutoff conflation, an unclosed additive inverse-coupling baseline, a missing 1+1-to-3+1 quantum lift, imported GR dynamics, and unconstrained quantum-scale API inputs as anti-patterns. No result or comparator from that branch is reused as authority.

The frozen proposal manifest passed repository validation before source exploration. Accepted claims read at registry and source level are `C-MED-003`, `C-SG-011`, `C-SG-012`, `C-SG-017`, `C-GRV-001`, and `C-GW-001` through `C-GW-004`; neighboring ceiling claims `C-GOR-001`, `C-GOR-002`, `C-STG-001`, `C-STG-002`, `C-MOM-001`, and `C-PDE-001` were also audited. Their canonical owners and consumers were traced with GitNexus before implementation. No predecessor migration unit is assigned by the issue.

## Invariants, Conventions, and Allowed Imports
The construction preserves the accepted sine-Gordon equations, scale meanings, authority order, claim statuses, and current accepted registry without reinterpretation.

Permitted framework inputs are accepted claims in `v0.159.0` whose exact statements and assumptions are verified at source after the freeze gate. Permitted external imports are standard differential geometry; the local-Rindler Clausius-to-Einstein theorem and Unruh temperature for candidate A; teleparallel differential identities for candidate B; and the massless-spin-2 consistent-self-coupling theorem for candidate C. Every import must expose its hypotheses and coefficient interface. Observed gravitational constants, Planck scales, solar-system values, and textbook benchmark outputs are comparators, not derivation inputs.

The final selected route must use one constrained action scale rather than independent `action_quantum`, `action_scale`, and `hbar` inputs. It must derive or explicitly construct the 3+1 degrees of freedom, measure, source conservation, dynamical metric equation, and physical radiation flux it uses.

## Candidate Preregistration
The candidate set spans thermodynamic, geometric, and collective-spin-2 mechanisms before any observational comparator is opened.

| Candidate | Description | Assumptions | Parameters | Framework-fit prediction | Decisive test |
| --- | --- | --- | --- | --- | --- |
| A | Quantized sine-Gordon link cells supply a counted boundary entropy; local causal-horizon entanglement equilibrium then yields Einstein dynamics. | A derived transverse phase-cell lattice, a constrained action quantum, local equilibrium, and the explicit hypotheses of the Clausius-to-Einstein theorem. | Accepted medium scales plus discrete cell degeneracy fixed by the sine-Gordon topological sector. | Strong if the cell Hilbert space, area density, Lorentzian lift, and conserved heat flux follow without fitted constants. | Derive entropy per area and the Einstein coefficient from the same microscopic object, with no bare additive coupling and no GR equation inserted as an action premise. |
| B | Four coupled sine-Gordon phase-frame sectors form a non-integrable coframe whose continuum link action reduces uniquely to the teleparallel equivalent of GR. | A phase-frame lattice, local Lorentz redundancy, and a microscopic rule selecting the torsion invariant combination. | Accepted stiffness and length/action scales, with no free TEGR coefficient ratios. | Strong only if the TEGR combination is forced rather than tuned. | Derive the torsion-scalar coefficients and coupling from the lattice expansion; reject if the Einstein-equivalent combination is an unexplained input. |
| C | A 3+1 sine-Gordon defect condensate has a protected massless helicity-2 collective mode whose consistent nonlinear self-coupling closes to GR. | A transverse defect network, exact symmetry protecting helicity two, and the hypotheses of the spin-2 bootstrap theorem. | Accepted defect energy and spacing; any condensate order parameter must be fixed upstream. | Strong only if unwanted scalar/vector modes are gapped and the coupling normalization is computable. | Exhibit the mode measure, two helicities, universal conserved-source coupling, nonlinear closure, and a derived coupling scale. |

## Selection Criteria and Blinding
Structural compatibility is evaluated before any empirical closeness.

The frozen ordering is: accepted-invariant compatibility; explicit 1+1-to-3+1 dependency closure; generation rather than import of metric dynamics; absence of additive or fitted coupling freedom; assumption and parameter economy; correct symmetries, dimensions, conservation laws, topology, and low-energy limits; shared upstream source for static and radiative arms; numerical robustness; and reusable implementation cost. Observed values and benchmark outputs remain blinded until a candidate comparison record freezes the mechanism, equations, conventions, and structural tests.

## Proposed Claim Delta
The campaign reserves a four-rung provisional claim ladder after a repository-wide collision search found no existing `C-EGR-*` identifier.

- `C-EGR-001`: the selected sine-Gordon microscopic construction produces a finite 3+1 collective geometry and a derived gravitational coupling with no independent additive baseline.
- `C-EGR-002`: the construction's low-energy collective equation is the Einstein equation, with explicit source conservation and two radiative tensor degrees of freedom.
- `C-EGR-003`: the same upstream model yields a static weak-field observable through the derived coupling.
- `C-EGR-004`: the same upstream model yields a radiative observable and physical flux through the derived coupling.

Dependencies, exact statements, evidence paths, and consumers will be narrowed after accepted-source exploration. No accepted claim is superseded at proposal scope; the new claims challenge only any explicit open ceilings they close.

## Implementation and Oracle Plan
The selected route will place canonical equations, constrained parameter objects, source maps, static solvers, and radiative solvers under `src/substrate_framework/`; proposal code will be thin orchestration.

SymPy or direct symbolic algebra will check exact coefficient matching, dimensions, conservation identities, weak-field limits, gauge content, and mutation-sensitive residuals. If candidate A survives, the theorem interface will be encoded so the entropy density and heat-flux normalization are independent inputs derived by separate microscopic functions before coefficient elimination. If candidate B survives, differential-form expansion will separately audit non-integrability, torsion invariants, boundary-term equivalence, local Lorentz behavior, and coefficient uniqueness. If candidate C survives, the quadratic mode spectrum and helicity projectors will be derived before applying the nonlinear bootstrap theorem.

SciPy will be used only for genuinely unresolved source dynamics or flux integrals, with stated equations, domain, data, discretization, tolerances, status checks, scale-relative error norms, refinement, conservation checks, and an independent route. Canonical sampled integration will use `trapezoid_integral`; executable syntax will be checked for direct, imported, dynamic, and eager-fallback `np.trapz` access. Static and radiative predictions will be derived from the same selected field equation and same conserved source, not pasted benchmark formulas. Load-bearing coefficients, scale relations, source conservation, gauge conditions, and selected microscopic degeneracies will each receive wrong-sign, wrong-normalization, or inconsistent-scale mutations that must fail.

Repository graph analysis and source-level import tracing will enumerate direct and indirect consumers. The final gate is targeted claim verification, affected package tests, `scripts/validate.sh`, and `git diff --check`, with validation and commit run as separate processes.

## Selected Construction and Result
Candidate A is selected. Candidates B and C are preserved as append-only scientific rejections: the accepted link dynamics do not select the TEGR torsion coefficients, and they do not supply a protected massless helicity-two mode or an evasion of the conserved-stress obstruction.

The selected collective weave uses four shared length-`ell` sine-Gordon links per tetrahedral cell. Their directions close and have covariance `I/3`; uniformly distributed coarse orientations give plane-crossing density `1/ell^2`. The exact parity-related `Q=+1` and `Q=-1` sectors form the declared equilibrium topological doublet, fixing entropy `log(2)` per crossing. A collective clock plus the four deformed link covectors reconstruct the Lorentzian coframe metric and volume measure.

Local horizon equilibrium then yields the nonlinear Einstein equation rather than inserting an Einstein action. The same upstream action `J=sqrt(lambda*T)`, speed `c=sqrt(T/lambda)`, and entropy density `log(2)/ell^2` give

`8*pi*G/c^4 = 2*pi/(mu*log(2))` and `G=T^2/(4*lambda^2*mu*log(2))`.

There is no cutoff, independently variable quantum action, supplied Newton constant, fitted dimensionless coefficient, or additive inverse-coupling baseline. The flat equilibrium branch fixes the cosmological integration constant to zero.

The dimensional sine-Gordon Noether stress is embedded with transverse delta density as an isolated, exactly conserved on-shell 3+1 source. One exact rest breather supplies both predictions. At inverse width `eta=0.03`, its profile compactness is `0.010387404294400536`. The full-retardation TT calculation, required because `omega/eta>30`, converges to `P/(mu*c)=0.002606646892787329`; the independent leading analytic transform gives `0.002608524438585714`, and the one-cycle loss fraction is `0.03413629282664175`.

## Attempts and Continuation
The append-only attempt ledger records each candidate's concrete derivation, command environment, evidence, verdict, failure mechanism, and next route.

Closed PR 2 is external attempt evidence and is not revalidated. New failed routes remain proposal evidence and trigger repair, reformulation, or the next materially distinct candidate; they never become the requested endpoint.

## Debt Ledger
The debt ledger tracks every new premise, import, parameter, residual, convention interface, and affected consumer introduced by the selected construction.

The implementation closes candidate selection, accepted-claim dependency tracing, the 3+1 collective metric and measure, the single-action constraint, the total coupling, isolated source conservation, and both oracle arms. The remaining review surface is epistemic rather than hidden implementation debt: maintainers must adjudicate the explicitly declared tetrahedral-link, Haar-isotropy, quantum-doublet, unique-action, and local-equilibrium hypotheses. No empirical constant selected or normalized the construction.

## Review and Promotion Plan
The pull request will expose each provisional claim, raw derivation artifact, verifier, mutation evidence, API, tests, dependency impact record, and candidate comparison for maintainer review.

Claims remain proposed and are reviewed individually. Axis will not self-promote them into the accepted registry or release before maintainer adjudication. If accepted, promotion must extract canonical APIs, update the claim graph and pinned release, render generated documentation, synchronize accepted memory, and preserve rejected candidates as immutable attempt evidence. No predecessor migration disposition is part of this issue unless accepted-source exploration establishes one.

## Done Gate
The campaign closes only when one selected construction satisfies the complete positive contract and leaves no unresolved dependency, coupling, source, dynamics, validation, consumer, or narrative debt.

If a gate fails, the next concrete attempt remains active. The current delivery target is a reviewable pull request, not an unearned accepted-claim or release signal.
