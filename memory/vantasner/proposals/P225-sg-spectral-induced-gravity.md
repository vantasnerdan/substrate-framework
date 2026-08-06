---
description: Construct the sine-Gordon-to-general-relativity chain by deriving the gravitational coupling from the accepted sG mass tower via curved-background one-loop induction
author: argus
created: '2026-08-06T19:40:00Z'
updated: '2026-08-06T19:40:00Z'
tags:
- substrate-framework
- campaign-proposal
- spectral-induced-gravity
category: proposals
confidence: exploratory
status: active
---
# P225 sine-Gordon Spectral Induced Gravity

## Question and Positive Deliverable

Campaign issue vantasnerdan/substrate-framework#1 asks: can the framework
construct a coherent path from sine-Gordon dynamics to general relativity?
The positive object is a single end-to-end, framework-native construction that
starts from the accepted sine-Gordon microscopic dynamics (C-MED-003 physical
coefficients; C-SG-001..019 breather sector; C-SG-007 action lattice), develops
the effective gravitational structure (C-GOR-001 Gordon metric; C-STG-001/002
Einstein-scalar sector), derives the gravitational coupling and scale from
upstream model parameters (closing the dimensionless coordinate C-GRV-001
proved unreachable by dimensional analysis and P074/AS3 left free), reaches a
regime that reproduces GR's defining low-energy dynamics (Newtonian limit and
linearized TT radiation), and demonstrates one static weak-field prediction and
one radiative prediction from the same upstream model. A no-go, residual,
ledger-only audit, or reduced-scope report does not complete this campaign.

## Base Release and Provenance

Accepted base release `v0.159.0`, source baseline `substrate@6d1f4e0`, working
commit `e38ca619` (Close framework migration with MR6 audit). Modules read at
source: `sine_gordon.py`, `dimensional_sine_gordon.py`, `gordon_metric.py`,
`gordon_scalar_compatibility.py`, `optical_geometry.py`, `einstein_scalar.py`,
`spherical_einstein_scalar.py`, `induced_gravity.py`,
`gravity_scale_confrontation.py`, `scale_provenance.py`, `effective_actions.py`,
`action_scales.py`. Campaigns read: P074 (induced gravity scaling ledger),
P004/P010 (optical dilaton), P036-P039 (GW quadrupole machinery),
P141/P142/P178/P179 (dilaton, Gordon, Einstein-scalar audits) by registry
summary. Prior gravity campaigns established conditional ledgers only; every
existing module docstring explicitly disclaims deriving a Newton constant,
Einstein dynamics from upstream physics, or the dimensionless coefficient.

## Invariants, Conventions, and Allowed Imports

The campaign preserves: mostly-plus signature (C-GOR-001, C-STG-001); the
C-GRV-001 dimensional ceiling (G = q*a^2*c^3/hbar with q dimensionless and
free under dimensional analysis alone — the campaign must fix q by supplying
physics, not by re-running dimensional analysis); C-SG-007's lattice premise
status (h is a declared action increment, not derived); C-MED-003's scale
identifications; the C-GW-001..010 TT/quadrupole conventions (A=2G, B=1/(32
pi G) normalized; triple-convention equivalences); the optical sector's
explicit ceiling that its 1+1 source equation does not normalize kappa to G.
Allowed imports: all accepted claims of v0.159.0; the heat-kernel one-loop
effective-action structure for a real scalar (standard result, role explicit:
the Seeley-DeWitt a2 coefficient supplies the induced R term; its coefficient
(1/6-xi)/2 with the (4*pi)^-2 measure is imported as the regulator-relative
induction formula); the Schwarzschild exterior solution and its weak-field
predictions (standard imports with explicit role); comparator values (the
observed Newton constant, Planck scale) remain blinded and unused.

## Candidate Preregistration

Three plausible concepts are registered before any coefficient is evaluated:
spectral Sakharov induction on the accepted mass tower, the optical-dilaton
direct route, and elastic collective-coordinate gravity.

| Candidate | Description | Assumptions | Parameters | Framework-fit prediction | Decisive test |
| --- | --- | --- | --- | --- | --- |
| A | Spectral Sakharov induction: the one-loop vacuum polarization of the accepted sG mass tower (phonon + C-SG-007 breather levels + kink) on the emergent Gordon geometry induces the Einstein-Hilbert term; the tower edge E=8pi is the spectrum's own cutoff, making log(Lambda^2/m^2) positive and finite per level; q becomes a pure computable number of the model. | heat-kernel induction formula; J_scale=hbar quantum premise; Gordon metric as the emergent geometry | h (declared), xi=0 (canonical), no new dimensionless inputs | composes with C-GOR-001, C-STG-002, C-GW-004 without new conventions | q is a closed-form tower sum; static Newtonian limit and TT quadrupole power both follow with the same derived G |
| B | Optical-dilaton direct: the 1+1 optical metric (C-OG-001..003) is already an exact Einstein-like system; coupling from the constitutive map and C-MED-002 medium dictionary. | optical metric is the spacetime metric | none new | structurally 1+1: no TT projector, no quadrupole | cannot produce a radiative prediction; fails the deliverable's radiative arm by dimension, not by tuning |
| C | Elastic collective-coordinate gravity: metric from lattice strain; coupling from elastic moduli (T, lambda). | strain metric identification | new collective-coordinate map | vector/phonon gravity; no TT-2 tensor structure from scalar elasticity | fails the spin-2 TT requirement (C-GW-002 projector) structurally |

## Selection Criteria and Blinding

Frozen before any coefficient is evaluated:
1. consistency with accepted invariants and conventions (mostly-plus; C-GW conventions; C-GRV-001 ceiling);
2. explanatory reach: must yield both a static weak-field and a TT radiative prediction from one upstream model;
3. assumption and parameter economy: no undeclared dimensionless inputs; every premise named;
4. correct symmetries, dimensions, topology, and limits (q dimensionless; correct C-GRV-001 monomial; sensible h->0 and tower-edge limits);
5. compatibility with accepted sectors (Gordon, Einstein-scalar, GW machinery reuse, not redefinition);
6. numerical robustness and implementability (exact symbolic ledger plus double-precision evaluation with refinement).
Comparator blinding: no empirical G, Planck scale, or astrophysical benchmark is inspected or used until the construction is frozen; confrontation is out of scope for this campaign.

## Proposed Claim Delta

Provisional identifiers (collision-searched against governance/claims.yaml,
campaigns/, and durable memory; no C-SIG identifiers exist anywhere in the
repository history):
- C-SIG-001: the sG mass-tower induction coefficient — the exact symbolic
  ledger giving 1/G_ind as the tower sum with regulator, xi, and per-species
  masses explicit, in the C-GRV-001 normalization (q extracted as a pure
  number).
- C-SIG-002: the end-to-end static weak-field claim — vacuum exterior of the
  Einstein-scalar sector with derived kappa is Schwarzschild; Newtonian
  potential and light-deflection angle with G from C-SIG-001.
- C-SIG-003: the end-to-end radiative claim — breather longitudinal-moment
  quadrupole power (C-GW-004/C-SG-009/C-SG-010 machinery) with the declared
  A=2G input replaced by the C-SIG-001 derived G; cycle-averaged power in
  upstream model units.
Status: proposed (review: unaudited; verification: symbolic/formal plus
numeric evidence; compatibility: compatible extension; epistemic: proposed).
No existing claim is challenged; C-GRV-001 is composed, not superseded.

## Implementation and Oracle Plan

New importable modules under `src/substrate_framework/`:
- `sg_spectral_induction.py`: exact sympy ledger for the per-species induced
  inverse coupling, the tower sum (phonon + breather levels from
  `sine_gordon.breather_action_lattice_energy`, kink), the C-GRV-001 q
  extraction, and h-dependence; numeric evaluation with mpmath-grade exact
  evaluation of the symbolic sum.
- `sg_gravity_chain.py`: the composition API — upstream (lambda, T, mu, h) to
  scales (C-MED-003) to derived G/kappa to (i) static exterior metric and
  deflection prediction, (ii) breather quadrupole power with derived G.
Oracles: exact symbolic identities (dimensional closure vs C-GRV-001
monomial; q independence of unit standard); mutation sensitivity (drop the
kink, change xi, shift h, remove the top level — q must move by the predicted
amount; wrong-sign xi flips sign); limiting cases (h->0 continuum divergence
diagnosed, not hidden); numeric cross-check of the symbolic tower sum by
direct float summation; static arm checked against the exact Schwarzschild
coefficients; radiative arm checked against the accepted C-SG-010
quadrature bracket (379.4646380687 < <mu'''^2> < 379.4646380688).
Tests under `tests/test_sg_spectral_induction.py` and
`tests/test_sg_gravity_chain.py`. Campaign verifier
`proposals/P225-sg-spectral-induced-gravity/verify.py` runs with
`PYTHONPATH=src` and imports only package APIs. Global replay:
`scripts/validate.sh`. NumPy compatibility preflight: canonical integration
via `trapezoid_integral`; mutable scripts use `np.trapezoid`; no legacy
`np.trapz` or eager nested getattr fallback anywhere in campaign code.

## Attempts and Continuation

Append-only attempts under `proposals/P225-sg-spectral-induced-gravity/attempts/`
with candidate, command, stdout/stderr, verdict, and diagnosis. Failure of A
(self-cancellation, wrong limits) routes to a declared multi-sector field
content or a foundational-revision proposal; B and C failures are recorded as
structural-rejection evidence, not campaign failure.

## Debt Ledger

Tracks: the quantum premise (J_scale = hbar) and its provenance scope; the
regulator choice (tower-edge cutoff) versus continuum-cutoff alternatives;
xi=0 as the canonical-scalar premise; the breather-as-loop-species
semiclassical premise; any convention conversions (kappa = 8*pi*G/c^4) named
explicitly. Each must be discharged or declared in-claim before promotion
review.

## Review and Promotion Plan

Claim-level independent review by a fresh reviewer agent receiving raw
artifacts and criteria only (delegation per the repository skill). Claims are
proposed in `governance/claims.yaml` only to the proposal's proper scope; the
PR presents the construction for maintainer adjudication. Acceptance, release
promotion, generated documentation, and accepted-memory synchronization are
the maintainer's gates after review; this campaign's PR records claims as
proposed with full verification evidence and leaves the registry's accepted
set untouched. Validation and commit run in separate invocations.

## Done Gate

The campaign closes only when the end-to-end construction exists in importable
framework code, q is a computable pure number of the upstream sG model, both
predictions evaluate from the same upstream parameters, all proposed-claim
evidence is durable, and `scripts/validate.sh` passes in full.
