---
description: 'P229: machine-checked reproduction of the Preparata essential-instability QCD-vacuum program, extracting reusable chromomagnetic-vacuum primitives and conditional claims (Luca request, tracker issue #44)'
author: giuliano
created: '2026-08-11T15:45:00.000000+00:00'
updated: '2026-08-12T11:36:38+02:00'
tags:
- substrate-framework
- campaign-proposal
category: proposals
confidence: exploratory
status: active
---

## Question and Positive Deliverable
This campaign asks whether the Preparata program on the essential quantum
instability of the perturbative Yang-Mills vacuum reproduces under
independent machine-checked recomputation. The positive deliverable is a
set of importable framework primitives plus conditional claims covering: the
one-loop Savvidy effective potential in a covariantly constant
chromomagnetic background; the two-loop effective potential as a function of
the source variable Lambda/gB (exact convention fixed from the primary
sources); the fluctuation-stability verdict including the Nielsen-Olesen
unstable mode and the imaginary part; the color-by-color and transverse-mode
decomposition of the energy density; and an improved classical background
with requantization and effective-potential comparison against Savvidy. The
consumer-facing report (LaTeX + PDF, issue #50) cites only merged artifacts.
A no-go or an honest obstruction is attempt evidence, not completion.

## Base Release and Provenance
The accepted base is release v0.159.0 at framework commit 6d1f4e0, verified
against governance/releases/current.yaml on 2026-08-11. External sources:
Savvidy Phys. Lett. B 71 (1977) 133; Nielsen-Olesen Nucl. Phys. B 144
(1978) 376; Preparata Nuovo Cim. A 96 (1986) 366 (DOI 10.1007/BF02833896);
the two-loop source is being fixed by an independent extraction (agent
PreparataLit) before comparator values open. Sources actually read will be
recorded per attempt with URLs and equation numbers. Related accepted
surfaces read: C-VAC-001..004, C-NAG-001, C-NVP-001/002, C-GAU-001 (via
claims.yaml), and the unpromoted conditional modules twisted_casimir.py,
lattice_gauge.py, flux_tube_ensemble.py, lorentz_orbits.py.

## Source Inventory and Access Gate
Retrofitted 2026-08-11 after Dan's process review: this gate did not exist
when the campaign started; the campaign proceeded with the primary source
unchecked and while the extraction agent (PreparataLit) was still running
(AP-14/AP-15, skill quantitative-verification). The verified inventory:

| Source | Access status | What is extracted from it |
| --- | --- | --- |
| Preparata 1986, Nuovo Cim. A 96, 366 (DOI 10.1007/BF02833896) | IN HAND since 2026-08-11 17:26 UTC: `~/sources/preparata-qcd-vacuum/preparata1986-nuovo-cim-a96-366.pdf` (28 pp, text layer OK, full text extracted to `preparata1986.txt`); two byte-identical provenances — Dan's email (UID 68, Processed) and the DOI mirror | full paper available for equation-level verification (follow-up slice; the merged-pending PRs #51-#55 predate access and cover secondary sources) |
| Savvidy 1977, Phys. Lett. B 71, 133 | paywalled; open restatement in hand: arXiv:1910.00654 | one-loop potential, minimum, V_min formulas (equation level) |
| Nielsen-Olesen 1978, Nucl. Phys. B 144, 376 | paywalled; quoted via open review | unstable mode E^2 = p_z^2 - gB; Im V = (gB)^2/8 pi (independently re-derived, both routes agree) |
| Nielsen-Olesen 1978, Phys. Lett. B 79, 304 | paywalled; quoted via open review | electric vortex lines construction (motivation for C-CMV-004) |
| Ambjorn-Olesen 1980, Nucl. Phys. B 170, 60 and 265 | paywalled; quoted via open review eqs. (146)-(152) | spaghetti condensate ansatz, lattice minimization, eq. (152) minimum (convention mismatch documented) |
| Bordag-Skalozub 2022, EPJC 82, 390 | open, in hand: arXiv:2112.01043 | eq. (3) definitions, eqs. (57)-(58) T=0 result — the two-loop source |
| Bordag 2023, Symmetry 15, 1137 | open, in hand: MDPI full text | eqs. (146)-(152) flux-tube lattice construction |
| Cea 2023, arXiv:2311.14791 | open, in hand | SU(3) tachyon sectors (gH, gH/2, gH/2) |

Consequence recorded honestly: the objective as written ("the Preparata
program") is only partially verifiable until the Preparata PDF arrives; all
merged-pending artifacts are labeled conditional and secondary-source.

## Invariants, Conventions, and Allowed Imports
Frozen conventions, declared once for the whole campaign: Lorentzian
mostly-plus signature; background-field gauge with the gauge parameter
explicit; SU(N) with the background in a declared Cartan direction; the
covariant combination gB (never bare B) as the background variable; the
renormalization scale and the Lambda definition (scheme) written into every
formula before any coefficient comparison. Allowed imports: SymPy/SciPy
exact and numeric machinery under repository conventions; the accepted
claims listed above as convention cross-checks only; the conditional
unpromoted modules listed above as refactor inputs with their conditional
status declared; the external literature above as source objects under
audit, never as authority. Anything added later is explicit debt.

## Candidate Preregistration
Two independent computational routes are frozen before comparison.

| Candidate | Description | Assumptions | Parameters | Framework-fit prediction | Decisive test |
| --- | --- | --- | --- | --- | --- |
| A | Schwinger proper-time representation of the one-loop determinant in the covariantly constant background; zeta/proper-time renormalization; two-loop from published master integrals once sourced | background-field gauge, covariantly constant F | g, B, mu, gauge parameter xi | composes with C-VAC determinant conventions | V_1 coefficient and log structure match the source-independent Landau sum of B |
| B | Direct Landau-level mode sum E^2 = p_z^2 + gB(2n+1-2 s_3) with heat-kernel/smooth regularization, numeric-symbolic hybrid in the twisted_casimir style | same spectrum encoding; regulator declared | same plus regulator shape | reuses twisted_casimir machinery patterns | same target, evaluated independently; disagreement voids both until diagnosed |

## Selection Criteria and Blinding
Structural criteria, in order: convention compatibility with the frozen
declarations; exactness and independence from regulator artifacts; mutation
sensitivity; reuse value as importable primitives; maintenance cost.
Comparator blinding: the published two-loop coefficients are opened only
after both routes' one-loop outputs are frozen and agree.

## Proposed Claim Delta
Proposed provisional identifiers (registry, campaign, and memory collision
search run 2026-08-11; C-CMV prefix unclaimed): C-CMV-001 (fluctuation
spectrum and one-loop effective potential in the chromomagnetic background),
C-CMV-002 (two-loop effective potential in the source convention),
C-CMV-003 (unstable-mode spectrum, imaginary part, and the stability
verdict), C-CMV-004 (improved-background effective potential comparison).
Dependencies: convention claims C-GAU-001/C-VAC-001 as cross-checks only.
Consumers: the report of issue #50, future substrate gauge-sector work.
No supersedes before acceptance.

## Implementation and Oracle Plan
Reusable APIs land in src/substrate_framework/ (planned:
chromomagnetic_background.py for the spectrum and one-loop potential;
further modules per child issue), each with tests/test_*.py oracles and
mutations: flipped spin term, wrong degeneracy, wrong beta-function
coefficient must each break a check. Numeric oracles record algorithm,
precision, refinement, and stopping status per the oracle reference.
Campaign attempts live under campaigns/P229-preparata-qcd-vacuum-audit/
attempts/NNNN with manifests. Replay: PYTHONPATH=src .venv/bin/pytest
selectors per PR, scripts/validate.sh --full at any claim-promotion
boundary.

## Attempts and Continuation
Failed routes append with diagnosis and the next materially different
candidate; the effort contract memory/giuliano/efforts/
preparata-qcd-vacuum.md carries the attempt ledger.

## Debt Ledger
This ledger tracks new assumptions, residual unverified source claims
(including any paywalled-equation gaps), convention conversions, and broken
consumers introduced by this campaign. Empty before closure.

## Review and Promotion Plan
Each child PR (#45-#49 slices) gets independent review before merge; claim
promotion of C-CMV-001..004 happens only at claim-level adjudication with
independent rederivation, followed by release manifest update, rendered
docs, and memory synchronization. The report (issue #50) merges through the
docs pipeline and is delivered to Luca only after merge.

## Issue 59 Claim Narrowing
The 2026-08-12 post-merge audit separates reusable exact units from the wider
campaign narrative. Author-side review packets now cover only C-CMV-001 and
C-CMV-003. C-CMV-001 is narrowed to the conditional pure-SU(2) one-loop
spectrum, heat-kernel coefficient, declared-scheme potential/minimum,
imaginary-part magnitude, and an exact typed bridge to the accepted Pauli-half
and beta-function conventions. C-CMV-003 is narrowed to the exact one-loop
sector cancellation and canonical SU(3) T3 root magnitudes. Neither packet is
an acceptance decision.

C-CMV-002 remains proposed because the `mpmath` calculation shares the same
zeta construction rather than independently deriving the two-loop master
formula, and because the running-coupling two-loop expression is singular at
the one-loop minimum. C-CMV-004 remains proposed because the ansatz comparison
does not compute the inhomogeneous-background fluctuation determinant required
for requantization. The primary Preparata PDF is now in hand as recorded above;
the earlier access-gap consequence is retained as chronological attempt
evidence rather than current state.

## Done Gate
Closes only when every child issue is merged, the claims are adjudicated,
the report is delivered, and the debt ledger is empty.
