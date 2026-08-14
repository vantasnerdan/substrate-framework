---
description: 'Review paper Fulceri (2026-08-14) on breather-einbein effective metric as a campaign; map each claim to existing or new substrate-framework claims. After the overlap audit only two genuinely new candidate claims remain: C-OG-006 (flowing 1+1 effective metric) and C-IMP-001 (impedance-matching identity).'
author: vantasner
created: '2026-08-14T13:00:00Z'
updated: '2026-08-14T14:30:00Z'
tags:
- substrate-framework
- research-arc
- campaign-P226
- sine-gordon
- einbein
- effective-metric
- optical-geometry
- philosophical-out-of-scope
- overlap-audit
- harvest-completed-already
category: efforts
confidence: exploratory
status: active
---
# P226 — Fulceri Breather Einbein Effective-Metric Review

## Source
Fulceri, T. (2026-08-14), "From 1D Sine-Gordon Energy Lumps to Geodesic Propagation in an Effective Metric Induced by Substrate Inhomogeneity and Flow", v1.0, 14 August 2026, 13:51 CEST. AReTech — Advanced Research Technology Consulting. Self-archived v1.0 of the original sine-Gordon paper (v2.0, 10 May 2026, Zenodo 20110139) and the einbein tutorial (v1.0, 11 August 2026, Zenodo 21879560). PDF attached to issue #62.

## Goal and Success Contract
Classify each of the paper's 19 distinct claims into one of:
- **matched_existing**: paper claim is a verbatim or parameterized restatement of an accepted substrate-framework claim (C-SG-001..015, C-OG-001..002, C-MED-003).
- **new_candidate**: paper claim is a new exact algebraic statement with a sympy verifier; recorded as a candidate proposal for a future batch promotion.
- **textbook_citation**: paper claim is a standard textbook identity (einbein reparametrization, EOM-to-geodesic, Christoffel formula, Picard-Lindelof).
- **out_of_scope_interpretive**: paper claim is philosophical; not derivable; not eligible for accepted claim status.

Completion requires: (1) every claim classified; (2) every new candidate has a sympy verifier that closes within seconds; (3) every mutation test confirms the gate's sensitivity; (4) the interpretive-claim out-of-scope disposition is recorded explicitly; (5) no claim is accepted in this attempt; (6) the integrator identifies the outstanding review packets and merger.

## Authority and Prior Work
The accepted release is v0.159.0 at commit 640c47f, with 202 accepted claims. The source modules are C-SG-001..018 (sine-Gordon breather family), C-MED-003 (dimensional 1+1 sine-Gordon), and C-OG-001..005 (static 1+1 optical metric, TF constitutive map, dilaton, collective acceleration, complex-scalar action). Memory searches reference C-OG-001 as the closest existing analogue for the paper's flowing effective metric. The campaign evidence ledger includes P004 (optical dilaton), P012 (boosted breather kinematics), P044 (radial sine-Gordon oscillon), and P095 (dimensional sine-Gordon). The genuine unresolved objective is to identify exactly which paper claims are new and which are existing; for each new claim, build a sympy verifier.

The cross-check with already-merged harvest infrastructure identified the following overlap:
- src/substrate_framework/relativistic_particle.py and tests/test_einbein_3plus1d_tutorial.py, test_einbein_2plus1d_tutorial.py, test_relativistic_particle.py collectively cover 60 tests; this corpus was closed via PR #41, issue #40, P227.
- src/substrate_framework/pseudo_riemannian.py is the canonical Levi-Civita connection machinery.
- src/substrate_framework/constitutive.py asserts epsilon = rho*Theta/c^2 and mu^-1 = rho*Theta as a constitutive mirror.
- src/substrate_framework/boundary_scattering.py is the passive half-line scattering with z = zeta/c impedance.
- src/substrate_framework/gordon_metric.py is the 3+1 Gordon effective metric family (different parameterization from the paper's 1+1 lab-time form).

## Disposition Summary (per attempt 0001, after overlap audit)

Out of the paper's 19 distinct claims, 8 reproduce accepted substrate-framework claims exactly, 5 are textbook citations (3 of which are already covered by existing harvest infrastructure), 2 are new candidate claims with sympy verifiers, 1 is a consequence of a new candidate claim, and 2 are interpretive claims recorded as out-of-scope. The original 4 candidate claims reduced to 2 after the overlap audit. Each row below names one paper claim, the substrate-framework claim it matches (or "none" for new candidates), and the verified disposition.

| Paper claim | Substrate status | Disposition |
| --- | --- | --- |
| rest_breather_form | C-SG-001 | matched_existing |
| breather_energy | C-SG-002, C-MED-003 | matched_existing |
| boosted_breather | C-SG-008 | matched_existing |
| four_vector_dynamics | C-SG-008 | matched_existing |
| de_broglie_relation | C-SG-008 | matched_existing |
| two_frequencies | C-SG-015 (rest frame only) | matched_existing_partial |
| light_in_a_box_factorization | C-SG-013, C-SG-015 | matched_existing |
| square_root_lagrangian | already-harvested via PR #41 (relativistic_particle.py) | cites_existing_infrastructure |
| einbein_form | Polyakov 1987 textbook; already imported in relativistic_particle.py | textbook_citation_only |
| einbein_eom | textbook derivation; already covered in tests/test_einbein_3plus1d_tutorial.py::test_eq053_055_printed_christoffels_match_general_formula | textbook_citation_only |
| impedance_matching | none | candidate_for_new_claim_C-IMP-001 (S7 PASS); no-reflection propagation theorem is a separate deferred claim |
| effective_metric_static | C-OG-001 | matched_existing |
| effective_metric_flowing | C-OG-001 limit (V=0 only); the V != 0 part is new | candidate_for_new_claim_C-OG-006 (S3/S4/S5/S6/S8/S9/S10/S11/S12/S13 PASS) |
| metric_to_minkowski_limit | C-OG-001 limit (n=1, V=0) | matched_existing |
| christoffel_symbols | C-OG-006 consequence via pseudo_riemannian.metric_christoffel_from_derivatives | consequence_of_C-OG-006 (S12/S13 PASS) |
| residual_massless | C-OG-006 consequence (massless null constraint) | consequence_of_C-OG-006 (S9 PASS) |
| gravitational_equiv_principle | none | out_of_scope_interpretive |
| constructive_relativity | none | out_of_scope_interpretive |
| smooth_mass_zero_limit | Polyakov 1987 textbook property; already covered in tests/test_einbein_3plus1d_tutorial.py::test_eq057_060_massless_chain | textbook_citation_only |
| full_ivp_well_posedness | Picard-Lindelof textbook | textbook_citation_only |

## Verifier (attempt 0001)
The sympy gate script is at `campaigns/P226-fulceri-sg-breather-einbein-effective-metric/verify.py` (13 gates, 6 mutation tests). The substrate-metric module is at `src/substrate_framework/substrate_metric.py` with the pytest suite at `tests/test_substrate_metric.py` (13 mutation-sensitive gates). The 13 gates pass (G1 through G13) and the 6 mutation tests pass (M1, M2, M3, M3b, M4, M5). The mutations confirm gate sensitivity to load-bearing inputs. The targeted pytest scope (115 tests total) all pass.

## Overlap Audit Findings

The cross-check with already-merged infrastructure revealed that several pieces of the paper were already canonicalized through prior harvest PRs that closed their canonical issues without, however, promoting any claim to the accepted registry. This is consistent with the framework's separation of artifact merge from claim promotion. The substrate_metric module was refactored so its Christoffel-symbols path delegates to pseudo_riemannian.metric_christoffel_from_derivatives; the canonical connection machinery is the single source of truth.

The C-IMP-001 caveat is important: the no-reflection propagation theorem is not a consequence of the impedance identity alone. The candidate is restricted to the algebraic identity; the propagation theorem is a separate deferred claim requiring additional boundary/scattering analysis (and does not follow from rho*Theta = Z_0^2 by symbolic elimination alone).

## Interpretive-Claim Disposition
The paper's Sec. 6.8 makes two philosophical claims that the framework's discipline (claim = exact algebraic statement with mutation-sensitive verifier) does not accommodate:
- "The equivalence principle is recovered as a theorem of continuum mechanics."
- "Lorentz invariance is dynamical, not fundamental."

These are interpretive positions, not derivations. They are not eligible for accepted claim status. They are noted in the paper's reference list as belonging to the Brown 2005 / Bell 1976 / Selleri 1996 tradition and as philosophically aligned with Schmelzer's GLET. The framework takes no position on these interpretive claims; it acknowledges the paper's self-identification as a Lorentzian ether theory and does not contest that self-identification.

## Outstanding Work
The remaining work for P226 is: an independent review packet per new candidate claim (2 claims: C-IMP-001, C-OG-006); a decision on batch promotion under P226 vs separate campaigns per claim; explicit acknowledgment that the no-reflection propagation theorem is a separate deferred claim not a consequence of C-IMP-001 alone; the final acceptance gate is a distinct merger (non-self-merge rule); and the boundary that no claim may be promoted that identifies a Standard-Model particle, recovers GR as a fundamental theory, or asserts Lorentz invariance is derived.

## Open Questions
The following questions remain open for the distinct reviewer to decide. First, what is the precise disposition of the impedance-matching identity vs the impedance-via-c^2 relation already in constitutive.py? The two are distinct statements: constitutive.py asserts epsilon = rho*Theta/c^2 and mu^-1 = rho*Theta as a constitutive mirror; C-IMP-001 asserts the invariance under the constitutive map. Both should coexist, with explicit cross-references. Second, should the Christoffel-symbols evidence (verifier S12/S13) be cited as part of the C-OG-006 evidence rather than as a separate gate? Likely yes, but this is a reviewer decision.
The previously open questions have been resolved: C-OG-008 (Christoffel symbols) is a consequence of C-OG-006 via the standard Levi-Civita connection and drops to consequences of C-OG-006; the smooth M->0 limit (C-OG-009) is a textbook property of the einbein reparametrization and is cited as textbook; and C-IMP-001 stands as the algebraic identity only, with the propagation theorem as a separate deferred claim requiring additional boundary/scattering analysis.
