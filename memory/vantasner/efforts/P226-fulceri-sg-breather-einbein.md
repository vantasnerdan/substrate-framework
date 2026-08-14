---
description: Review paper Fulceri (2026-08-14) "From 1D Sine-Gordon Energy Lumps to Geodesic Propagation in an Effective Metric Induced by Substrate Inhomogeneity and Flow" as a campaign; map each claim to existing or new substrate-framework claims; reproduce the new candidate claims by sympy verifier.
author: vantasner
created: '2026-08-14T13:00:00Z'
updated: '2026-08-14T13:00:00Z'
tags:
- substrate-framework
- research-arc
- campaign-P226
- sine-gordon
- einbein
- effective-metric
- optical-geometry
- philosophical-out-of-scope
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
- Accepted release: v0.159.0 at commit 640c47f.
- Accepted claim count: 202.
- Source modules: C-SG-001..018 (sine-Gordon breather family), C-MED-003 (dimensional 1+1 sine-Gordon), C-OG-001..005 (static 1+1 optical metric, TF constitutive map, dilaton, collective acceleration, complex-scalar action).
- Memory searches: prior effort entries (this session) reference C-OG-001 as the closest existing analogue for the paper's flowing effective metric.
- Campaign evidence: P004 (optical dilaton), P012 (boosted breather kinematics), P044 (radial sine-Gordon oscillon), P095 (dimensional sine-Gordon).
- Genuine unresolved objective: identify exactly which paper claims are new and which are existing; for each new claim, build a sympy verifier.

## Disposition Summary (per attempt 0001)

Out of the paper's 19 distinct claims, 8 reproduce accepted substrate-framework claims exactly, 5 are textbook citations, 4 are new candidate claims with sympy verifiers, and 2 are interpretive claims recorded as out-of-scope. Each row below names one paper claim, the substrate-framework claim it matches (or "none" for new candidates), and the verified disposition.

| Paper claim | Substrate status | Disposition |
| --- | --- | --- |
| rest_breather_form | C-SG-001 | matched_existing |
| breather_energy | C-SG-002, C-MED-003 | matched_existing |
| boosted_breather | C-SG-008 | matched_existing |
| four_vector_dynamics | C-SG-008 | matched_existing |
| de_broglie_relation | C-SG-008 | matched_existing |
| two_frequencies | C-SG-015 (rest frame only) | matched_existing_partial |
| light_in_a_box_factorization | C-SG-013, C-SG-015 | matched_existing |
| square_root_lagrangian | none | candidate_for_new_claim_pending_review |
| einbein_form | none | textbook_citation_only (Polyakov 1987) |
| einbein_eom | none | textbook_citation_only |
| impedance_matching | none | candidate_for_new_claim_C-IMP-001 (G9 PASS) |
| effective_metric_static | C-OG-001 | matched_existing |
| effective_metric_flowing | C-OG-001 limit (V=0) | candidate_for_new_claim_C-OG-006 (G5/G6/G10/G11 PASS) |
| metric_to_minkowski_limit | C-OG-001 limit | matched_existing |
| christoffel_symbols | C-OG-002 (TF, V=0 only) | candidate_for_new_claim_C-OG-008 (G7 PASS) |
| residual_massless | none | candidate_for_new_claim_C-OG-007 (G8 PASS) |
| gravitational_equiv_principle | none | out_of_scope_interpretive |
| constructive_relativity | none | out_of_scope_interpretive |
| smooth_mass_zero_limit | none | candidate_for_new_claim_C-OG-009 (G2+G8 cross) |
| full_ivp_well_posedness | none | textbook_citation_only (Picard-Lindelof) |

## Verifier (attempt 0001)
Sympy gate script at `campaigns/P226-fulceri-sg-breather-einbein-effective-metric/verify.py`.
- 11 gates pass (G1 through G11).
- 5 mutation tests pass (M1 through M3b, M4).
- Mutations confirm gate sensitivity to load-bearing inputs.

## Interpretive-Claim Disposition (explicit)
The paper's Sec. 6.8 makes two philosophical claims that the framework's
discipline (claim = exact algebraic statement with mutation-sensitive
verifier) does not accommodate:
- "The equivalence principle is recovered as a theorem of continuum mechanics."
- "Lorentz invariance is dynamical, not fundamental."

These are interpretive positions, not derivations. They are not eligible for
accepted claim status. They are noted in the paper's reference list as
belonging to the Brown 2005 / Bell 1976 / Selleri 1996 tradition and as
philosophically aligned with Schmelzer's GLET. The framework takes no
position on these interpretive claims; it acknowledges the paper's
self-identification as a Lorentzian ether theory and does not contest that
self-identification.

## Outstanding Work
- Independent review packet per new candidate claim (5 claims: C-IMP-001, C-OG-006, C-OG-007, C-OG-008, C-OG-009).
- Decision: batch promotion under P226 vs separate campaigns per claim.
- Final acceptance gate: distinct merger required (non-self-merge rule).
- No claim may be promoted that identifies a Standard-Model particle, recovers GR as a fundamental theory, or asserts Lorentz invariance is derived.

## Open Questions
- Is C-OG-008 (Christoffel symbols) acceptable as a single claim, or should it be a derived consequence of C-OG-006?
- Is the smooth M->0 limit (C-OG-009) a separate claim, or a property of the affine-geodesic equation already implicit in the textbook einbein reparametrization?
- Does the impedance-matching claim (C-IMP-001) require a corresponding no-reflection Newtonian propagation theorem, or is the rho*Theta = Z_0^2 condition sufficient as declared?
