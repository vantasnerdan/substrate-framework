---
description: 'Coherencev1 mobile-ion screening: conditional retarded drift-diffusion estimate and 100eV gate'
author: giuliano
created: '2026-11-18T00:05:00Z'
updated: '2026-11-18T01:01:54Z'
tags:
- substrate-framework
- research-arc
- coherence
category: efforts
confidence: exploratory
status: active
---

## Objective and Authority
Issue237 and branch research/237-ionic-drift-screening answer Axiom1830's ionic-first discriminator: can mobile D in the actual PdD run-end state, or a standard proton-conducting oxide, add100eV at Ecm0.25..6.5keV? The positive deliverable is a checked conditional closed form and source-bound numbers, not a universal material no-go. No accepted claim is changed; base90323480/releasev0.183.0 contains no imported physical ionic-screening claim. Work shares g5_pairs' single800core-minute whole cap and20pure reserve; no new Native or budget.

## Definitions and Analytic Construction
The preregistered model is an infinite isotropic stable linear drift-diffusion medium and an externally moving unit charge. With number susceptibility chi=dn/dmu_J and collective diffusion D, kappa²=e²chi/(epsilon0 epsilon_r). Source equations imply delta_n=-e chi Dq² phi/(Dq²-iomega), omega=q dot v. Subtract the same background-only potential; the ionic reduction kernel is kappa²/[q²(q²+kappa²-i(v/D)q cos(theta))]. This is a response assumption at small perturbation, not a licensed atomic collision model.

The candidate radial integral is U=(e/(4pi epsilon0 epsilon_r)) kappa asinh(s)/s, s=v/(2Dkappa), with its continuous s=0 value1. To establish it: verify the Fourier elimination and angular rational integral; differentiate J(alpha)=integral_0^infinity atan(alpha q/(q²+kappa²))dq/q, use a*b=kappa²,a-b=alpha partial fractions, establish J'(alpha)=pi/sqrt(alpha²+4kappa²), and anchor J(0)=0. Independent log-wavenumber quadrature of the angular-integrated kernel must agree to1e-8relative. Check static/fast limits and sensitivity to erroneously using the static value at finite speed. All constants are SI, reported energy eV. Deuteron lab speed uses E_lab=2Ecm and the CODATA deuteron mass.

## Discovery, Candidates and Structural Selection
FamilyD uses the sourced diffusive response and earns a cheap conditional number. FamilyI uses an inertial elastic/phonon susceptibility, documented by Rappolt etal DOI10.1140/epjb/s10051-025-00941-2, but has no source-qualified PdD/proton-defect parameters here. Only familyD is executed first under the explicit100eV continuation rule. Published mean-field diffusion excludes molecular vibration/rotation times and strong fields; failure of that applicability is not repaired by a successful integral or called physical exclusion. A passive real electronic background>=1 reduces the ionic contribution in this linear model; its omission supplies a formal best-case estimate, not a quantum bound.

## Claim Delta and Frozen Review
Primary Janssen/Bier arXiv1802.02777v2, SectionII/Eqs1-5 and fast-timescale limitation, and Braun etal arXiv1106.1924v1, hydration Eq2/transport Eq3, are saved/extracted and read. Braun reports3mol%H and cubic reference a4.206A matching sample XRD; the oxide estimate uses the crystallite-local number density0.03/a³ and ideal-carrier susceptibility. RT QENS was assumed immobile; no unmeasured RT mobility or H-to-D mobility is inserted. The oxide static value bounds its finite-speed value only WITHIN this ideal drift-diffusion closure. Literature inputs are cited by URL/DOI and consumed values, without packet attestations.

## Dependency Achievements
The chain is source equations to retarded kernel to checked integral to source-bound PdD/oxide evaluation to100eV allocation decision. Passing algebra licenses the conditional evaluation, not atomic-q applicability or measured fusion yield. The latter still needs the separate electronic response and closest-pair model. A below-trigger result does not establish exhaustion of all ionic mechanisms.

## Implementation and Supervision
Canonical API is substrate_framework.ionic_diffusion. The thin Coherencev1 driver uses the actual7200s PdD/Ecm2.5 profile at301.404969K, its mobile x(z), local metal density from saved layers, pinned853c377 isotherm and Darken mobility. It excludes trapped xt from mobile susceptibility and does not equilibrate v. Source mu' is evaluated with two central-difference relative steps; require relative step change<=1e-5 and positive slope. The source mobility maps D_collective=D_tr x(1-x/xmax)mu'/kT and kappa²=e nM/(epsilon0 mu'_eV). Report top100nm and whole-foil maxima, not a total-H substitution. No new Native starts until this estimate is consumed per1830.

## Attempts and Campaign Checkpoints
The preregistered evaluation completed: result66374bytes SHA6658365e2045f292d7251586f3b1f28e6e0f8930c0ee528931821f2a18d6df5e and proof2445bytes SHA0a6d39c9a541ea76f4ed4f9545fcf7354ff238bdbfad87ef6733ce86a66d6d5e under coherence-pass5/pairs/ionic-first-20261118. Source derivative refinement2.4999e-7, spectral/material relative errors<=3.50e-16/7.67e-15, Darken-Maxwell identity<=5.46e-16. Runtime0.0665+0.4066CPU-s,0Native. Stopped Gamma charged16.1170231833core-min; density readout36.281CPU-s remains inside the20pure reserve.

## Fit, Oracle and Consumer Evidence
The exposing proof checks independent Fourier-kernel quadrature over small/large speed ratios, symbolic elimination and partial fractions, the source Darken flux map and finite-difference stability. A static-at-fast-speed mutation must fail; units and lab/cm speed factor are explicit. One scoped independent review follows the complete deliverable, with one bounded correction. No broad accepted-claim replay.

## Strongest Result and Next Achievement
Within the declared diffusion closure, at0.25keVcm PdD Umax is4.115861e-5eV in the top100nm and0.107213eV anywhere in the foil; both decrease over the requested energy range. BZY10 ideal static model upper is24.19036eV, with no invented RT mobility. Model100eV gate is not reached, not an atomic-response or universal exclusion. Independent IonicReview is complete: full findings array empty, conditional scope supported; original10839bytes SHAe939a52fb2703282d74320e542169f8e8fccd59f9b00235597202f284ee044ea was admitted before acknowledgment. No rerun or correction needed. Focused PR remains before canonical delivery; no ionic Native.

## Promotion, Continued Investment and Operational Stop
A focused PR advances237 after one scoped check/review. Conditional API acceptance and physical applicability remain separate. A below100eV result closes only this model's escalation gate; no deeper ionic Native is launched. Preserve source equations, model assumptions and missing atomic-response construction.

## Cross-References
https://github.com/vantasnerdan/substrate-framework/issues/237; Fleet Coherencev1/g5_pairs entries1830/1832; parent effort coherence-pass5-recomb-pairs-dynamic-2026-11-17; source packet banked privately under coherence-pass5/sources; Native binding obligation236 is separate.
