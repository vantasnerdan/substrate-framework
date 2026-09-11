---
description: "Referee verdict on the 2026-09-11 Newton sign-flip announcement (d186 comment 18406566, rev 294): no sign was flipped; form-level algebra verified, energy-level claim pending OpenWave R19-1"
author: giuliano
created: '2026-09-11'
updated: '2026-09-11'
confidence: established
status: active
category: decisions
tags: [referee-audit, P252, newton-sign, so13, fermez-R19]
---

## Claim Under Review

The audit refereed the announcement "Finally flipped Newton sign, finite nonzero frequency!" (discussion #186, comment 18406566, report rev 294), checking its form-level claims against independent oracles in PR #212 (branch research/P252-newton-sign-external-audit, canonical issue #211).

- Final tally: ALL 37 CHECKS PASS, ALL 13 MUTATIONS BREAK (exact SymPy + mpmath; static-source derivation plus independent field-equation route).
- Verified: sector symmetries, bracket spin content, the signature flip (signed structure triples s(K,K) = -s(J,J) in so(1,3), + in so(4); Maurer-Cartan field check), exchange sign rule, trace/eps identities, multipole hierarchy 1/d, 1/d^3, 1/d^5, composition numbers (Z^2 liquid-drop, better than 1 percent).

## Verdict and Open Items

The subject line is contradicted by the report's own content: rev 294 states the sign was never wrong and that range, universality, and strength still fail.

- The energy-level static Newton read is pre-registered as OpenWave R19-1 (gates NEWTON_SIGN_REVERSED / CANDIDATE_REFUTED); the flip is established at constraint level only.
- Debts: inc-side index convention and the 2mp counting convention are unpinned from the comment (need Zenodo 22714918); the 48^3 lattice run is provenance_only.
- "Finite nonzero frequency" reduces to omega = K/I with K topologically fixed (massless clock) - consistent with the certified P249/P250 clock structure, no new mechanism.
