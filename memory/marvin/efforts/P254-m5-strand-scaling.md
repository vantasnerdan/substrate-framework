---
description: Conditional straight-strand tension theorem comparing two distinct M5 potentials
author: marvin
created: '2026-09-25T16:40:29+00:00'
updated: '2026-09-25T17:20:48+00:00'
tags:
- substrate-framework
- research-arc
- M5
category: efforts
confidence: working
status: active
---

## Objective and Authority

The bounded result for issue #218 is an independently derived theorem for a regular straight transverse-pair radial strand, relating potential order to tension scaling. Release v0.183.0 at source commit 90323480 accepts a distinct auxiliary-frame M5 action; P252/issue #211 inventories older discussion #186 claims without promoting the later R25 or Report 018. The 61 top-level comments from 18303935 through 18601274 (plus thread replies) form the parent harvest, not 61 accepted claims.

## Definitions and Analytic Construction

The two frozen-slot ansatzes have the same transverse matrix shape but different contractions and potentials. Write q=b², s=rho²/2 and derive the matrix curvature and the one-dimensional energy independently for each action; then complete a positive square and integrate between q(0)=0 and q(infinity)=b0². A bound proved on the fixed-axis chart licenses neither arbitrary slot relaxation nor a stable 3D strand or pair force. No numerically irreducible proposition is needed to establish this conditional analytic bound; lattice rows remain applicability evidence only.

## Discovery, Candidates and Structural Selection

The motivating tension is that OpenWave R25 reports a delta-fourth-power line cost while Mikulski Report 018 reports a beta-cubic cost. Both can be true under different potentials; compare the actual matrix density and potential before comparing fitted energies. Candidate A is the exact q,s square completion for a generic positive potential specialized to both actions. This is a fixed theorem, not an open mechanism search. A different, unbuilt route—positive two-derivative doublet stiffness with a mass—would require its own model and asymptotic-source audit and is not a shortcut to Newton attraction.

## Claim Delta and Frozen Review

C-M5L-001 now has a scoped accepted review and an on-branch registry/release transaction, with no accepted dependency on the two external actions. It is a conditional compatibility extension with frozen slots, a separate constructive R25 director counterexample, and source exposure declared before P254 freeze; no fresh predictive numerical claim or blinding is asserted. `proposals/P254-m5-strand-scaling/reviews/C-M5L-001.md` records the initial false assertion, governance errors and independent correction checks.

## Dependency Achievements

The regular core and fixed slots license a sharp radial infimum, not a full static lower bound: an R25 diagonal director-gradient commutes with the pair curvature and a positive annular bump strictly lowers its trace potential. The winding-two equality profile is non-C1 at the axis, but smooth cores recover its infimum with O(epsilon²) excess. These proofs do not certify natural localization, uniaxial mass, a stationary ring or charge-pair confinement. R25 charge-strand rows remain mostly unconverged; Report 018 pair charges melt for tested pinned flows.

## Implementation and Supervision

The importable `m5_strand` SymPy API and six focused independent regressions own matrix curvature, trace potential, sharp bound and director variation; no third-party code is copied. Three read-only scouts split the discussion into R14–R19, R20–R24 and R25–latest, and two independent reviewers inspected the mathematical and governance slices. Initial math review found a real error about director-gradient cost; initial governance review found winding/exterior/target-kind mismatches. Both were repaired, and each reviewer reported a correct scoped correction check without running tests. Their read-only reviews are not merger authorization.

## Attempts and Campaign Checkpoints

Initial external executions: Report 018 `strand_bogomolny.py` at 4823e5e exited 0, matrix density True and printed beta=.3,k=.5 bound 0.03998594644342529; OpenWave R25-0 `m5_32_r25_0_form.py run` at a73f6fc exited 0 with all its `a,a2,b,d,e,f` groups reporting PASS; Report 017 `symbolic.py` exited 0 with velocity degrees [0,2], trace blindness True/True and principal-symbol identity True. The independent exploratory matrix probe first failed an assertion comparing algebraically equal expressions structurally, then passed using a simplified difference; at beta=.3 a floating director slot lowered the local trace-power potential from 0.00248529516 to 0.00142035015. These are exploratory receipts, not a proof of a full-field minimum. The proposal's first schema validation exposed a missing v2 `exhaustion_certificate` key; the corrected proposal passed `scripts/validate_repository.py` (271 existing accepted claims, 15 proposals).

## Fit, Oracle and Consumer Evidence

The six-test `tests/test_m5_strand.py` suite passes and the affected `tests/test_governance.py`, `tests/test_repository_validation.py` and `tests/test_public_contribution_surfaces.py` pass with it: **42 scoped tests** under `scripts/validate.sh --pytest-scope ...`, including fixed repository/generated/memory/skill checks (272 accepted claims; 1072 memory entries valid with 44 warnings). Its independent 4x4 mixed-matrix trace, 2x2 curvature, spectral matrix potential and free-director commutator oracles catch wrong normalization, regular-core boundary loss and potential stiffness freezing. A throwaway mutation probe returned exact residual `32 b²(b')²/rho²` for half norm, `9/16` for a free core and `1.00369078` variable/frozen stiffness at delta=.3. Installed `.venv` without `PYTHONPATH` loaded the new API and staged release/memory. A full run was started then canceled on user direction before exit; only the scoped completed run is evidence. The per-route receipt is `proposals/P254-m5-strand-scaling/attempts/0001/receipt.md`; generated v0.184.0 consumers are staged on the branch, not deployed.

## Strongest Result and Next Achievement

The strongest supported result is the dual-potential exact radial infimum plus a constructive same-R25-action lower-energy static competitor. It corrects both the apparent cubic/quartic conflict and the full-winding smooth-attainment overclaim; `PYTHONPATH=src .venv/bin/python -m pytest tests/test_m5_strand.py -q` passed 6 tests, and the source scripts gave the narrower outputs recorded above. Newton, clock and pair-confinement mechanisms remain open as issue #211's separate frontier.

## Promotion, Continuation and Pause

On the branch, P254 is terminal_success/complete and C-M5L-001 review, registry, v0.184.0/current pointer, generated claim index and accepted framework memory agree. The 42 scoped affected tests and installed consumer smoke passed; a full suite was canceled without a verdict because its scope and runtime were not warranted by this change. The review PR's remote state still requires verification; approval, merge, installation and physical acceptance are not implied. Issue #211 remains the separate broader live harvest, and the requested PR must remain open for a distinct merger.

## Cross-References

Primary sources: discussion #186 comments 18581233, 18583989, 18601274; OpenWave a73f6fcd71d6b45bf479e76a3999123b39e6d9d6; Mikulski 4823e5ec2368420f28229ac73067108ec01e6e23; `proposals/P254-m5-strand-scaling/proposal.yaml`; issue #218; P252/issue #211.
