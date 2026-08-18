---
description: Deliver the exact-mass one-loop proper-time determinant coefficients with preregistered regulator comparison as the next Route 1 (induced-gravity) rung for issue #76
author: prime-agent
created: '2026-08-18T00:00:00Z'
updated: '2026-08-18T00:00:00Z'
tags:
- substrate-framework
- effort
- route1-induced-gravity
- issue-76
category: efforts
confidence: working
status: active
---

## Goal and Success Contract
This effort delivers one importable, target-blind reusable unit for the #76 goal: the exact-mass one-loop proper-time coefficients of the accepted Route 1 fluctuation operator, compared across preregistered regulators. It is complete for this rung when a new module under `src/substrate_framework/` returns closed forms for the tau^-2 class (curvature / inverse-Newton-matching integral) and the tau^-3 class (vacuum/cosmological integral) for three declared regulators (sharp proper-time cutoff; smooth essential-singularity proper-time weight; power-divergence-subtracted zeta scheme), each cross-checked against independent numerical quadrature; when the sharp massless limit reproduces the accepted leading coefficient `s = N*(1-6*xi)/(12*pi)` from `scalar_induced_newton` exactly; when the exact-mass sharp factor `exp(-z) - z*E1(z)` named load-bearing by the PR #13 review is exposed as a derived, mutation-sensitive API rather than a footnote; and when the scheme spread among the three regulators is returned as an exact quotable ledger. It does NOT identify a cutoff with a substrate scale, choose a renormalization condition, derive a total Newton constant, produce a sourced nonflat solution, compose worldline probes, or open any empirical comparator; those are later rungs and open campaign frontier for #76. A no-go, residual, or obstruction is attempt evidence, not completion.

## Accepted Baseline
Work starts from accepted release `v0.160.0`, framework commit `1b00c3a` (branch `research/p230-exact-mass-regulator-rung` off updated `origin/main`, which contains the P227 closure `2072648`). Source audited at this commit (with landed-conditional status noted per the PR #77 review): `AGENTS.md`; `.agents/skills/physics-erdos-loop/SKILL.md` and its `references/governance.md` and `references/oracles.md`; `.agents/skills/research-pr-harvest/SKILL.md`; `tools/agent-memory/skills/shared_agent-memory-usage/SKILL.md` (CLI sections); `governance/releases/current.yaml` (207 accepted claims); `governance/claims.yaml` entries `C-GRV-001`, `C-GOR-001/002`, `C-STG-001/002`, `C-WLN-001..003`, `C-LOR-001/002`, `C-GW-001..010`; issue #76 body; issues #9 and #12 including all four comments each; `src/substrate_framework/scalar_induced_newton.py` and `src/substrate_framework/covariant_sine_gordon_action.py` (both landed conditional unpromoted APIs; the PR #14 and PR #25 harvests promoted no claims); `src/substrate_framework/einstein_scalar.py` docstring conventions; `src/substrate_framework/governance.py` (proposal schema). Chronology and memory prose are not authority; every reused fact was verified at source.

## Constraints and Invariants
Conventions are exactly the landed conditional `scalar_induced_newton` ones (accepted C-GRV-001 supplies only the conditional dimensional and additive-baseline ledger): `D_E = -nabla_E^2 + xi R_E + m^2`, boundaryless Gilkey-Seeley-DeWitt weights `a_2 = ((1/6 - xi) R_E - m^2)` under the `(4*pi)^-2` prefactor, one-real-scalar determinant weight `1/2`, and the Euclidean Einstein-Hilbert matching factor `16*pi`. The mass is retained exactly: no expansion of `exp(-tau*m^2)` inside the regulated proper-time integral. Target-blindness is a tested invariant: no observed `G`, `M_Pl`, or benchmark value may appear; the cutoff stays a declared formal scale; the power-subtracted scheme requires an explicit `renormalization_scale` argument (regulator non-defaulting, mirroring the accepted module). No new fitted constant. Permitted imports: `sympy` (including `expint`/`E1` and Bessel `K` closed forms), `scipy` quadrature in tests only, framework `exact_symbolic`, `scalar_induced_newton`, and the declared composition slot `covariant_sine_gordon_action` (`m^2 -> V''(phi_bg)`). Write boundary: `src/substrate_framework/scalar_one_loop_mass.py`, its `__init__` exports, `tests/test_scalar_one_loop_mass.py`, `proposals/P230-exact-mass-regulator-rung/`, and this record. Author-does-not-merge: open a PR naming #76 with `Advances #76`; a distinct reviewer/owner harvests and merges.

## Decomposition
1. [x] Recall and source verification (authority, #76 gates, #9 no-go, #12 ladder and PR #13/#25 handoffs, accepted rung conventions).
2. [x] Candidate and selection-criteria preregistration in `proposals/P230-exact-mass-regulator-rung/proposal.yaml` (validated WORKFLOW VALID at `1b00c3a`).
3. [x] Coordination comment on #76 (branch, proposal id, claim ids after collision search, slice, write surfaces, comparator-blinding point).
4. [x] Importable implementation of the closed forms for the tau^-2 and tau^-3 classes under all three regulators, with exact symbolic APIs.
5. [x] Verifier and sensitivity audit: quadrature cross-check of every closed form, massless-sharp reduction to the accepted `s*Lambda^2`, small-z asymptotic recovery, mutations (xi sign, prefactor, special-function branch, unknown regulator, non-defaulted scale), and scheme-spread ledger.
6. [x] Framework-fit and downstream replay (`scripts/validate.sh --full`: 2235 passed, ALL REPOSITORY WORKFLOW CHECKS PASS; `git diff --check` clean; GitNexus re-indexed at current commit, impact LOW, detect-changes maps only the new module and `__init__` exports; generated AGENTS.md injection from `gitnexus analyze` reverted per AGENTS_START_HERE) (`scripts/validate.sh` boundary appropriate to the diff; `git diff --check` separately).
7. [~] Harvest PR advancing #76 opened as #77 (issue-first, `Advances #76`, author does not merge). Independent harvest review at head `04f26f2`: CHANGES REQUIRED (three blocking findings, all verified independently and repaired in attempt 0002 — vacuum double-count, zeta sign semantics, authority wording). Repaired head pushed; re-review requested; merge authority remains with the distinct reviewer.

## Candidate Preregistration and Selection Criteria
Selection criteria frozen before any comparison value is computed: (1) convention compatibility with the accepted rung APIs; (2) closed-form exactness in standard special functions with no fitted constants; (3) exact-mass integrity inside the regulated integral; (4) mutation sensitivity; (5) reuse value for the later cutoff-identification and renormalization-condition rungs; (6) regulator explicitness (scheme dependence stays visible).

- Candidate A (sharp proper-time, exact mass): regulated integral `integral_{1/Lambda^2}^infty tau^-2 exp(-m^2 tau) dtau = Lambda^2*(exp(-z) - z*E1(z))`, `z = m^2/Lambda^2`. Directly discharges the PR #13 load-bearing finite-mass factor.
- Candidate B (smooth essential-singularity weight `exp(-1/(Lambda^2 tau))`, exact mass): `integral_0^infty tau^-2 exp(-m^2 tau - 1/(Lambda^2 tau)) dtau = 2*Lambda^2*sqrt(z)*BesselK_1(2*sqrt(z))`. One convergent integral family for all weights; no hard step.
- Candidate C (power-divergence-subtracted zeta/Mellin finite part with declared scale `mu`): `m^2*(ln(m^2/mu^2) + gamma_E - 1)`. Exhibits the logarithmic-running alternative in which no `Lambda^2` term survives.
All three are retained as parallel schemes in the module (the rung's deliverable is the comparison itself, so the selection is "which schemes to expose", not "which one wins"); empirical comparators may not select among them.

## Attempts
Attempts are append-only and individually reproducible; each row names the diagnosed mechanism and the next materially different attempt.

| Attempt | Candidate or repair | Artifact and command | Verdict | Mechanism | Next attempt |
| --- | --- | --- | --- | --- | --- |
| (pending) |  |  |  |  |  |

## Validation
Validation covers the actual objective, verifier sensitivity, limits, conventions, and dependency replay, not merely an exit code. Targeted oracle at the repaired head: `PYTHONPATH=src .venv/bin/python -m pytest tests/test_scalar_one_loop_mass.py -q` = 25 passed (5.84 s), including the review-driven additions: the mass-resummed determinant-integrand oracle, the d I_3/d m^2 = -I_2 derivative-identity oracle in all three schemes, the double-count rejection mutation, and both zeta sign counterexamples with undecidable and positive branches. Earlier 20/20 boundary at the first head (5.20 s). Full boundary: `bash scripts/validate.sh --full` = 2235 passed in 309.16 s, `ALL REPOSITORY WORKFLOW CHECKS PASS (full pytest suite)`; `git diff --check` clean in a separate invocation. GitNexus: index refreshed at the working commit; `impact leading_scalar_newton_shift_coefficient` = LOW risk, single new direct caller (`exact_mass_inverse_newton_shift`); `detect-changes` maps only `src/substrate_framework/scalar_one_loop_mass.py` (new) and `__init__.py` `__all__`/imports; zero affected execution flows. The `gitnexus analyze` run auto-appended a tracked-instruction block to `AGENTS.md`; it was reverted per AGENTS_START_HERE section 3 and not committed.

## Debt Ledger
Every new assumption, import, parameter, residual, broken consumer, or narrative inconsistency is recorded here and discharged; this table must be empty at close.

| Debt | Introduced by | Why it is real | Discharge artifact | Status |
| --- | --- | --- | --- | --- |
| (none yet) |  |  |  |  |
