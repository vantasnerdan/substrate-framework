---
description: 'Lessons learned from PR #63 closure (P226 substrate_metric broad PR). Six blockers identified in the harvest review: misread paper conventions, paired-wrong-definition gates, tests that didn''t test the claim, sys.path host-injection, docstring/implementation mismatch, and over-broad PR scope. PR #64 is the focused harvest that addresses all six.'
author: vantasner
created: '2026-08-14T13:30:00Z'
updated: '2026-08-14T13:30:00Z'
tags:
- substrate-framework
- process-lessons
- harvest-review
- verifier-design
- dimensional-analysis
- pr-scope
- campaign-P226
- self-optimization
category: decisions
confidence: established
status: active
---
# Lessons from PR #63 closure (P226 substrate_metric harvest review)

## Source

PR #63 (https://github.com/vantasnerdan/substrate-framework/pull/63, "P226 attempt 0001: substrate_metric module + overlap audit (Advances #62)") was closed unmerged after a harvest-review comment identified six blocking findings. The exact head branch `p226-fulceri-sg-einbein-review` is preserved as closed-unmerged research provenance. PR #64 (https://github.com/vantasnerdan/substrate-framework/pull/64, "Selective harvest: conditional flowing 1+1 metric (Advances #62)") is the focused 2-file harvest that addresses all six blockers.

## Six blockers and the lessons

The harvest review identified six distinct blocking findings. Each one is recorded below with the specific failure mode and the corrected approach.

### 1. Read paper definitions before dimensional analysis

The dimensional validation in `dimensional_validation.py` misread three paper conventions. The lesson is that dimensional review must extract definitions from the paper text first, then check formulas against those definitions.

The paper's equation (1) defines `mu0^2` as the on-site potential coefficient in an energy-per-length Lagrangian density, and equation (5) defines `omega0 = mu0/sqrt(rho0)`. Under those definitions plus C-MED-003, `mu0^2 * c0 / omega0 = sqrt(Theta_0 * mu0^2)` is the correct energy scale. The D1 gate in the closed PR instead hardcoded `mu0` as 1/length and reported a false dimensional inconsistency.

The paper's equation (14) prints `P = gamma * (v/c0^2) * E0`. The closed PR's D2 gate tested `v/c0` instead of `v/c0^2` and reported a paper typo that is not present.

For coordinates with `[t] = T` and `[x] = L`, the metric components `g_tt`, `g_tX`, `g_XX` necessarily have different component units; multiplying by `dt^2`, `dt dx`, `dx^2` makes every term of `ds^2` have units `L^2`. The closed PR's D6 reported this as a metric inconsistency when it is in fact a coordinate-units property.

The dimensional_validation script returned status zero despite printing D1/D2/D6 failures, so automation treated the failed review as passing. Scripts that print "FAIL" must not return zero.

### 2. Verifier gates must not pair wrong definitions

The closed PR's G2 used `gamma = 1/sqrt(1 - v^2 * c0^2)` and `P = gamma * E0 * v`. The paper uses `gamma = 1/sqrt(1 - v^2/c0^2)` and `P = gamma * E0 * v/c0^2`. The residual closes only because the wrong definitions are paired. The lesson is that any pair of test expressions must be checked against an independent reference, not against each other. AGENTS.md already names this as the "PASS through double error" anti-pattern.

The closed PR's G3 solved for a free `h_eff` to fit the desired relation instead of checking the paper's independently defined effective action scale. Tests must encode the claim's stated definitions, not solve for parameters that make the test pass.

### 3. Tests must test the claimed claim

Four tests in the closed PR did not exercise the property they claimed to test. G10 compared two literals to themselves rather than comparing against `optical_metric_1d`. G13 passed a coordinate-independent symbol `n`, so every derivative is zero and the test never exercised the claimed nontrivial static connection. G8 checked only that the quadratic has two roots, not that they equal either sign convention. G12 / G13 injected `/home/dan/substrate-framework/src` into `sys.path`, making the campaign verifier host-specific.

The lesson is that test names must be auditable: the assertion should fail under any mutation that violates the claimed property, and the inputs must be general enough that the test actually exercises the non-trivial case.

### 4. Module docstrings must match implementation

The closed PR's module docstring claimed `g^XX = V^2 - n*V^2/c0^2` while the implementation computed the correct `g^XX = 1/n - n*V^2/c0^2`. Docstrings and code must agree, and the public documentation should expose the actual closed-form result.

### 5. PR scope must be focused

The closed PR included 11 files: the campaign proposal, adjudication, attempts manifest, evidence placeholder, reviews placeholder, effort memory, dimensional validation, verifier script, module, and test. The reviewer's recommendation was to harvest only the smallest independently correct and reusable unit. The focused PR (#64) is 2 files: the module and the test. The campaign narrative, dimensional findings, and impedance product were left in the closed PR's history.

### 6. Use precise authority language

The closed PR's module header did not distinguish between "conditional unpromoted infrastructure" and "candidate claim". PR #64's header explicitly states "Authority status: conditional, unpromoted infrastructure linked to issue #62 and source PR #63". The body of PR #64 also explicitly excludes any claim promotion.

## What the focused PR (#64) does right

The focused harvest addresses all six blockers through a narrower scope and corrected verifiers. The key properties of PR #64 are listed below.

- 2-file scope: module + test
- 9 explicit tests covering inverse, static rescaling, Minkowski limit, general Christoffels against an INDEPENDENT reconstruction (not just the package helper), static connection match, load-bearing mutations (timelike deletion, off-diagonal sign flip), input guards for 9 invalid input types, chart-helper constant-c0 enforcement, and a public-API surface check
- Massless speeds are explicitly `-V - c0/n, -V + c0/n` (corrected sign from the paper)
- Static limit test asserts `g_source(V=0) = c0^2 * optical_metric_1d(n, c0)` (the actual rescaling)
- `__all__` is explicit and narrow (4 symbols)
- `impedance_product` is dropped (trivial, not a reusable constitutive construction)
- `dimensional_validation.py` is dropped (misread paper)
- Campaign / adjudication / evidence / review placeholders are dropped (no durable reviewed object)

## Cross-references

The PRs, comments, and campaign artifacts referenced by these lessons are listed below.

- Closed PR: https://github.com/vantasnerdan/substrate-framework/pull/63
- Focused PR: https://github.com/vantasnerdan/substrate-framework/pull/64
- Review comment: https://github.com/vantasnerdan/substrate-framework/pull/63#issuecomment-5293561189 (approx, harvest review at 13:19:27Z)
- Issue #62: https://github.com/vantasnerdan/substrate-framework/issues/62
- Campaign: P226 in campaigns/P226-fulceri-sg-breather-einbein-effective-metric/

## Generalization

For future paper reviews and campaign work:

1. Before any dimensional analysis, extract the paper's variable definitions from the source PDF (especially Lagrangian densities, energy density definitions, constitutive maps). Do not assume SI conventions.
2. For any verifier gate that pairs two test expressions, check at least one of them against an independent reference. If both are needed only to cancel each other in a residual, the gate is not testing the claim.
3. Test inputs must be general enough to exercise the non-trivial case. Coordinate-independent symbols produce zero derivatives and test nothing.
4. Docstrings must match implementation. Closed-form algebraic results should appear in the docstring AND the code, with the same expression.
5. PR scope: prefer 2-file focused harvests over broad 10+ file bundles. Place narrative, dimensional reviews, and campaign state in the closed PR history, not in the merge unit.
6. Authority language must be precise: "conditional unpromoted infrastructure" is not the same as "candidate claim".
