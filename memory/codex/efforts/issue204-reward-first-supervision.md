---
description: Compress research instructions around reward-first supervision while preserving scientific authority and user control
author: codex
created: '2026-09-07T08:49:20+02:00'
updated: '2026-09-07T08:49:20+02:00'
tags:
- substrate-framework
- workflow
- reward-engineering
category: efforts
confidence: working
status: active
---

## Objective and Boundary

The user requested shorter AGENTS/skills, emphatic achievement/reward framing,
bounded Sol 5.6 high workers, supervisor-led approach/alternative review using
GitNexus or AST plus mathematical reasoning, and PR creation and merge.
Issue #204 predates the PR. Base is main b6fc902a; branch is
workflow/reward-first-supervision. This independent workflow task changes no
scientific APIs, tests, accepted claims/releases, or P253 artifacts. P253 stays
paused; the original worktree's unrelated GitNexus addition remains untouched.
The user explicitly authorized author self-merge of this bounded change; that
is operational permission, not independent scientific acceptance.

## Compression and Semantic Map

The refactor removes repeated policy rather than deleting technical evidence.
Whitespace-delimited word counts compare the original base with this branch;
the counts below exclude this review record. AGENTS owns policy, onboarding owns commands,
skills own specialized methods, and templates collect concrete evidence.

| Surface | Before words | After words | Reduction |
| --- | --- | --- | --- |
| AGENTS and onboarding | 12738 | 3200 | 74.9% |
| Four main skill files | 11011 | 6519 | 40.8% |
| Nine memory templates | 9187 | 4846 | 47.3% |

| Preserved effect | Authoritative location after compression |
| --- | --- |
| Reward before decision; recurring negative-language rewrite; user objective | AGENTS reward and supervision sections |
| Four levels, route verdicts, minimum continuation, certified exhaustion | AGENTS; campaign-proposal manifest keys unchanged |
| Ten achievements, authority, review/merge/completion distinction | AGENTS; governance reference |
| Bounded workers, actual-context provenance, early approach/alternative review | AGENTS; subagent-task; skill execution steps |
| Typed object/ensemble/space/observable; full variations; analytic remainder | problem-deconstruction; physics skill |
| Exact vs numerical/formal evidence, sensitivity, FFT, independent maps | oracles reference; physics skill |
| Numerical floors, full error budget, zero modes, crossed refinement, precision | small-ratio skill (technical methods retained) |
| One substantive review plus bounded correction; strongest useful scope | review templates; harvest and synthesis skills |
| Issue-before-PR, external harvest, terminal-close, user-authorized self-merge | AGENTS; harvest; onboarding |
| Proportionate replay, first-run receipts, immutable history, explicit pause | AGENTS; physics skill; effort/arc templates |

Two precision repairs preserve the intended evidence boundary: an ordinary
script proves only its evaluated predicates, not an entire prose derivation;
a negative Hessian direction refutes an energy minimum, while dynamical
instability needs its separate operator/dynamics bridge. These are workflow
interpretation corrections, not new accepted scientific claims.

## Structural and Behavioral Review

The registered GitNexus graph belongs to another worktree and was stale, so the
supervisor inspected the exact local AST/source instead. choose_validation_scope
maps these policy/skill/template paths to test_public_contribution_surfaces.py
and test_repository_validation.py. AST inspection traced those functions'
file-read/assertion consumers; relative Markdown links resolved. No changed
scientific symbol/caller graph exists in this documentation-only diff.
Graph/AST evidence does not certify mathematical completeness.

A bounded read-only Sol 5.6 high review was requested for AGENTS, onboarding,
physics/harvest/synthesis skills and the worker template. Six scenarios checked:
compatibility abort versus scientific refutation; lemma versus parent completion;
explicit pause versus persistence; self-merge versus scientific independence;
worker ansatz/objective drift; stale graph versus mathematical proof.
The review found two inconsistencies, both repaired: supervisor stewardship
rather than ownership of the user's objective, and user-or-owner authorization
for self-merge. Worker PR/registry return ownership was also made explicit.
This scenario review supports semantic preservation, not a measured claim about
future agent performance or token-pressure resistance.

## Validation and Next Achievement

The local scope selector chose fixed checks plus test_public_contribution_surfaces.py
and test_repository_validation.py. Fixed repository/generated/memory/skill/import/
compile checks passed (1066 memory files valid; 43 existing advisory warnings).
The first pytest run had 17 passes and one literal-text failure: a line wrap split
"strongest meaningful" in AGENTS. The phrase was joined, preserving its meaning
and the existing test. Only the affected contribution-surface tests were replayed:
11 passed, alongside the unchanged seven repository-validation tests from the
first run. All 18 affected tests pass at the final boundary; no tests were modified.
The new memory record validates without warnings; git diff --check passes.

Commands were run from the isolated worktree, using the existing interpreter:

```sh
PYTHON=/home/dan/substrate-framework/.venv/bin/python scripts/validate.sh --pytest-scope tests/test_public_contribution_surfaces.py tests/test_repository_validation.py
PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python -m pytest -q tests/test_public_contribution_surfaces.py
memory validate --base /home/dan/substrate-workflow-rewards /home/dan/substrate-workflow-rewards/memory/codex/efforts/issue204-reward-first-supervision.md
git diff --check
```

The first command's literal-text failure and bounded repair are preserved above,
not relabeled as an initial passing suite. The bounded reviewer correction check
also passed. Final PR/landed commit and closure are tracked in issue #204.
Unchanged scientific oracles remain untouched; P253 requires user resumption.

## Cross-References

Canonical issue: https://github.com/vantasnerdan/substrate-framework/issues/204
Paused campaign: https://github.com/vantasnerdan/substrate-framework/issues/203
