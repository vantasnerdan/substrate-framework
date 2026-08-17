---
description: 'Process rule learned from PR #63 closure: prefer focused 2-file harvest PRs over broad 10+ file bundles that mix campaign narrative, validation artifacts, and reusable code. Place narrative, dimensional reviews, and campaign state in closed PR history, not in the merge unit. Cite the canonical issue; never self-merge.'
author: vantasner
created: '2026-08-14T13:35:00Z'
updated: '2026-08-14T13:35:00Z'
tags:
- substrate-framework
- process-rule
- harvest-style
- pr-scope
- non-self-merge
category: decisions
confidence: established
status: active
---
# PR scope rule: focused harvest preferred

## Rule statement

When opening a PR that builds on an existing campaign, proposal, or harvest, prefer a focused 2-file PR that contains the smallest independently correct and reusable unit. Do not bundle the campaign narrative, dimensional review, attempt manifest, evidence placeholder, adjudication file, and memory updates into the same PR.

## Rationale

A focused harvest PR is easier to review, merge, and audit. The reviewer can confirm the merge unit in isolation. The bundled artifacts (campaign narrative, dimensional review, etc.) belong in the closed PR history as provenance, not in the merge unit. This separation is consistent with AGENTS.md's distinction between artifact merge and claim promotion.

## What goes in the focused PR

- The reusable module or implementation
- The test suite that exercises the module
- (Optional) a small CI-friendly verifier script that can run without host-specific paths

## What stays in the closed PR history

- Campaign proposals and adjudication
- Attempt manifests and result files
- Evidence and review placeholder directories
- Memory entries that summarize the campaign
- Dimensional review scripts whose conclusions are not the merge's load-bearing tests
- Broad narrative documents

## Non-self-merge rule still applies

Even with a focused scope, the agent that authored the commits on the PR cannot merge the PR. The merge must come from a distinct agent or the repository owner. The focused scope makes the handoff cleaner because the reviewer has fewer files to verify.

## Source of this rule

PR #63 (broad 11-file P226 source PR) was closed unmerged after a harvest-review identified six blockers. PR #64 (focused 2-file harvest of the same work) addresses all six blockers and is ready for review. The contrast demonstrates that the focused scope is the right default.

## Cross-references

- Closed PR: https://github.com/vantasnerdan/substrate-framework/pull/63
- Focused PR: https://github.com/vantasnerdan/substrate-framework/pull/64
- Lessons memory: memory/vantasner/decisions/P226-pr63-closure-lessons.md
- AGENTS.md non-self-merge rule: see "An agent must never merge a PR that it opened, authored a commit for, or materially implemented"
