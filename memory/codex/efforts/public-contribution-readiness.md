---
description: Prepare substrate-framework for safe public contributions with explicit rights, CI, security policy, and protected-main governance
author: codex-public-readiness
created: '2026-08-18T09:40:10+02:00'
updated: '2026-08-18T11:10:00+02:00'
tags:
- substrate-framework
- effort
- public-contributions
- repository-security
category: efforts
confidence: working
status: active
---

## Goal and Success Contract
This effort delivers a safely public, contribution-ready `vantasnerdan/substrate-framework` repository under the owner-selected Apache-2.0 license, with reviewable contribution/security/conduct policies, least-privilege pull-request CI, dependency maintenance, enabled GitHub security features, and a protected `main` branch. It is complete only after every public reachable ref is free of unredistributable or secret material, the repository becomes public, the protection and security settings are verified through GitHub APIs, a fresh anonymous clone passes the declared checks, and the debt ledger is empty. A readiness report, blocked visibility change, or unmerged documentation PR is progress rather than completion.

## Accepted Baseline
Work starts from `main` commit `1b00c3a` and accepted release `v0.160.0`. The normative local sources are `AGENTS.md`, `AGENTS_START_HERE.md`, the existing issue and pull-request templates, GitHub repository settings, and exact external rights metadata. Zenodo record `10.5281/zenodo.21879560` states that `incoming/einbein_1plus1D_tutorial.pdf` is open access under `GPL-3.0-or-later`. Crossref metadata for `10.1007/BF02833896` supplies Springer text-and-data-mining terms, not verified redistribution permission for the publisher PDF or its full extracted text.

## Constraints and Invariants
Publication must not expose secrets or material lacking redistribution permission. The user selected Apache-2.0 and designated `vantasnerdan`, `axis-marbell`, and `mlops-kelvin` as the only merge-authorized maintainers whose CODEOWNER approval may satisfy the protected-main review gate. The user has not authorized a history rewrite/force-push; that materially different owner decision remains gated. Existing accepted scientific authority, immutable campaigns, generated documentation, and PR #77 are outside the implementation write boundary. Every file-change PR requires issue #78, uses `Advances #78` while publication remains incomplete, and is merged only by a distinct reviewer/owner. The original worktree's untracked `.claude/` and `CLAUDE.md` are preserved.

## Decomposition
Work proceeds through these dependency-ordered steps and continues after failed attempts.

1. [x] Audit GitHub settings, community profile, all refs, contributor metadata, bundled documents, secrets, and rights provenance.
2. [x] Create canonical issue #78 and claim an isolated branch/write boundary.
3. [ ] Add contribution, conduct, security, CI, dependency, and scanner configuration without changing scientific authority.
4. [ ] Validate the documentation/workflow boundary and open a non-self-merged PR.
5. [ ] Land the owner-selected Apache-2.0 license and obtain the remaining history-sanitization decision; execute the approved rights-safe strategy.
6. [ ] Make the repository public, apply protected-main/security settings, and verify through a fresh anonymous clone.

## Publication Strategy Alternatives
Selection is blocked pending explicit owner authority because the alternatives have materially different provenance and collaboration effects.

| Candidate | Construction | Benefit | Cost or blocker | Selection evidence |
| --- | --- | --- | --- | --- |
| A | Rewrite this repository's reachable history to purge the Preparata publisher PDF/full extraction, repair durable provenance references, and force-push the verified replacement refs | Preserves the requested repository URL and full sanitized history | Changes commit IDs, invalidates pinned source hashes and collaborator clones, disrupts open PR #77, and requires explicit owner authorization plus a migration map | Fresh-clone all-ref absence proof, registry/release/provenance replay, contributor coordination |
| B | Create a new sanitized public repository from an approved snapshot with a provenance manifest while this repository stays private | Avoids destructive rewriting of the governed private history | Does not make the exact requested repository public and splits discovery/provenance surfaces | Owner approval of a public mirror and explicit private/public authority boundary |
| C | Keep this repository private until redistribution permission is documented for the source artifacts | Preserves current history and scientific pins exactly | Does not satisfy the requested public outcome unless permission is obtained | Publisher/rights-holder license or written permission covering redistribution |

Selection criteria are rights safety, preservation of scientific provenance, impact on active contributors/open PRs, auditability, and satisfaction of the requested repository URL. Convenience cannot override redistribution rights or the non-force-push boundary.

## Attempts
Attempts are append-only and individually reproducible.

| Attempt | Candidate or repair | Artifact and command | Verdict | Mechanism | Next attempt |
| --- | --- | --- | --- | --- | --- |
| 0001 | Redacted all-history and filesystem secret scan | `ghcr.io/gitleaks/gitleaks:v8.30.1 detect` on 706 commits/all refs and the detached main worktree | Qualified pass | One `generic-api-key` hit is a plain 28-character identifier list under YAML key `public_api`, with no URL/assignment/credential structure; all other rules clean | Add a narrow scanner allowlist with documented false-positive provenance and run CI scan |
| 0002 | Rights and community-surface audit | GitHub API community/settings queries, file inventory, Zenodo/Crossref metadata | Blocked publication | No project license; paywalled Preparata PDF/full extraction are reachable; CI/security/community files and branch protection are absent | Implement reversible readiness files, then obtain owner legal/history choices before visibility mutation |

## Validation
- Rights/secret oracle: redacted Gitleaks all-ref history scan plus filesystem/archive scan; exact Zenodo/Crossref metadata and tracked-source inventory.
- Workflow safety: static inspection for least-privilege permissions, untrusted fork execution, mutable action tags, and secret access; GitHub Actions execution after PR creation.
- Dependency replay: repository fixed checks and contribution documentation link/config validation.
- Targeted tests during implementation: documentation/workflow parsers plus the smallest justified repository test scope.
- Final repository validation: `scripts/validate.sh --full` because the final diff changes package/build metadata and the public CI boundary.
- `git diff --check` runs separately at the unchanged PR boundary.
- Public-state verification: GitHub APIs for visibility/rules/security plus a fresh anonymous clone and all-ref rights/secret scan after the approved history transaction.

## Debt Ledger
Every unresolved item inside the requested public outcome remains debt until discharged.

| Debt | Introduced by | Why it is real | Discharge artifact | Status |
| --- | --- | --- | --- | --- |
| Owner-selected project and documentation license absent | Repository baseline | Public availability alone grants no reusable software rights and cannot define contribution licensing | Apache-2.0 `LICENSE`, package metadata, contribution terms, and third-party notices | in progress |
| Preparata publisher PDF and full extracted text reachable in history without verified redistribution permission | Commit `81df095` and descendants | Public Git history would redistribute the complete paywalled work | Written permission or owner-approved verified purge/mirror strategy | open |
| History rewrite/mirror authority absent | Safety boundary | Purging all refs is destructive, changes commit IDs, and affects open PRs/collaborator clones | Explicit owner decision and coordination/migration plan | open |
| Contribution/security/CI surfaces missing | Repository baseline | External contributors lack policy and automated feedback | Issue-backed PR with validated files and workflow | in progress |
| `main` unprotected, nondesignated collaborators retain write, and security features are disabled | Private/free GitHub state | Direct pushes, undesignated merge authority, and unscanned dependencies/secrets remain possible | Post-public protection, collaborator-permission reduction, and security API verification | open |
| Contributor email metadata will become public | Existing Git history | Commit metadata contains personal/domain email addresses | Owner confirmation or explicitly approved author-rewrite strategy | open |

## Results
The owner selected Apache-2.0 and limited qualifying CODEOWNER review and merge authority to `vantasnerdan`, `axis-marbell`, and `mlops-kelvin`. Automatic deletion of exact same-repository heads after merge is enabled; a branch/PR reconciliation found no retained merged heads, while open and closed-unmerged/failed heads remain preserved. Discussions are enabled. Actions are restricted to GitHub-owned actions, immutable SHA pins are enforced, the default workflow token is read-only, and workflows cannot approve pull-request reviews. The `dependencies` label required by Dependabot now exists. Preflight established that the repository is private, `main` is unprotected, security scanning/update features are disabled, and no CI workflow exists on `main`. The all-history scan found no credential leak after adjudicating one identifier-list false positive. Rights review identified the exact publication-blocking source artifacts and confirmed the independent Zenodo tutorial's open GPL provenance.

## Canonicalization
This effort changes no scientific claim, release manifest, campaign, migration disposition, or generated scientific documentation. Durable coordination lives in issue #78, this effort record, the eventual contribution-readiness PR, and GitHub settings evidence. Any history transaction must separately migrate commit-pinned provenance before execution.

## Done Gate
The effort remains active. Contribution files, the owner history decision, a distinct merge, public visibility, branch protection/security/collaborator settings, and anonymous-clone verification are all still outstanding.

## Cross-References
- Canonical issue: https://github.com/vantasnerdan/substrate-framework/issues/78
- Coordination comment: https://github.com/vantasnerdan/substrate-framework/issues/78#issuecomment-5325114494
- Existing contribution workflow: `AGENTS_START_HERE.md`, `.github/pull_request_template.md`, `.github/ISSUE_TEMPLATE/agent_task.yml`
- Rights-sensitive sources: `incoming/einbein_1plus1D_tutorial.pdf`, `proposals/P229-preparata-qcd-vacuum-audit/sources/`
