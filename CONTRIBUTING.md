# Contributing to Substrate Framework

Thank you for helping improve the framework. Contributions are welcome from
people and from transparently identified software agents. Scientific authority
comes from reviewed claims and pinned releases, not from confidence, chronology,
or a passing test alone.

## Before opening a pull request

1. Read [`README.md`](README.md), [`AGENTS_START_HERE.md`](AGENTS_START_HERE.md),
   and the normative [`AGENTS.md`](AGENTS.md).
2. Search existing issues, pull requests, claims, campaigns, and repository
   memory for overlap.
3. Create exactly one canonical issue before opening even a draft pull request.
   Use an issue form and state the positive objective, write boundary,
   dependencies, success gates, and coordination/merger plan.
4. Comment on the issue with the branch name and exact files or symbols you
   intend to change. Do not overwrite another contributor's active work.

Use `Advances #N` while any part of the issue objective remains. Use `Fixes #N`
only when the complete positive objective and all applicable governance gates
are satisfied.

## Development setup

The supported development baseline is Python 3.11 or newer. From a clone:

```bash
scripts/bootstrap.sh
scripts/validate.sh --full
```

Bootstrap creates `.venv` for the package and installs the repository's
`memory` command with `pipx`. For a bounded change, use the affected pytest
selectors justified by impact analysis:

```bash
scripts/validate.sh --pytest-scope tests/test_affected_module.py
git diff --check
```

Run the full validator for public exports, shared numerics or verification,
dependencies, conventions, governance semantics, claim promotion, releases,
or an uncertain dependency boundary. Do not duplicate an equivalent full run
at an unchanged boundary.

## Scientific and code contributions

- Put reusable APIs under `src/substrate_framework/` with focused tests.
- Keep exploratory composition under `proposals/` until adjudication.
- Never hand-edit `docs/generated/` or `migration/source-claims.yaml`.
- Declare every new public symbol as accepted-claim-backed, conditional and
  unpromoted, or non-scientific utility.
- Choose the oracle that matches the claim. Exact, formal, numerical, and
  simulation evidence are different verdicts.
- Demonstrate verifier sensitivity with meaningful mutations or
  counterexamples. A green tally that survives a broken input is not evidence.
- Preserve failed routes as append-only attempt provenance; do not present them
  as completion.

Documentation and workflow changes still require the canonical issue and
proportionate validation, but scientific oracle rows may be marked `N/A` with
an explanation.

## Pull requests and review

Use [the pull-request template](.github/pull_request_template.md). Keep artifact
merge, claim promotion, and goal completion as separate decisions. Describe the
smallest coherent units, exact validation commands, GitNexus impact, authority
status, remaining frontier, and any debt inside the proposed merge scope.

A person or agent must not merge a pull request that they opened, authored a
commit for, or materially implemented. A distinct reviewer with write access
must reproduce the load-bearing evidence and perform the merge. Direct pushes,
force pushes, and deletion of protected or unverified branches are prohibited.

The protected `main` gate requires approval from at least one designated
CODEOWNER: `@vantasnerdan`, `@axis-marbell`, or `@mlops-kelvin`. Public reviews
and comments are welcome but do not satisfy that code-owner gate, and public
contributors do not receive merge access. GitHub deletes the exact
same-repository head branch automatically after merge; open and
closed-unmerged/failed branches remain preserved unless the owner explicitly
retires them.

## Rights and sensitive material

By submitting a contribution, you represent that you have the right to provide
it under the repository's [Apache License 2.0](LICENSE) and any separately
stated documentation or data terms. Retain copyright and license notices for
third-party material. See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)
for the repository's current license exceptions and attributions.

Do not commit:

- credentials, tokens, private keys, personal data, or private correspondence;
- full paywalled papers or extracted full text without explicit redistribution
  permission;
- copied code, figures, datasets, or prose whose license is absent or
  incompatible; or
- generated host-specific files, local paths, caches, or virtual environments.

Prefer stable citations and public source links. If direct source retention is
scientifically necessary, record its checksum and access status outside Git
until redistribution rights are verified.

## Reporting problems

Use a bug or documentation issue for ordinary problems. Do not disclose a
security vulnerability in a public issue; follow [`SECURITY.md`](SECURITY.md).
Participation must follow [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).
