# Agent Contributors: Start Here

Read [AGENTS.md](AGENTS.md) as the normative reward and scientific contract.
This guide owns operational commands; skills own specialized methods. The elite
opening anchors the user's real objective and makes the next useful action clear.

## 1. Locate authority and ownership

Read `governance/releases/current.yaml`, relevant `governance/claims.yaml`
entries, source modules and proposal/issue. Search before allocating IDs.
Inspect the working tree and preserve unrelated edits:

```bash
git status --short --branch
git log -5 --oneline --decorate
gh issue list --state open --limit 30
gh pr list --state open --limit 30
```

Each PR uses one canonical issue created first. Comment with the positive
deliverable, branch, write surface, dependencies and ownership. Use
`research/<proposal>-<topic>`, `harvest/<issue>-<topic>`, or
`fix/<issue>-<topic>`. The user's existing authority governs actions;
honor a requested pause immediately.

## 2. Recall and prepare

```bash
scripts/bootstrap.sh
memory --version
scripts/check_lean.sh
memory search --base "$PWD/memory" "<issue, claim, mechanism>"
memory grep --base "$PWD/memory" -F "<claim-id>"
memory validate --base "$PWD" "$PWD/memory"
```

Use the existing environment when ready; bootstrap is for setup. Verify memory
facts against sources. A positional absolute memory target prevents host
configuration from redirecting validation. Start every filled section with a
prose description for memory indexing.

| Task | Skill / contract |
| --- | --- |
| Physics, claims, framework reconciliation | `.agents/skills/physics-erdos-loop/SKILL.md`; `campaign-proposal.md` or `research-arc.md` |
| Open mechanisms, new constructions, or a stalled route | `.agents/skills/physics-discovery/SKILL.md`; existing contract's discovery/candidate section |
| Higher theorem from accepted claims | `theorem-synthesis`; `memory-templates/theorem-synthesis.md` |
| Near-floor force, spectrum or splitting | `small-ratio-numerics`; `numerical-check.md` |
| Authorized bounded worker | `subagent-task.md` |
| Proposed/changed claim review | `claim-review.md` |
| Changed attachment roles | `evidence-attachment-review.md` |
| Existing external or terminal campaign PR | `research-pr-harvest`; `delegated-continuation-pr-template.md` |
| Other durable work | `effort-contract.md` |

Validate the matching scientific proposal before its frozen source/comparator
access: `PYTHONPATH=src .venv/bin/python scripts/validate_repository.py`.
Primary-source availability and previously exposed values are recorded honestly.

On a new or resumed research session, load `physics-discovery` before selecting
an open mechanism or replacing a route. Refresh the compact discovery brief in
the existing contract from the actual branch frontier; use its next informative
construction to begin work. A paused campaign resumes only on the user's
authority. Native skill discovery and direct file reading use the same local
SKILL.md; this routing also serves agents without native skill support.

## 3. Supervise structure and strategy

Give a worker one bounded question, exact sourced inputs, disjoint output and
a return condition. Keep approach selection with the supervisor. Use the user's
model/effort/permission settings and verify actual context when exposed.

Before consuming a worker result or extending its route, trace its equations
and implementation to actual consumers. The structural view is `ripwire`
(installed by `scripts/bootstrap.sh`). It reads the working tree directly, so
its answer is always current and byte-identical across machines, and it sees
the module-level calls that campaign attempts make into the framework. These
calls earn the view in under a second:

```bash
ripwire . --impact=src/substrate_framework/<file>.py:<symbol>  # blast radius, campaign callers included
ripwire . --situ=<file>        # what an edit to this file reaches, tests to run, co-changed files
ripwire . --test-gate          # the same for the whole working diff, before a commit
ripwire . --top-k=0 --expand=<file>:<symbol>  # one body with callee signatures, not the whole file
ripwire . --grep=<claim ID or mechanism> --grep-in=any  # hits with their enclosing YAML key or heading
ripwire . --recall="<claim ID or question>"   # the ranked pages themselves
```

Name a symbol as `file:name`; a bare name merges every same-named local copy
in attempts into one answer. Read counts as floors: call edges are name-based,
campaign verify scripts run as subprocesses and never show as tested, and a
zero means "none found". `governance/claims.yaml`,
`migration/dispositions.yaml`, and `migration/source-claims.yaml` sit above
ripwire's fixed YAML ceiling; read the registry with `rg` or through
`memory/framework/claims/`. When the tool is unavailable, inspect AST
definitions, calls, imports, branches and test assertions against current
source, and name that limitation in the result. Documentation-only changes can
have no symbol impact; their links, schema consumers, and process tests are the
relevant structure.

The supervisor separately checks mathematical logic, omitted perturbation
classes, alternative representations and the parent-objective payoff. Great
supervision redirects an unproductive lemma chain before another worker grows it.
A graph/AST cannot prove mathematics or completeness of the approach.

## 4. Implement and validate the useful object

Reusable definitions belong in `src/substrate_framework/` with exposing tests;
campaign adapters call them. Capture attempt evidence on first execution and
reuse unchanged receipts. Scientific details live in the physics references.

Use current-source/graph impact and the diff to select:

```bash
PYTHONPATH=src .venv/bin/python scripts/validate_changed.py --base <base> --head HEAD --print-only
scripts/validate.sh --pytest-scope tests/test_affected_module.py
git diff --check
```

Run validation and commit as separate invocations. Use `--fixed-only` when
no pytest scope is affected; `--full` for the cross-cutting/uncertain cases
in AGENTS.md. An additive public export with bounded known consumers can remain
scoped. Prose/template changes receive structural checks and affected process
tests. Record base/head or tree, selectors, command, result and limitations.
The reviewer reuses this receipt at an unchanged boundary; a duplicate CI or
post-merge suite is not an extra scientific achievement.

## 5. Deliver and review within the authorized boundary

Bank active scientific campaign work in commits until full success or certified
exhaustion. Independently requested process/software work uses its own issue
and PR. Use [.github/pull_request_template.md](.github/pull_request_template.md).
State separately artifact mergeability, claim promotion, campaign terminal state
and goal completion. `Advances #N` keeps an incomplete goal open; `Fixes #N`
records complete positive success.

Use one substantive review and one correction check. Preserve the strongest
supported result with the minimum evidence-led repair. Direct falsifiers,
missing/circular steps, incorrect dependency use and broken affected consumers
earn corrections; adjacent observations remain follow-ups. Existing accepted
dependencies retain authority at their exact scope.

A valuable source PR remains in `request changes`, `active refactor` or
`active harvest` while its repair can land. Apply AGENTS.md's terminal-close
test before unmerged closure. New scientific routes found during terminal
review return to campaign execution. Existing external PR harvesting follows
its own skill and never supplies an early campaign PR.

## 6. Merge and preserve continuity

Use a distinct merger unless the user/owner authorized self-merge of this named
PR or bounded change. Record the authorization and validation; it is operational,
not independent scientific acceptance. Merge the reviewed boundary, update the
canonical issue with landed/correct-next/history-only units and their evidence,
then sync memory and verify exact merged-head cleanup. Preserve unrelated and
closed-unmerged branches unless their owner directs retirement.

Promote claims through reviewed registry/release/generated-state changes.
For a pause, save branch, worktree, pending constructions, exact artifacts,
model/session recovery information and receipts. Resume on the user's authority.
A clean operational disposition can complete an external PR task while its
scientific parent stays open.
