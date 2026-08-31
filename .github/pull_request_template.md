## Objective and issue

Describe the positive deliverable and the smallest coherent merge boundary.

Canonical issue: <!-- Required. The issue must exist before even a draft PR is submitted. -->

Issue relationship: <!-- Use `Advances #N` while work remains; use `Fixes #N` only when the full goal is complete. -->

PR eligibility: <!-- Non-campaign change; externally supplied existing PR; terminal campaign success; or terminal campaign exhaustion. An active campaign rung/milestone is ineligible. -->

Campaign terminal evidence: <!-- N/A, complete positive-success receipt, or scientific-exhaustion certificate. -->

Active campaign relationship: <!-- none, or terminal campaign id; work that enables an active campaign is campaign work regardless of file type. -->

Authoring agent:

Intended merger: <!-- Distinct by default; otherwise link explicit owner/user self-merge direction. -->

- [ ] The canonical issue existed before this PR was submitted.
- [ ] The merger is distinct, or explicit owner/user self-merge direction is recorded.

## Change classification

- [ ] Reusable implementation or tooling
- [ ] Documentation or workflow
- [ ] Scientific proposal or campaign evidence
- [ ] Claim promotion transaction
- [ ] Harvest from an externally supplied incomplete research PR
- [ ] Terminal scientific campaign: positive completion
- [ ] Terminal scientific campaign: certified exhaustion
- [ ] Compatibility-only repair

## Authority and scope

Base release and commit:

Accepted claims and canonical sources used:

Declared imports, assumptions, and conventions:

Files or concerns intentionally out of scope:

New or materially changed public interfaces:

| Symbol or module | Authority status | Accepted claim IDs | Owning issue/PR | Conditional boundary |
| --- | --- | --- | --- | --- |
|  | accepted claim / conditional unpromoted / non-scientific utility | none or exact IDs |  |  |

An exported symbol is not accepted scientific authority by itself. For every
conditional unpromoted API, state its explicit inputs, assumptions, exclusions,
and why it is still independently reusable.

## Useful-unit disposition

List only independently mergeable units or materially different risk
boundaries. Do not atomize every file, theorem entrypoint, or evidence row. A
mergeable unit must remain useful and true without an unearned headline claim.

| Unit | Local claim or purpose | Dependencies | Evidence and tests | Proposed disposition |
| --- | --- | --- | --- | --- |
|  |  |  |  | merge / minimum correction / history only |

Answer each decision independently:

- Artifact merge: <!-- yes/no and why -->
- Claim promotion: <!-- none, or exact claim IDs -->
- Campaign terminal state: <!-- N/A / positive completion / certified exhaustion -->
- Goal completion: <!-- yes/no and the remaining gate -->

## Validation receipt

Record the base/head or tree hash, frozen impact surface, exact commands, and
results once. Include only claim-appropriate sensitivity and affected consumers.
Reuse this receipt while its inputs are unchanged; review prose and counts do
not stale it.

```text
# command
# exit status and meaningful verdict
```

## GitNexus impact

Record index freshness, pre-change symbol/API impact when applicable, and final
diff change detection. For documentation-only work, state that no graph-mapped
symbol or process was expected or found.

## Memory, governance, and generated state

Durable contract or decision entry:

Proposal/campaign/claim/release changes:

Generated outputs and synchronization commands:

Debt remaining inside the proposed merge unit: <!-- must be empty for acceptance -->

Campaign frontier outside this merge unit:

## Validation boundary

Mark a nonapplicable scientific row `N/A` and explain why in the verification
section; do not claim it passed. An additive public export may stay scoped when
impact analysis shows a bounded sector, no changed existing contract, known
consumers, and targeted API coverage. Record that rationale rather than treating
the export alone as a full-suite trigger.

- [ ] Targeted tests and named scientific verifiers pass.
- [ ] Load-bearing mutations, counterexamples, or wrong-convention probes fail as expected.
- [ ] Affected downstream consumers replay.
- [ ] `scripts/validate.sh --pytest-scope ...` passes with the exact selectors recorded, `--fixed-only` is justified by no affected pytest scope, or the full-suite trigger is explained and `scripts/validate.sh --full` passes; the recorded local `scripts/validate_changed.py` decision agrees.
- [ ] The pytest scope remains valid against the merge base; an equivalent unchanged validation is not duplicated.
- [ ] `git diff --check` passes in a separate invocation.
- [ ] No unrelated, generated-by-hand, or host-specific artifacts are included.

## Author handoff

State what the reviewer should independently reproduce, the riskiest assumption,
and the next decisive action if this PR advances rather than completes the goal.

---

## Reviewer disposition

Review the actual head commit and complete this section or provide the same
fields in a formal review.

- Canonical issue predates PR: <!-- yes/no -->
- Authoring or implementing agent:
- Merger and basis: <!-- distinct, or explicit owner-authorized self-merge -->

- Artifact merge: <!-- yes/no with unit-level rationale -->
- Claim promotion: <!-- none or exact claim IDs and review evidence -->
- Campaign terminal state: <!-- N/A / positive completion / certified exhaustion, with evidence -->
- Goal completion: <!-- yes/no and still-open gate -->
- Merge as written: <!-- yes/no -->
- Correct or harvest: <!-- exact unit, minimum repair, and one landing check -->
- Leave in PR history: <!-- exact units and rationale -->
- Source PR lifecycle: <!-- request changes, active refactor, active harvest, merged, or terminal closed -->
- Terminal-close evidence: <!-- N/A, or qualifying reason and landed replacement links -->
- Next decisive action: <!-- one concrete step -->
- Head branch disposition: <!-- auto-delete after merge; otherwise preserve unless owner retirement follows a recorded terminal-close decision -->

### Reviewer checks

- [ ] Base release, issue, accepted boundary, diff, memory, and discussion inspected.
- [ ] The canonical issue existed before the PR was submitted.
- [ ] The merger is distinct or explicit owner/user self-merge direction is linked.
- [ ] The load-bearing result and existing validation receipt were checked once at the frozen boundary.
- [ ] Dependencies, conventions, imports, consumers, and GitNexus impact audited.
- [ ] Every new public symbol has an explicit authority status and owning issue.
- [ ] Applicable verification sensitivity and numerical/formal limits were audited without duplicate replay.
- [ ] Merge, claim-promotion, and goal-completion decisions kept independent.
- [ ] A scientific campaign PR is terminal for the whole campaign, not a rung, milestone, partial goal, utility, or clean checkpoint.
- [ ] Any work that enables an active campaign stayed on its campaign branch until terminal eligibility opened.
- [ ] Review found no active obligation, untried plausible route, broken license edge, or invalid coverage argument that revokes terminal eligibility.
- [ ] A finite correction has the minimum repair and one landing check; ownership is named only when handed off.
- [ ] Any proposed unmerged closure passes the terminal-close test and records landed replacements.
- [ ] Canonical issue handoff is posted or preserved ready to post.
- [ ] The exact same-repository head will be deleted after merge; open or closed-unmerged/failed heads are preserved unless an owner explicitly retires them.
