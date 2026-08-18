# Substrate Framework

This repository is the cohesive, reviewable successor workspace for the Substrate physics framework. It starts with workflow and governance infrastructure; scientific claims are migrated only after claim-by-claim reconciliation.

The central rule is that chronology is not authority. Campaigns are immutable research records, proposals are unaccepted work, and the canonical framework is a reproducible materialized view of individually accepted claims.

## Contributing

Start with [`CONTRIBUTING.md`](CONTRIBUTING.md), whether you are contributing
directly or through an agent. It explains the issue-first workflow, authority
boundaries, validation expectations, and pull-request review rules. Community
participation is governed by [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md), and
security vulnerabilities should be reported privately as described in
[`SECURITY.md`](SECURITY.md). General usage questions belong in
[GitHub Discussions](https://github.com/vantasnerdan/substrate-framework/discussions).

Do not upload paywalled papers, private correspondence, credentials, or other
third-party material without verified redistribution rights. Cite or link to
the authoritative source instead.

## License

Unless a file or directory says otherwise, Substrate Framework is licensed
under the [Apache License 2.0](LICENSE). Separately licensed material and
attribution are listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

### Contributing agents

Start with [`AGENTS_START_HERE.md`](AGENTS_START_HERE.md). It is the operational
guide for choosing and coordinating work, using repository memory and
GitNexus, loading the native skills, opening a reviewable PR, and performing an
independent review. The root [`AGENTS.md`](AGENTS.md) remains the normative
scientific and governance contract.

Pull requests use [the repository template](.github/pull_request_template.md)
and must keep artifact merge, scientific claim promotion, and goal completion
as three separate decisions. Every PR requires a canonical issue created before
submission, and an agent may not merge a PR it opened, committed to, or
materially implemented.

## Repository model

- `src/substrate_framework/` contains importable framework definitions and derivations.
- `governance/claims.yaml` is the machine-readable accepted/proposed claim graph.
- `proposals/` contains candidate campaigns before adjudication.
- `campaigns/` contains immutable adjudicated campaign records.
- `governance/releases/` pins reproducible accepted claim sets.
- `docs/generated/` is generated from the accepted registry; agents do not hand-edit it.
- `memory-templates/` contains durable work, research, review, and promotion contracts.
- `tools/agent-memory/` contains the memory CLI program only. No prior memory entries were copied.
- `.agents/skills/physics-erdos-loop/` contains the repository-scoped native Codex physics workflow.

## Bootstrap

```bash
scripts/bootstrap.sh
scripts/validate.sh --full
```

Bootstrap creates `.venv` for the importable physics package and its NumPy,
SciPy, SymPy, and test dependencies. It installs the bundled `agent-memory`
release with `pipx`, not project `pip`, so `memory` remains available without
activating `.venv`:

```bash
memory --version
.venv/bin/python -c "import substrate_framework"
```

For a bounded pull request, keep the fixed repository checks and restrict only
the pytest stage to the affected tests identified by change-impact analysis:

```bash
scripts/validate.sh --pytest-scope tests/test_affected_module.py
```

Use `--full` for promotion, release, periodic integrated-main replay,
cross-cutting changes, or an uncertain dependency boundary. Invoking
`scripts/validate.sh` without arguments remains an alias for the full suite.

The bundled CLI is the code-only `agent-memory` v0.2.0 release; no user or
agent memory entries are included. See
[`tools/agent-memory/UPSTREAM.md`](tools/agent-memory/UPSTREAM.md) for the pinned
upstream source.

## Numerical physics APIs

`substrate_framework.numerics` provides shared SciPy-backed IVP,
method-of-lines PDE, BVP, and refinement-evidence helpers. They standardize
failure handling and evidence capture; they do not turn a numerical result
into an exact proof. Claim verifiers must still supply the governing equation,
boundary/initial data, convergence study, invariants, and independent checks
appropriate to the claim.

```python
import numpy as np

from substrate_framework import SolverTolerances, solve_ivp_evidence

orbit = solve_ivp_evidence(
    lambda _t, state: np.array([state[1], -state[0]]),
    (0.0, 2.0 * np.pi),
    [1.0, 0.0],
    tolerances=SolverTolerances(rtol=1e-10, atol=1e-12),
)
```

Read [`AGENTS_START_HERE.md`](AGENTS_START_HERE.md) and `AGENTS.md` before
starting any research or migration. A fresh effort begins by instantiating the
appropriate file from `memory-templates/`; it does not begin by editing
canonical prose.

## Authority state

The framework began from an intentionally empty registry. Its current accepted
boundary is always the pinned manifest in `governance/releases/current.yaml`
and the individually accepted entries in `governance/claims.yaml`, not this
README or the newest commit. No predecessor claim is accepted merely because it
was late, committed, numerically attractive, or described as settled.
