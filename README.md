# Substrate Framework

A physics research framework for humans and software agents: discover useful
constructions, verify their consequences, and turn accepted results into
reproducible mathematics and importable code.

The goal is the advancement and full promotion of meaningful validated claims.
A scientific result earns acceptance through explicit assumptions, an exposing
proof or numerical check, individual review, and a pinned release. New ideas,
failed routes and useful intermediate constructions remain available for the
next attempt.

[Accepted claim index](docs/generated/claim-index.md) ·
[Latest release manifest](governance/releases/current.yaml) ·
[Agent starting guide](AGENTS_START_HERE.md) ·
[Contributing](CONTRIBUTING.md)

## What the framework has established

**Release snapshot: [v0.183.0](governance/releases/v0.183.0.yaml), September 5,
2026 — 271 accepted claims**, up from 207 in v0.160.0. This is a summary of
reviewed statements within their declared models and hypotheses. Acceptance
does not make an imported physical premise a derived result.

The registry's primary verification labels are **224 symbolic, 10 formal,
31 numerical and 6 simulation** claims. Individual claims can carry additional
evidence with its own scope. The [claim registry](governance/claims.yaml)
records statements, assumptions, dependencies, review status and evidence;
the generated index provides a readable entry point.

Recent additions include:

| Area | Accepted result | Scope and evidence |
| --- | --- | --- |
| Localized M5 clocks | A nonempty family of finite-energy, fixed-charge relative equilibria with intrinsic localization and orbital stability of the minimizing set. [C-M5C-003](docs/generated/claim-index.md#c-m5c-003), [C-M5C-004](docs/generated/claim-index.md#c-m5c-004). | Theorem-backed results for the declared tensor-plus-scalar action and auxiliary-frame quotient; particle identity and gravity are separate questions. |
| Clock walls and bags | Exact wall first integrals and thin-wall selection laws, with resolved wall tension and a seven-rung rotating spherical bag family. [C-M5W-001](docs/generated/claim-index.md#c-m5w-001), [C-M5W-005](docs/generated/claim-index.md#c-m5w-005), [C-M5W-006](docs/generated/claim-index.md#c-m5w-006), [C-M5W-007](docs/generated/claim-index.md#c-m5w-007). | Symbolic identities plus numerical evidence on the declared slice and thin-wall window. Fixed-frequency and fixed-charge stability are distinguished. |
| Optical variables and gravity | An invertible optical/ADM metric map and an equivalent optical-variable presentation of a supplied covariant action, including the accepted Einstein–scalar and worldline sectors. [C-GOT-002](docs/generated/claim-index.md#c-got-002), [C-GOT-006](docs/generated/claim-index.md#c-got-006). | Exact conditional field redefinition; it preserves the supplied equations without deriving microscopic gravity or Newton's constant. |
| Induced-coupling ledgers | Exact constant-mass one-loop coefficient families, regulator dependence, a total inverse-coupling composition, and conditions for its attractive sign. [C-IGR-001](docs/generated/claim-index.md#c-igr-001), [C-IGR-004](docs/generated/claim-index.md#c-igr-004), [C-GRV-002](docs/generated/claim-index.md#c-grv-002). | Declared operators, cutoff conventions and an explicit additive baseline remain inputs. |
| Confined-clock spectra and self-gravity | Kinetic-normalized fluctuation spectra, a zero-point mass correction, stress/constraint calculations and compactness obstructions. [C-M5S-010](docs/generated/claim-index.md#c-m5s-010), [C-M5S-011](docs/generated/claim-index.md#c-m5s-011), [C-M5S-012](docs/generated/claim-index.md#c-m5s-012), [C-M5S-013](docs/generated/claim-index.md#c-m5s-013). | Exact identities and numerical evidence for the specified confined backgrounds and coupling assumptions. |
| Euler to Cosserat/micropolar response | Prepared acoustic and rotational response, inherited action and physical-current matching, culminating in a compact-ring ensemble claim. [C-CST-015](docs/generated/claim-index.md#c-cst-015), [C-CST-016](docs/generated/claim-index.md#c-cst-016), [C-CST-017](docs/generated/claim-index.md#c-cst-017), [C-CST-018](docs/generated/claim-index.md#c-cst-018). | Prepared linear response and second variation through second spatial order on fixed finite time windows. The compact-ring claim has a follow-on audit question described below. |
| Lean formalization | Kernel-checked finite implications for declared polarization counts, gauge/chirality encodings, rotor bands and scalar identities. [C-GW-011](docs/generated/claim-index.md#c-gw-011), [C-EW-001](docs/generated/claim-index.md#c-ew-001), [C-ROT-002](docs/generated/claim-index.md#c-rot-002). | The formal theorem proves its encoded implication; physical interpretations and imported formulas keep their stated assumptions. |

The earlier foundation remains part of this release:

| Area | Examples |
| --- | --- |
| Sine-Gordon and localized profiles | Exact breathers, action and energy identities, stationary Q-ball profiles and fluctuation operators: [C-SG-001](docs/generated/claim-index.md#c-sg-001), [C-QBL-002](docs/generated/claim-index.md#c-qbl-002), [C-QBL-003](docs/generated/claim-index.md#c-qbl-003). |
| Gauge, topology and representation theory | U(1) identities, SU(3) algebra and representations, anomaly constraints, running-coupling and topological constructions: [C-U1-001](docs/generated/claim-index.md#c-u1-001), [C-IRR-001](docs/generated/claim-index.md#c-irr-001), [C-ANO-001](docs/generated/claim-index.md#c-ano-001). |
| Fields, moments and geometry | Optical/dilaton relations, stress moments, conditional radiation kinematics, vacuum polarization and long-range dipole interactions: [C-MOM-001](docs/generated/claim-index.md#c-mom-001), [C-SKY-001](docs/generated/claim-index.md#c-sky-001). |
| Units and classical worldlines | Dimensional ledgers, Lorentz orbits and stabilizers, massive/massless worldline mechanics: [C-DIM-001](docs/generated/claim-index.md#c-dim-001), [C-LOR-001](docs/generated/claim-index.md#c-lor-001), [C-WLN-001](docs/generated/claim-index.md#c-wln-001). |

Each link names a specific claim, not acceptance of an entire physical theory.
Use the [generated index](docs/generated/claim-index.md) for the complete set,
including oscillon simulations, chiral models, BPS/resolvent results, coherence
and kinetic models.

## Current research frontier

The current particle program is
[derive electron and neutrino mechanisms from the Euler substrate, electron
first (#203)](https://github.com/vantasnerdan/substrate-framework/issues/203).
It asks for localization, persistence, inertia, spin, charge or neutrality,
interactions, and a shared quantum/relativistic mechanism on the same physical
construction. These remain research objectives; a classical vortex or a prepared
continuum response does not by itself establish either particle.

The accepted Euler-to-Cosserat statements retain their finite-window,
prepared-response scope. A
[follow-on geometry audit](https://github.com/vantasnerdan/substrate-framework/blob/c6e950df1e5b97c40ea8e8c631a7725ee60c9980/proposals/P252-nonlinear-euler-objective-audit/attempts/0001/geometry-audit.md)
records a missing nonradial free-boundary inverse argument in the compact-ring
existence proof used by `C-CST-018`. The claim remains accepted in the pinned
registry; the audit is an open evidential challenge, not a counterexample to
Euler rings. [Issue #203](https://github.com/vantasnerdan/substrate-framework/issues/203)
records that dependency question and alternative particle routes.

Accepted canon is reviewable. A proposed replacement earns authority through
review and promotion; chronology, merge status and pass counts do not decide
the science.

## How discovery and validation work

Agents start with the physical question, sourced equations and constraints,
the unresolved tension, and useful mathematical structures. They distinguish
user-owned invariants from accepted premises and optional ansatz choices.
Structural analogies, representation changes, controlled limits and exploratory
calculations help generate candidates. A promising seed records its useful
consequence, missing bridge and next informative check.

The five repository skills give this work a shared practice:

| Skill | Purpose |
| --- | --- |
| [Physics discovery](.agents/skills/physics-discovery/SKILL.md) | Develop materially different mechanisms and constructions using focused context and informative exploratory checks. |
| [Physics Erdős loop](.agents/skills/physics-erdos-loop/SKILL.md) | Carry the scientific objective through derivation, verification, review and promotion. |
| [Small-ratio numerics](.agents/skills/small-ratio-numerics/SKILL.md) | Resolve soft spectra, weak forces and tiny splittings with measured error budgets and suitable representations. |
| [Theorem synthesis](.agents/skills/theorem-synthesis/SKILL.md) | Prove missing bridges between accepted claims and promote a useful higher theorem. |
| [Research PR harvest](.agents/skills/research-pr-harvest/SKILL.md) | Review existing eligible research PRs while separating reusable artifacts, claim acceptance and goal completion. |

Exploration generates hypotheses. Exact proofs establish their stated
implications; numerical evidence carries its resolution and error limits;
empirical tests address applicability to nature. Observations used to build a
candidate remain visible as inputs, with fresh consequences reserved for testing.

Useful intermediate results and failed routes are preserved on the campaign
branch. Scientific completion means the user's full objective has earned the
[ten achievements in AGENTS.md](AGENTS.md#the-ten-completion-achievements).
The terminal campaign PR records full success or evidence-backed scientific
exhaustion. Independently requested software and documentation improvements
have their own bounded objectives.

## Get started

The development baseline is Python 3.11 or newer, with `pipx` available.

```bash
scripts/bootstrap.sh
memory --version
.venv/bin/python -c "import substrate_framework"
```

Bootstrap creates `.venv`, installs the package and development dependencies,
installs the memory CLI through `pipx`, and sets up the pinned Lean/mathlib
environment. Use `scripts/check_lean.sh` when working on the formal surface.

For integrated repository validation:

```bash
scripts/validate.sh --full
```

For a bounded change, select checks from the actual diff and affected consumers:

```bash
PYTHONPATH=src .venv/bin/python scripts/validate_changed.py --base origin/main --head HEAD --print-only
scripts/validate.sh --pytest-scope tests/test_affected_module.py
git diff --check
```

The scoped example uses a placeholder test path. Follow the selector's result:
`--fixed-only` when no pytest scope is affected, the named selectors for bounded
changes, or `--full` for shared machinery, changed claims/contracts or uncertain
impact. The [operational guide](AGENTS_START_HERE.md) owns the complete workflow.

## Using the numerical APIs

`substrate_framework.numerics` supplies SciPy-backed initial-value,
boundary-value, method-of-lines PDE and refinement helpers. Claim-specific
code supplies equations, data, error norms and the scientific interpretation.

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

Canonical definitions live in [src/substrate_framework](src/substrate_framework);
tests and claim evidence document their accepted or conditional use.

## Repository map

| Path | Role |
| --- | --- |
| `src/substrate_framework/` | Importable equations, transformations, derivations and numerical tools |
| `governance/claims.yaml` | Claim statements, status, dependencies and evidence |
| `governance/releases/` | Pinned accepted claim sets; `current.yaml` selects the latest release |
| `proposals/` | Research candidates, active work and their evidence |
| `campaigns/` | Immutable adjudicated campaign records |
| `formal/SubstrateFramework/` | Lean definitions and theorem developments |
| `docs/generated/` | Documentation generated from the registry |
| `memory/` and `memory-templates/` | Sourced working state, recall and reusable research contracts |
| `.agents/skills/` | Discovery, execution, numerical, synthesis and review guidance |
| `tools/agent-memory/` | Source for the bundled memory CLI |

## Contributing

People and transparently identified agents are welcome. Begin with
[CONTRIBUTING.md](CONTRIBUTING.md) and
[AGENTS_START_HERE.md](AGENTS_START_HERE.md); [AGENTS.md](AGENTS.md) is the
normative scientific and governance contract.

Open one canonical issue before a PR, state the positive objective and ownership,
and work on a branch. Keep artifact merge, claim promotion and goal completion
explicit. Use a distinct merger by default; the user or owner can explicitly
authorize author self-merge of a named PR or bounded change. That operational
exception does not supply independent scientific review.

Preserve prior attempt evidence and adjudicated campaigns, and generate canonical
documentation from its source. Share only material you have the right to publish;
link to restricted papers rather than copying them into the repository.

- [Code of conduct](CODE_OF_CONDUCT.md)
- [Security policy](SECURITY.md) — private vulnerability reporting
- [Support](SUPPORT.md)
- [GitHub Discussions](https://github.com/vantasnerdan/substrate-framework/discussions)

## License

Apache License 2.0. See [LICENSE](LICENSE) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for attribution and exceptions.
