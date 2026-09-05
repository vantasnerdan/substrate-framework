# Research continuation: actual smooth Euler to coupled Cosserat continuum

## Decision and status

On 2026-09-05 the repository owner explicitly approved separating the
unfinished research from the bounded review and repair of PR #199. This
issue (#200) owns that continuation; #198 remains the original scientific intake
and PR #199 retains #198 as its pre-existing canonical issue. The separation
is a user-directed scope decision, not a proof of success or exhaustion.
Research execution is paused by that decision. The scientific objective
remains open, with no reliable completion estimate.

**Bottom line:** there is a useful conditional continuum construction and
substantial actual-Euler dynamical progress. There is not yet one completed
proof joining the actual dynamics, physical observations, action, and
stationary compact-tube ensemble into the full requested continuum.

The durable research branch is
[`research/pr199-completion`](https://github.com/vantasnerdan/substrate-framework/tree/research/pr199-completion).
Its complete attempt archive, including unsuccessful routes and their
corrections, lives under `proposals/P251-cosserat-from-vortex-euler/attempts/`.
The last checkpoint before this handoff is
[`e88763a`](https://github.com/vantasnerdan/substrate-framework/commit/e88763a07f4e9a1ad67a9ae849bea237b0119c75).
The final preservation commit will be linked in this issue. Later research
is **not** being merged into main through the bounded PR #199 repair.

## Purpose: why this research is worth doing

The scientific prize is to explain how a fluid governed only by Euler can
support an effective medium with independently rotating microstructure,
internal torque, couple stress, and acoustic/optical branches. Ordinary
displacement elasticity cannot describe those independent rotations. A
successful construction would identify which vortex organization creates
the extra variables and which measured microscopic integrals determine
their coefficients, rather than inserting the desired elastic law.

That would connect topology and vortex dynamics to continuum mechanics in
a falsifiable way: geometry supplies persistent structures; dynamics must
supply the actual rotational response and its reaction on translation.
Topology alone is not a constitutive law. The work also explains when
decorrelation removes a coherent response without eliminating fluctuations.

Practical value now is a research toolkit and a family of conditional
models and exact benchmarks, not a validated model of a particular material
or experiment. Possible future uses include vortex-based effective-medium
models, tests of micropolar closures, and reduced models that retain angular
momentum correctly. No experimental prediction, universal material constant,
Standard-Model application, or engineering design claim has been earned.

## What can already be used, and how

| Reusable result | Use | Boundary to retain |
|---|---|---|
| Repaired Rankine mode API and direct Cartesian/Bessel checks | Derive isolated-vortex mode equations and precision-controlled branch residuals | Conditional, unpromoted PR #199 infrastructure; not an elastic-modulus derivation |
| Micropolar energy and Fourier operator | Differentiate stresses, compute coupled dispersion, check conventions and singular limits | Caller supplies justified coefficients; importing the API does not prove their Euler origin |
| Compact Euler KKS/action constructions | Compute genuine phase, energy, core-angle and spin integrals on declared cells | Full physical observation differs from an arbitrary canonical coordinate |
| Conditional Euler-material-action continuum | Study a declared phase Cauchy–Born closure with computed coefficients and reaction terms | Restricted variational closure, not proof that unrestricted Euler trajectories stay in that family |
| Actual acoustic and optical preparations | Benchmark finite-window Euler response against reduced equations; test pressure, current and preparation corrections | Respect each theorem's cell, observation, time scale, amplitude order and spatial-error limits |
| Finite-core geometry and local persistence calculations | Supply/test actual invariant tubes and return-map twist on the stated field | A periodic noncontractible tube is not a compact Euclidean knot or a stationary array of rings |
| Full material torque and boundary-current identities | Keep ambient reaction, initial angular charge and boundary flux when changing continuum representatives | Equality of conservation laws alone does not identify a constitutive law |

Start with the module docstrings and direct tests, not the attempt tally.
On the research branch, relevant modules are `micropolar.py`,
`euler_fourier.py`, `euler_phase.py`, `euler_orbit.py`, `euler_compact.py`,
`euler_core_packet.py`, `euler_acoustic.py`, `euler_optical_response.py`,
`euler_displacement_preparation.py`, and `euler_passive_control.py` under
`src/substrate_framework/`. Their corresponding `tests/test_*.py` files
are executable usage examples. `rankine_modes.py` is also in PR #199.

For an isolated checkout with project dependencies installed:

```sh
git fetch origin research/pr199-completion
git worktree add --detach ../p251-research origin/research/pr199-completion
cd ../p251-research
PYTHONPATH=src python -m pytest tests/test_micropolar.py tests/test_euler_passive_control.py -q
```

Run only the tests/verifier for the result being reused. Read its exact
claim and exclusions in `governance/claims.yaml`, then its pinned source
and review. Full validation is not a prerequisite to inspecting this work.
Research-only modules are not promised to exist in main after #199 merges.

## Authority: what is accepted, and where

The research branch records nine individually reviewed claims,
C-CST-008 through C-CST-016, in branch-local releases v0.175.0–v0.181.0.
The latest branch release contains 269 accepted entries overall. These
records are preserved research-branch governance state, **not a claim that
those releases or APIs have landed in main through #199**. The proposed
original C-CST-001..007 ladder and its historical completion assertions
must not be substituted for the later scoped statements.

| Branch claim | Positive result and essential scope |
|---|---|
| C-CST-008 | Actual compact core/reaction pairs, positive canonical forms and measured moments on the declared stationary Beltrami tube law |
| C-CST-009 | Exact conditional second-gradient micropolar action under the explicit phase Cauchy–Born closure; unrestricted dynamical invariance remains separate |
| C-CST-010 | Product-Haar coherent torque/stiffness cancellation, with fluctuations and the separate dynamical closure retained |
| C-CST-011 | Actual finite-window optical Kelvin packet in the same persistent tube, with physical tagged observations and controlled errors |
| C-CST-012 | Actual planar Euler acoustic-window theorem on its specific array and observation |
| C-CST-013 | Full physical acoustic second-jet/current identities and controlled same-background optical insertion; moving action retained |
| C-CST-014 | Smooth stationary column optical-response construction, within its registered preparation and observation scope |
| C-CST-015 | Actual fixed-cell positive displacement preparation and initial physical phase |
| C-CST-016 | Same-cell common displacement/velocity prepared acoustic second-order limit on fixed time windows |

Exact statements, dependencies, reviews and source paths are in the
[checkpoint registry](https://github.com/vantasnerdan/substrate-framework/blob/e88763a07f4e9a1ad67a9ae849bea237b0119c75/governance/claims.yaml).
The later results below are research evidence, not new promoted claims.

## Progress map for the entire continuation

The complete archive is retained; this map gives entry points so a successor
does not have to read hundreds of attempts sequentially. Earlier assertions
can be superseded by later append-only corrections. Use the registry and
the final receipt for a result, not chronology or a green tally, as authority.

| Attempt range / anchor | Contribution and present use |
|---|---|
| 0001–0025 | Original seven-node proposal, frame identities, mode/dispersion work, review defects and repaired verifiers |
| 0026–0034 | Corrected Euler/Coriolis and Bessel correspondence; closed six-point-vortex angle action; prescribed-flow elliptic patch; conditional long-wave map; corrected mutual-energy kernel. This is the existing PR boundary |
| 0035–0104 | Finite-core, physical-current and compact-reaction exploration leading to the actual compact-pair and conditional action construction; use 0085, 0095, 0097, 0098, 0102, 0103 as positive-result anchors |
| 0105–0113 | Frozen transaction, independent claim reviews (0108), population limits, and branch promotion of C-CST-008..010 |
| 0114–0159 | Actual Euler/Lin optical histories, finite packets and geometry/pressure transfer; 0142/0145/0147 are construction anchors, 0157/0159 the reviewed C-CST-011 transaction |
| 0160–0178 | Planar acoustic-window and full observed-cell/current continuation; C-CST-012 review at0165, implementation at0167; subsequent claims are indexed by their registry evidence |
| 0179–0208 | Smooth optical suppliers, positive fixed-cell displacement, common-velocity action and actual finite-control limits; key reviews0197/0207 and APIs0199/0208 underpin C-CST-015/016 |
| 0209–0219 | Fixed positive-tag clocks, actual field-changing Kelvin lift, and same-ring geometry; 0211 constructs a literal constant-curl closed core, independently reviewed0215; 0216/0218 separate a label clock from an actual velocity-changing mode |
| 0220–0230 | Full curved pressure feedback and a positive bending pole; actual Kelvin action normalizers; conditional joint branch-to-Cosserat algebra and mixed second jets. Reviews0224/0225/0229/0230 cover their stated boundaries |
| 0231–0232 | Actual acoustic angle/axial-spin repair and exact material pressure-torque, ambient reaction and boundary-equivalence identities |
| 0233 | Curved Kelvin first-carrier and geometric/current calculations; partial construction with its final pause receipt, not a completed common-field supplier |
| 0234 | Two positive tag fractions on one Euler field give actual optical angle/full-spin second-row controls and a stated ordered error construction; 34 exact checks, independent review not executed |
| 0235 | Full acoustic spin trace and exact boundary-action improvement preserving physical momentum, including initial charge and time-dependent memory; 39 exact checks, no independent review of this new construction yet |
| 0236 | Same-C016 periodic finite core, nonzero flux twist and fixed positive density (15 checks); periodic-image forcing expansion (8 corrected checks), not a steady ring-array existence theorem |
| 0237 / 0240 | Preregistered independent review contracts for0236/0234; source review was not executed before the owner-directed pause |

Some original verifier executions failed and were repaired. Both outputs and
diagnoses are preserved. Counts above describe individual exact algebra or
implementation checks; they neither certify every analytic estimate nor add
up to a proof of the parent theorem. No new full-suite run was performed
for this handoff.

## Decisive remaining scientific construction

Construct one common stationary smooth Euler ensemble and preparation whose
**actual measured histories, inherited action, full momentum/angular current,
and finite-core geometry** have the same controlled coupled Cosserat limit.
That is the missing bridge, independently of which route eventually wins.

The current route has three concrete open pieces:

1. **Joint physical/action/current compatibility.** Combine the actual
   acoustic and optical suppliers with the full mixed action normalizer;
   evaluate the reduced dynamical residual and compare the literal total
   current divergence. Retain initial G, the q-dot boundary memory, ambient
   reaction, directional pressure projection, and one compatible error
   ordering. Action normalization alone does not supply actual dynamics.
2. **Compact Euclidean geometry at usable stationary density.** The same
   periodic cell already has a finite invariant tube. The distinct compact
   ring route still needs a genuine stationary ensemble/density realization
   compatible with its dynamical scaling. Periodizing/superposing isolated
   rings creates an Euler residual. The positive-frequency inverse in0220
   is not the missing zero-frequency steady inverse.
3. **Final theorem and reuse transaction.** After the bridge is proved,
   individually review the new claims, extract reusable definitions/tests,
   replay affected consumers and reconcile the release/registry/docs. No
   full-parent promotion is justified merely by the existing local results.

The straight-tube formula alpha=L_v T/6 remains unsupported: exact free
director rotation cancels in Green–Lagrange strain. The positive conditional
joint-symbol relation alpha=j nu^2/4 is a different statement; its actual
microscopic suppliers must be joined, not assumed. Likewise a bending pole,
an optical packet, or an EPS existence result alone is not this continuum.

## Restart contract and coordination

This issue is a discoverable research backlog, not a request to resume an
unbounded run automatically. Obtain a new owner-approved effort boundary
before restarting. The useful first bounded deliverable is a joint supplier
statement with one dependency/error/current table and an evaluated residual:
success makes the same-cell periodic construction coherent; a mismatch names
the concrete source row to repair. Existing independent review contracts
0237/0240 can be resumed without re-reviewing accepted inputs.

Report at the agreed reassessment point what new mathematical object was
constructed, what remains, and whether the proposed next route materially
changes the prospects. Do not use another checkpoint or check count as an
implicit promise that full completion is near. No scientific-exhaustion
certificate exists; plausible routes remain.

Coordination: continuation code, attempts and branch-local governance stay
on `research/pr199-completion`; keep historical attempts immutable. Preserve
the original issue #198 and PR #199 as provenance. Future publication/harvest
requires an explicit bounded scope and appropriate review; this issue is not
authority to merge the entire research branch or silently narrow the theorem.
