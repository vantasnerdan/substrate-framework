# 0127 — reusable physical scalar chart

Owner /root. Extract the established exact0119 chart algebra for reuse
by0123/0124. This is unpromoted conditional infrastructure, not a new
Euler existence or physical-inertia-normalization claim. The original
and stronger campaign objectives and accepted v0.175.0 remain unchanged.

Frozen API boundary: append PhysicalScalarChart and physical_scalar_chart
to euler_phase.py; existing MovingPhasePullback and moving_phase_pullback
behavior is unchanged. Inputs are a constant nondegenerate2D symplectic
form, its compatible time-dependent generator and derivative, an actual
scalar observation and its first two derivatives, and measured spin row.
Output the full angle/rate transformation, kinetic and stiffness terms,
measured spin-rate inertia, spin connection, and normalized bracket.
The caller licenses actual derivatives and microscopic observation maps.

Oracle: exact differentiated coordinate equations and independent
Euler--Lagrange variation, physical/nonphysical frame counterexamples,
independent spin decomposition, and singular/incompatible domains. No
numeric spectrum, fit, or comparator. One strongest exact proof route
is sufficient for this fixed theorem extraction.

Impact: GitNexus959db88 reports no callers of moving_phase_pullback,
but rg supplies tests/test_euler_phase.py and attempts0114/0119/0126.
These known incoming edges are omitted by the graph and retained in the
impact boundary. No prior function body or accepted claim changes.
Replay the existing phase tests plus the new scalar-chart tests; select
the final fixed/scoped workflow after freezing this additive boundary.
The full2576-test receipt0122 remains valid for unchanged implementation.

## Implementation and first execution

The first targeted command `PYTHONPATH=src
/home/dan/substrate-framework/.venv/bin/python -m pytest -q
tests/test_euler_phase.py tests/test_euler_scalar_chart.py` passes all8
tests in4.11s, exit0; captured in first-pytest.stdout. Ruff passes.
The new tests derive coordinate dynamics by differentiation, check the
independent scalar action variation, and expose positive-but-mismatched
spin, retained connection, negative/degenerate observation winding and
incompatible/singular input domains.

Implementation SHA256454fea1db617e999e1aaeb39346670ad571e5da09241d58c57e31085c7639819;
new tests c3827abab233ae7fab62ae43733a380a018d0de402d477ed94e4995240673252.
GitNexus precommit detection sees only an indexed memory section because
the new symbols are not yet indexed; explicit source inspection fixes the
actual additive boundary. It is not treated as a complete call-graph result.

Final conservative workflow validation is complete below. Its selector
treats any existing framework-module edit as full, even this additive API;
the existing full0122 receipt is not misreported as validating new code.

Frozen candidate tree86bcf15ff674483c764fde968d7b33d03d3f4775 pins the
implementation/test boundary above. The first selector invocation supplied
that tree as --head and failed because its three-dot git comparison
requires a commit, not a tree (selector.stdout). This is an invocation
failure, not a scientific result. Source inspection of the selector's
changed_existing_modules rule selects full; full.stdout captures that
conservative execution. The commit-based selector will be checked after
the validated checkpoint, reusing the same full receipt.

Full command `PYTHON=/home/dan/substrate-framework/.venv/bin/python
scripts/validate.sh --full` completed with exit0:2580 tests pass in361.42s
and ALL REPOSITORY WORKFLOW CHECKS PASS.263 accepted claims,12 proposals,
1048 memory files and43 pre-existing warnings are unchanged. No canonical
source/test input changed during this run.0124's concurrently completed
proof and its20-check exact first receipt are attempt evidence only and
do not change that implementation boundary. Final record-sensitive
validation and diff checks accompany the checkpoint.

Route verdict: established for the typed conditional scalar-chart API.
Evidence scope: exact coordinate/action/spin algebra with independent
variation and exposed wrong-normalization/winding probes. This is reusable
infrastructure for the ongoing Euler construction, not parent completion.
