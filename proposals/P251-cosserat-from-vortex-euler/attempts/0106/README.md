# 0106 — importable exact compact-isovortical jet algebra

Owner `/root/orientation_construction`; write surface is this attempt,
`src/substrate_framework/euler_compact.py` and its matching unit test only.
Accepted base v0.171.0 is unchanged. Parent P251 remains active; this is
the bounded reusable-code obligation, not a scientific promotion.

Frozen target: extract the actual right-normal composition matrix for
`D xi=(div xi,div(xi cross omega))` and the separate generator/velocity
angular-adjoint rows from 0085/0103. Inputs are a declared derivative order
and actual exact vorticity jets, with no prototype, expected rank or desired
moment encoded in the API. Preserve all existing source verifiers.

The fixed proof route is exact differential-operator composition, checked
independently by differentiating polynomial fields and coefficient
functions at small orders. Angular rows are checked by independently
differentiating their physical test fields, including a displaced moment
origin. A new thin consumer may replay the existing modular fixture against
the importable builder, without repeating its expensive rational rank.
Exact arithmetic has no numerical soft-mode or floor claim.

Impact search found no existing right-normal jet API or consumers in src
or tests. No package-root export, accepted claim, convention, registry or
release changes are in this child boundary. The physics skill governance
and oracle references were read; a passing API test licenses coefficient
algebra, not rank, analytic field existence, stationary Euler backgrounds,
positive energy, physical inertia or parent completion by itself.

## Completed implementation and evidence

The new conditional module exposes `derivative_indices`, the immutable
`CompactJetSystem`, and `compact_isovortical_jet_system(jets, order=6,
radius=...)`. It accepts exact arbitrary-order actual derivatives, supports
an arbitrary moment origin, returns the full constraint and separate
generator/velocity angular matrices, and performs no rank selection.
Missing derivatives and floating inputs are explicit input errors.
No WKB helper was added: it was unnecessary for this coefficient API.

The first unit run passes 12 tests, including independent actual
differentiation of right-normal operators with variable coefficients at
orders 0/1/2; direct physical adjoint differentiation at two origins;
derivative-normalization/sign mutations; and domain/immutability checks.
The 35 adjacent Euler tests pass. Ruff and scoped diff checks pass.
The new thin canonical consumer reproduces every existing order-six
modular constraint and velocity-row coefficient, and the same six-row
modular rank 241, with 4/4 checks. No characteristic-zero rank was rerun.
The original source scripts remain immutable.

Commands, working directory `/tmp/pr199-completion`:

    /home/dan/substrate-framework/.venv/bin/python -m pytest tests/test_euler_compact.py -q
    /home/dan/substrate-framework/.venv/bin/python -m pytest tests/test_euler_affine.py tests/test_euler_orbit.py tests/test_euler_displacement.py tests/test_euler_fourier.py -q
    env PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P251-cosserat-from-vortex-euler/attempts/0106/canonical_consumer.py
    ruff check src/substrate_framework/euler_compact.py tests/test_euler_compact.py proposals/P251-cosserat-from-vortex-euler/attempts/0106

The first direct consumer command omitted `PYTHONPATH=src`; the shared
venv therefore looked at the main checkout rather than this worktree's
new module. Its native import error is preserved. The corrected command
selects the actual worktree source and passes without changing scientific
inputs or the environment installation. Pytest already selects this
worktree through the repository configuration.

`route_verdict: established`; `evidence_scope: SYMBOLIC_COEFFICIENT_ALGEBRA`.
The API and exact coefficient tests discharge this bounded importability
obligation. Parent claim review and completion remain the parent's active
transaction; this code extraction supplies no new scientific authority.
