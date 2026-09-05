# Exact vector-lift and reusable implementation receipt

Central registration precedes the derivation. schema.stdout/.exit report270
accepted claims and0. The first exact check ran in the separate Herdr pane:

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P251-cosserat-from-vortex-euler/attempts/0260/verify_vector_lift.py

verify.stdout/.exit record success and0. The original script is preserved as
verify_vector_lift_initial.py. A later Ruff check flagged the Python identifier
I as ambiguous (E741). The current script only renames that Python binding
to action; the symbolic coordinate, expressions, assertions and recorded
scientific result are unchanged. Targeted Ruff then passes. This was a style
repair, not a new scientific execution or a changed verdict.

The reusable coefficient inverse is extracted additively into
src/substrate_framework/euler_action_angle.py. It derives the shear from the
actual supplied angular velocity, exposes explicit resonances/nonconstant
curl inputs, and leaves physical chart/divergence assumptions explicit.
The earlier derivation script is frozen evidence of the algebra; new consumers
use the importable API.

The direct API tests independently construct a curl through exterior
differentiation with a variable volume Jacobian, then evaluate the actual
coordinate bracket and divergence of the returned lift. They include both
nonzero toroidal and mean-zero axisymmetric modes, and expose omission of
the angular shear. Resonance and nonconstant-curl cases are distinguished.

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python -m pytest tests/test_euler_action_angle.py -q

api-tests.stdout/.exit record the first execution:3 passed in2.36s, exit0.
No prior API, accepted claim or dependency is changed. Direct search finds
only the new tests as consumers of this new module.

route_verdict: established as stated for the exact two-direction
cohomological Kelvin lift, actual normal-displacement rows and the leading
nonzero distinct-polarization KKS pairing.

evidence_scope: exact analytic preparation, finite source bounds, symbolic
corroboration and tested coefficient inverse. No joint history rank, positive
Jacobi normalization, periodic action chart or compact geometry is inferred.

The failure-derived n=0 extension is proved from the actual zero-mean curl
and shear source, not by declaring the full cohomological operator invertible.
The direct axisymmetric API test corroborates its coefficient algebra.
The separate0250 executor computes the physical same-frequency gain and
full-pressure history implications.0259 accepted the narrower first
polarization; that review does not automatically review this new result.

Fixed repository validation passes exit0 after the direct API tests:
fixed-validation.stdout/.exit report270 accepted claims,1063 valid memory
files and43 unchanged memory warnings. Targeted Ruff and git diff --check
pass after the style-only repair. No unrelated source or test scope is
replayed for this additive module and proposal-evidence transaction.
