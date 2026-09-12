# Exact shear-memory validation

Four focused tests pass, exit0, in3.80s; command: PYTHONPATH=src with the project interpreter, pytest -q tests/test_euler_shear_memory.py. Complete output is memory-pytest.stdout/exit, produced in the dedicated scripts pane after the a7676b7 full-suite boundary.

The oracle differentiates a nontrivial finite-amplitude two-harmonic Euler field and checks all velocity components/divergence; integrates actual transported sine moments through sixth order to expose observed amplitude/initial-force parity and normalization; checks the full Laplace block identity with nonzero unresolved forcing and the short-time kernel limit; and verifies equal initial energies with opposite projected slopes. These are symbolic/exact checks, not numerical production evidence. The document proves the all-time bounded-operator elimination; the finite Taylor check alone is not that proof.

Canonical inventory inspected euler_time_control.py and euler_passive_control.py: they implement prepared finite-band/polynomial controls and wrapped-streamline observations, not the derived nonlinear triangular solution or exact complement kernel. The new additive module changes neither existing contract. No accepted claim or quantum/particle authority is attached to this API. Independent0011 review is active.
