# Validation and claim boundary

The load-bearing result is the analytic degree-of-freedom and applicability
audit in `construction.md`.  The importable API implements only the exact
Schwinger--Hopf algebra conditional on a supplied canonical doublet.

The initial focused repository-interpreter run was executed in the dedicated
Herdr scripts pane:

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python -m pytest -q tests/test_euler_schwinger_hopf.py tests/test_euler_quantum_two_state.py tests/test_euler_measurement_bridge.py

It reports `14 passed in 1.90s`, with empty stderr and exit `0`. After the
physical-normalization audit added the positive-frequency energy-action
scaling to `canonical_doublet`, the same focused inventory was rerun in that
pane. A strengthened symbolic replay then put unequal positive frequencies
directly into the Poisson-bracket oracle and reports `14 passed in 2.00s`,
empty stderr, exit `0`; append-only final receipts are
`normalization-correction-v2.*`. The tests
derive the Stokes Poisson algebra and Casimir symbolically, check invariance
under a nontrivial unitary mixing, expose the common Hopf phase, reject a
degenerate zero-action sphere and nonunitary mixing, and replay the directly
affected one-plane and detector APIs.  They do not establish an Euler mode
doublet, invariant projection, physical mixing operation, quantization,
measurement, or particle identity.

The pre-0056 physical-coordinate checkpoint SHA-256 values were:

- README: `321f8429919a3b6c201693ff2220ef00114fcf44958dc910d0a00be6420c448b`;
- source audit: `e355dd7046e93595ff1f8062f391bd9f28380e11391ebe240c382408d11db50b`;
- construction: `8df69d4c88f675af7fae4446656856d1b5f4b5069c544bff64bc3bdf8024d25d`;
- result: `a889c0fb325bdb08db640dcc35dd5244c8a54be8d2db8fd471db48ce145dbe07`;
- API: `3e3a9981e8be6ca7f21cc6e30087dca2924c624e40c4ac019f08d7dc5b6f9e4e`;
- tests: `9868f214c652d8629b27ea7fe524e6fdd344e9d7022cd31bb4137c3fac257ee6`;
- corrected command: `d8d75be024e6ecd2c6056ce7c1abcbd82e9fa2c7936b1c10172a25edde9a4f0c`;
- corrected stdout: `ef3f62b56e136ea2c6fe5868700df626f91d87ca9ec74e69f125119686e4eae7`;
- empty stderr: `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`;
- exit `0`: `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

Independent review is required before this exact conditional algebra and its
Euler nontransfer boundary are consumed as reviewed evidence.  P4 and the
parent campaign remain active.

The subsequent 0056 review found and repaired the physical action moment-map
factor described in `action-moment-map-correction.md`.  The corrected oracle
now derives `H=sum_a nu_a |z_a|^2`, `{z_a,J}=-i z_a`, `|S|=J/2`, and reduced
area `2*pi*J`, including unequal positive frequencies.  The focused direct
consumer inventory was rerun and again reports `14 passed in 2.00s`, empty
stderr, and exit `0`; its append-only receipts are
`action-normalization-correction.*`.  All route verdicts remain unchanged.
