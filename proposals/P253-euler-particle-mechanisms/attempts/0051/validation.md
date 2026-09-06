# Validation and claim boundary

The load-bearing result is the analytic degree-of-freedom and applicability
audit in `construction.md`.  The importable API implements only the exact
Schwinger--Hopf algebra conditional on a supplied canonical doublet.

The focused repository-interpreter run was executed in the dedicated Herdr
scripts pane:

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python -m pytest -q tests/test_euler_schwinger_hopf.py tests/test_euler_quantum_two_state.py tests/test_euler_measurement_bridge.py

It reports `14 passed in 1.90s`, with empty stderr and exit `0`.  The tests
derive the Stokes Poisson algebra and Casimir symbolically, check invariance
under a nontrivial unitary mixing, expose the common Hopf phase, reject a
degenerate zero-action sphere and nonunitary mixing, and replay the directly
affected one-plane and detector APIs.  They do not establish an Euler mode
doublet, invariant projection, physical mixing operation, quantization,
measurement, or particle identity.

Pinned SHA-256 values are:

- README: `3d31ab6820254ba09c60ab1fde20fa7501bc732df5a5a890fd33768869d46d4d`;
- source audit: `dc3d53c8b90f5bc83079896e885d97767c62294a39096ad5983c1e1894b70fe3`;
- construction: `78802a7206506a7daa8ee012a9362cb746695adcebc9dac098c94976f995879c`;
- result: `4fdb8fe64ae135869fdc6f21553836a94bf1f04903537e2e681c69674949e65b`;
- API: `7bc7ead5a9a8eb800028a17a408e08a7a042d23d2ccb4fe2c74afda1932bd3d9`;
- tests: `02ac131ea1015581badd0835636c3f18755997bf786683b384219faace97f402`;
- command: `d8d75be024e6ecd2c6056ce7c1abcbd82e9fa2c7936b1c10172a25edde9a4f0c`;
- stdout: `5fc73e80996797494093e06cbd2b935387b453b382c21a146e148922c2c52ed6`;
- empty stderr: `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`;
- exit `0`: `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

Independent review is required before this exact conditional algebra and its
Euler nontransfer boundary are consumed as reviewed evidence.  P4 and the
parent campaign remain active.
