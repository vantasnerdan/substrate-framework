# Exact-check receipt

Binormal/NLS: first run2 failed/2 passed in8.04s. The failures were structural SymPy equality and unspecified reality of derivatives of real functions. Original test source and output are preserved; subtract-and-expand plus explicit real-jet substitution repairs the verifier only. Repaired4 passed in8.28s, exit0.

Finite-core bend:2 passed in2.14s, first-run exit0. No numerical quadrature, soft eigenvalue, fitted comparator or computed continuum stability verdict is used.

Commands: `PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python -m pytest -q tests/test_binormal_hasimoto.py` and the corresponding `tests/test_euler_bent_tube.py`. Runs occurred in Herdr scripts pane w3:p2, with stdout/stderr/exit captured on first execution.

The curve/PDE/metric/torsion and density/IBP/quantile identities are exact symbolic checks. The analytic global decay, Sobolev construction and relative-helicity integration-by-parts arguments are in finite-core-initial-bridge.md; the finite symbolic tests do not prove these analytic statements.

Final code boundary:

- `src/substrate_framework/binormal_hasimoto.py` SHA256 `0f612e7dc47f5edc89eb477e237fe6ae0af6633c4e33ef6ab45366681a20cdfc`

- `tests/test_binormal_hasimoto.py` SHA256 `e3bfeab7e049b4b317e5ce8feb1f51e986c989797b65accdbe449c0168373fa0`

- `src/substrate_framework/euler_bent_tube.py` SHA256 `cb489cef46a6181264a7c91b2863b4809f6ad68dcd56734436623430e78f1362`

- `tests/test_euler_bent_tube.py` SHA256 `3af1f82cf4a3757c30f50d41d2b13cc1515f9502c0d6633e628d20277eed2a8c`
