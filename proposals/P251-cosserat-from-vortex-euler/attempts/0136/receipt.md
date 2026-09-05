# 0136 frozen construction receipt

Owner: /root. Baseline a93e78e; accepted v0.175.0 unchanged.
`route_verdict: established` for the stated shaped column, local closed
torus, global planar Bernoulli lift and exact invariant planar action.
`evidence_scope: exact analytic construction with symbolic regression`.
The parent same-field continuum and global EPS geometry remain active.

The strongest oracle is the full analytic curl/Euler construction and
Dirichlet invertibility argument in shaped-forcefree.md, not the tally.
The new API is explicitly conditional, unpromoted infrastructure. Its
only direct consumers are its new unit tests and this attempt verifier;
no existing API, claim, release or generated surface changes.

First execution retained in first.stdout: 20/20, exit0. Final current
verifier adds the global lift and full KKS contraction: final.stdout,
24/24, exit0. It replaces one equivalent last derivative predicate;
there was no failed scientific predicate. The independent Cartesian
API tests first passed6/6; after adding explicit NaN/complex-infinity
level rejection, final-pytest.stdout passes6/6, exit0 (2.47 seconds).
Ruff on the module, test and verifier passes. No production numerics,
soft eigenvalue or numerical error threshold is used.

Commands (PYTHONPATH=src; Python is the repository environment):

    python proposals/P251-cosserat-from-vortex-euler/attempts/0136/verify.py
    python -m pytest -q tests/test_euler_forcefree.py
    ruff check src/substrate_framework/euler_forcefree.py tests/test_euler_forcefree.py proposals/P251-cosserat-from-vortex-euler/attempts/0136/verify.py

SHA256:

- shaped-forcefree.md: b7df051351847f38bfc547896a5989cc4bc61539eaa807dce1eccce5abe35a86
- verify.py: 8d57843e1aa864e6e91b20e55d3a1583dc7d1fef81370d7c6d904c6d93c058ab
- src/substrate_framework/euler_forcefree.py: 8d98928b33e7f6f3702c6a469c0edb56579ccffaf3380327ec7d08e1eb614702
- tests/test_euler_forcefree.py: 40f232baddde3bcfa3cf27b1fe07b2ad218de841a4058ef95768040ed0d99602

Next: use the exact global lift on0139's stationary array, preserve its
actual planar acoustic branch and construct core optical response on
that same field in0141. Constant-factor/global stationary geometry
routes run separately; none is supplied by a generalized Beltrami label.
