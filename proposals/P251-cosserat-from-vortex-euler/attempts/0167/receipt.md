# 0167 direct implementation receipt

The reusable triangular Euler field and full-pressure axial Kelvin
initial phase are implemented in euler_acoustic.py. Five direct tests
passed on first execution in3.98s, exit0, captured in first-pytest.stdout.
Ruff passes. No full test suite or numerical spectral computation ran.

Tests independently differentiate the Cartesian field and pressure,
integrate the actual complete oblique-cell covariance, check the
nonzero separatrix dual, derive the finite-k pressure factor, expose
the unprojected mutation, and test a different two-mode spectrum and
input domains. The API supplies no prescribed acoustic trajectory.
General-spectrum commensurability/stationarity remain explicit caller
hypotheses;0161 supplies the triangular instance's analytic license.

Scientific hashes:

- euler_acoustic.py: 3e71921ea1e4986e9fa9d07870a49e589b70aa37313a80a3801ffbe9eece792e
- test_euler_acoustic.py: efb2748bcd74ecd73ef7dfacd714fb5c1d104d98dc12fc5641a7392937553b5d
- first-pytest.stdout: 72d29ea82839de4ed76e4d684bbd689ec686e46ab47cf2711a4ee6fcd114acdb

route_verdict: established

evidence_scope: exact field and initial-action extraction with independent
integral/differential/projector tests. Actual acoustic-time dynamics
comes from0161's reviewed analytic proof, not these five tests alone.
0168 materializes the claim. Parent coupled continuum remains active.
