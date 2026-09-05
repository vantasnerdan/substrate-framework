# 0115 — time-dependent Euler phase restriction with its physical gauge

Owner /root; stronger same-cell realization continuation after the reviewed
conditional008--010 promotion. This does not change their accepted statements.
The fresh0112 actual EPS-core Floquet construction supplies a new route;
0114 independently constructs its microscopic KKS/Krein and moment rows.

Frozen positive target: derive the full time-dependent symplectic action
pullback, moving projection/complement and physical-observable map. Expose
Floquet frame winding as a representation change rather than evidence of
positive physical energy. Exact finite-dimensional algebra is the oracle;
no numerical energy difference, eigenvalue tolerance or comparator is used.
The infinite-dimensional Euler application keeps its Hilbert/domain and
residual bounds explicit; a finite matrix identity does not prove them.

For the signed convention Omega(q,s)=B, L=B*s*qdot-H/2, the symmetric
one-form is -z^T Omega zdot/2. Freeze a constant ambient symplectic matrix
Omega, symmetric H, and a differentiable full-rank embedding E(t). The
derived quantities are Omega_E=E^T Omega E, Q=E^T Omega Edot,
H_eff=E^T H E+sym(Q), dotOmega_E=Q-Q^T. Variation must retain the latter
term. The actual projected generator is obtained from that same action,
not a copied Floquet frequency. A periodic gauge E->E F changes the
coordinate Hamiltonian and also the physical observable O E F.

Implementation surface: new unpromoted euler_phase.py plus its own tests;
no existing canonical definition or consumer changes. Direct imports and
GitNexus index show no prior symbol with these names; the only new direct
consumer will be its test file and this attempt. Tests independently vary
the original time-dependent action and include a winding-sign example,
the missing-connection mutation, and a full moving-complement identity.
One first-run stdout is retained. The source proof is the analytic ladder,
not a numerical invariant-subspace assertion.

## First evidence receipt

The first execution `PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python
-m pytest -q tests/test_euler_phase.py` returned zero: four tests passed in
2.97 seconds (`first-pytest.stdout`). Ruff and `git diff --check` passed.
The independent action variation retains the nonconstant symplectic-form
term; the winding example changes the coordinate energy sign while leaving
the physical trajectory unchanged. No numerical spectral claim is involved.

SHA-256: implementation `b4fcc265b171062569dd811801441a97e5477a578746b71c6eae72eaff503306`;
tests `333a69c7fe348755e39064282620a1c40c40e8e40da13f98271db52615720204`;
proof `b38c59670838386d3c627f0e40ddd7041b5f6bff3fb67c4f7051ea829690982b`.

`route_verdict: established` for the exact moving-action and complement
identities. `evidence_scope: exact finite-matrix identities with explicit PDE
domain obligations`. The next application is0114's actual EPS-core packet
energy and its physical observable map; this receipt does not establish an
autonomous Euler-invariant Cosserat continuum.
