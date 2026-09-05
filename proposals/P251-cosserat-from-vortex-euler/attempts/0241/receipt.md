# 0241 receipt — exact joint residual and complete material momentum

Base integration a2fcb46; no accepted claim or old API changed. New additive
conditional API: src/substrate_framework/euler_joint.py. Tests are direct
finite material Fourier differentiation and actual angular/band integration.
The exact joint residual is joint-residual.md, equations(1)-(8).

The full moment expansion retains spin, shape rate, absolute-velocity
quadrupole and moving centroid phase. The branch residual contains the actual
history and observation errors; its explicit model mismatch begins at degree3.
The current superpotential retains q_dot and the independent initial charge.

Commands and first outcomes:
- Historical .venv pytest: exit1, no pytest installed; first-pytest.stdout.
- Main-workspace interpreter, initial tests: exit1, three pass and two
  Piecewise antiderivative tests fail; repaired-pytest.stdout.
- Corrected analytic-extension tests with the same interpreter: exit0,
  five pass; second-pytest.stdout.
- Initial joint verifier: exit1 after one pass, unsimplified structural
  equality at Eq6; joint-first.stdout.
- Second joint verifier: exit1 after three passes, immutable mutation copy;
  joint-repaired.stdout.
- Corrected joint verifier: exit0, nine checks pass; joint-third.stdout.

Exact working commands (from research root):

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python -m pytest tests/test_euler_joint.py -q
    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P251-cosserat-from-vortex-euler/attempts/0241/verify_joint.py

Initial source snapshots and diagnoses are preserved. The unused tautological
current predicate was replaced before it could supply evidence. Nine lexical
check calls and nine executed checks belong to the final joint verifier;
five runtime pytest cases include the two parameterized parity cases.

Route verdict: established exact interface. Evidence scope: material calculus
and exact symbolic algebra. An arbitrary finite material sample does not prove
Euler existence, and model-map substitution does not supply actual histories.
Next route executed in0243: both-parity same-cell passive controls for the
full hybrid acoustic acceleration correction. Parent campaign remains active.
