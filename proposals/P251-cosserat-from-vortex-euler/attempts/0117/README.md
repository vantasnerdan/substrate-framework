# 0117 — actual material-tag spin and coarse translational response

Owner /root; parent issue198/P251 and the stronger same-core realization
remain active. Accepted conditional C-CST-008..010 are unchanged. Frozen
base release v0.175.0 at ef41e48. No empirical comparator or production
numerical remainder is involved.

Positive target: derive the actual moving-tag Fourier momentum map with
both intrinsic spin and symmetric shape-rate retained, then combine it
with the symmetric full Euler momentum stress and isotropic coherent
response. This constructs a first-gradient centroid/ambient translational
response to an actual changing tag spin, even when the complete Euler
point-filter dipole is constant. It does not assume canonical momentum
equals physical spin or infer an autonomous constitutive law from balance.

Candidate provenance:0098 already distinguishes hybrid and point momentum,
and0112 constructs actual same-tube material moments.0116 investigates
periodic background coherence independently. The failure-generated question
is whether complete packet cancellation really removes the physical
centroid response, or instead supplies the compensating orbital current.
The exact moving-domain identity is the selected proof route; arbitrary
tag truncation without ambient/shape/current terms is the rejected
representation, not a refutation of the physical mechanism.

Oracle: direct Taylor expansion of the defining Fourier integrals and
Reynolds transport, with an independent finite-mass polynomial fixture
for signs, the factor one-half, and the omitted shape-rate mutation.
SO(3) averaging is applied to the actual linear response tensors, not to
uniform marginals or an assumed oscillator. The analytic remainder is a
finite-radius Taylor bound, not a fitted error. Euler/Bloch regularity,
nonconstant observed spin and the shared optical window remain explicit
application inputs. Nonzero local spin alone does not establish them.

Implementation: a new unpromoted euler_observation.py with its own tests;
reuse conserved_moments.discrete_mass_moments for existing mass moments.
No existing API changes. Direct consumers are these tests and this attempt;
the other workers retain disjoint attempt directories.

## First evidence and route result

The first execution `PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python
-m pytest -q tests/test_euler_observation.py` returned zero: four passed
in3.30s, captured in `first-pytest.stdout`. Ruff initially found only the
ambiguous loop variable `l`; renaming it to `axis` changed no equation.
Ruff and diff checks then passed. The final scoped workflow will replay
the test at its edited source boundary.

SHA-256: module `21adf738d3ef2b4afbee4aa0b68931f37cffaa21e320536648d5841ec6d608ca`;
tests `e4f25050d881e4f4226a6f46dd61cc33e8384253372f5e930a0507a1a5d435bd`;
proof `0afe99b6cdb95ece60882f937b15dd79c8b1be6a9f654adcf65be2ffa9a932f2`.

`route_verdict: established` for the material Fourier identity and its
conditional isotropic Euler first-jet implication. `evidence_scope: exact
moment algebra and analytic finite-radius bound, with the actual Euler
response/ensemble inputs explicitly retained`. The result supplies a
physical centroid-spin coupling route, not the remaining autonomous
full-Euler constitutive identification.0114's nonconstant actual spin
and0116's periodic/Bloch construction are its immediate application inputs.
