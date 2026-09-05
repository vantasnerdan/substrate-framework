# 0100 — reusable exact Euler Fourier algebra

Owner `/root/construction_review`; bounded implementation task under P251.
Frozen target: extract general finite-Fourier operations from immutable
0040/fourier_orbit.py into `substrate_framework.euler_fourier`, with exact
density and nonzero Beltrami-eigenvalue handling. Keep coadjoint energy
Hessian H and material Jacobi stiffness K as DIFFERENT explicit APIs.
Compute K directly from physical pressure Hessian and background transport,
and test H-K=rho||A xi||² with the independent material Kelvin operator.

Authority: base release v0.171.0. This is a conditional unpromoted utility;
no accepted claim, release, old attempt, or package-root export changes.
Read physics skill plus governance/oracle references completely before
editing. Impact search found no existing Euler finite-Fourier or Leray API
in src/tests; old0040 consumers remain immutable. Existing unrelated sampled
Fourier/radial and micropolar matrix APIs are different mathematical objects.

Fixed proof/implementation route, no empirical comparator or numerical
approximation: coefficient convolution, exact Fourier derivatives and
Leray projection, rational smooth fixture, independent physical-space
pressure/transport identity, density/lambda/sign and zero-mode mutations.
Wavevectors are rational and commensurate; outputs are common-cell volume
AVERAGES, not unnormalized integrals. Default Leray preserves the mean;
an explicitly zero-mean projector is a different available convention.
Targeted new tests and lint are the complete affected consumer replay.

## Implementation and evidence receipt

The new module supplies exact scalar convolution, derivatives, vector
products, transport, volume-average pairing and explicit Leray mean-mode
selection. `coadjoint_matrices` checks the background's actual nonzero
Beltrami eigenvalue and generator incompressibility, and returns the full
Eulerian velocity/H/KKS matrices with density. `material_jacobi_matrix`
instead consumes PHYSICAL pressure and independently computes its Hessian
and the convective quadratic form. `material_kelvin_operator` directly
forms P[(u.grad)xi+(Dxi)^T u].

With F=xi cross u and v=lambda P F, curl F=(u.grad)xi-(xi.grad)u and
A xi=curl F-v. Thus the exact identity H-K=rho||A xi||² follows by
expansion, with the periodic mean retained. The tests evaluate its matrix
polarization as well as the quadratic statement. Positive H is not silently
identified with positive K.

The independent pressure test differentiates rational real trigonometric
fields in physical space, checks stationary Euler directly, and performs
exact period integrals; it does not reuse Fourier derivatives or zero-mode
extraction for that side of the oracle. It also verifies H minus this
independently computed K equals the direct Kelvin Gram. Nonzero positive
and negative eigenvalues, rational density, nonzero KKS, time reversal,
Leray constants, and missing-density/helicity/mean-mode mutations are tested.

Initial run:8 passed,1 failed (`first-pytest.txt`). Diagnosis: the selected
phase fixture had zero averaged pressure-Hessian entry, so its
missing-pressure-density mutation was insensitive. Replacing that fixture
with a phase having a nonzero cross-pressure entry repaired the TEST;
no API equation changed. `repaired-pytest.txt` records9 passed in3.51s.
Adding the explicit H-minus-independent-pressure-K assertion produced the
final receipt `final-pytest.txt`:9 passed in3.50s, exit0. Every run is
preserved; no rerun was made merely to recreate an old stdout.

Final checks: Ruff and scoped diff whitespace pass. Impact was determined
by explicit rg over src, tests and the new attempt, not the GitNexus index
which the parent identified as stale (c9e1c0f versus861da15). Only the new
test imports the new namespace; the attempt is a provenance consumer.
Existing0040/campaign consumers and __init__.py are unchanged. No graph
reindex or unrelated downstream replay was represented as having run.

SHA256 at this scientific boundary:

- src/substrate_framework/euler_fourier.py:
  3d5ce4ac6923e0c6282dddb3bd8773ae9f67ebc9faaa7a0779c3d171e78e8aee
- tests/test_euler_fourier.py:
  e628b1ceddf31a788ef445b47951db7e0c0d70a06ae2cdea188fac4d3546710a

route_verdict: established as stated.
evidence_scope: exact conditional importable Fourier/Euler algebra and its
material/coadjoint distinction; no EPS geometry, positive-action closure,
accepted claim promotion, or completion of the parent objective inferred.
