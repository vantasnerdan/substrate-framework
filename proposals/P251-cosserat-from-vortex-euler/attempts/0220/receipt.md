# 0220 frozen construction receipt

Author: `/root/construction_review`. This is a new construction,
not its independent claim review. Parent confirmed central/schema
activation with 269 accepted claims before substantive work. The
accepted base remains v0.181.0; current provenance when freezing is
4bae7090c9d068423517c160bed0d1bb9a2d95a1. The parent campaign and
all unchanged dependencies retain their own authority/status.

## Route verdict and proposed meaningful result

`route_verdict: established`

`evidence_scope: proposed actual positive-energy low finite-harmonic
Euler bending pole in the invariant compact-active-vorticity sector
of the same smooth global ring; individual review pending`

Candidate A closes the curved-reaction route as stated in the proof.
For the fixed new 0211 smooth profile and each fixed n>=2, sufficiently
large finite R has a simple positive-frequency bending pole with

    omega_n²=[Gamma_delta log(R/a_ref)/(4pi R²)]²
                     n²(n²−1)[1+o(1)].

Gamma_delta is the fixed limiting seed circulation; finite-R
circulation and finite-core corrections are included in o(1). The
mode is an actual smooth full-Euler velocity eigenfunction with the
complete exterior pressure. Its isolated-contour assertion is in the
declared invariant active-vorticity sector, not a gap for arbitrary
independent exterior passive-vorticity perturbations.

The inherited phase energy is positive, including the eliminated
fast response and its frequency derivative. The actual n=2 global
material quadrupole is nonzero on the same inner invariant torus.
Its linear global material-spin vector is EXACTLY zero by symmetry.
The circular reference covariance is not assigned a fictitious
nonzero planar eigenvalue gap. The n=1 Euclidean Jordan family stays
separate and unchanged. These distinctions leave the parent
microrotation/current/common-continuum obligations active.

Candidate B, an optical-window normal-form history without a pole,
was not selected because candidate A supplies the stronger actual
spectral object within its declared sector. This is not exhaustion
of response-specific alternatives or completion of the parent.

## Load-bearing construction sequence

1. The exact toroidal Piola, force, pressure, kinetic and KKS
   pullbacks retain their actual volume factors.
2. Actual material action-angle coordinates give a unit volume
   form and the precise covariant identities a_y=G(H),
   (average a_theta)'=G'(H)W. The first averaged curvature vanishes
   by the complete centered m=1 equation; second-order means remain.
3. The full mean-metric forced Sturm equation includes angular-axial
   coupling and its actual particular pressure. Its positive-sector
   form is a relative perturbation of the Hardy form: the exact
   cross metric carries G, and the J' term is integrated before it
   is estimated. The full flat-tail cancellation survives.
4. Compact pressure-gauge subtraction normalizes a force by its
   PHYSICAL velocity, uniformly at small axial k, while retaining
   the ikC axial current. The explicit scalar R³ ring Green kernel
   then controls the actual exterior return. No globally small
   toroidal metric is assumed.
5. In the same volume variables both metrics have the same steady
   background u and omega. Their Kelvin difference acts only
   through the physical Leray response, avoiding the spurious large
   gauge-parallel generator. The reference fast inverse plus this
   small pressure comparison sums the curved feedback.
6. A finite augmented KKS system changes from the reference
   translation complement to the EXACT geometric bending complement.
   The actual 0213 residual makes its correction small at the
   low-carrier scale. The full self-energy, and by Cauchy its
   frequency derivative, are lower order than the bare phase rows.

This sequence, rather than the number of checks, is the proposed
proof. The scalar Green comparison, full fast-domain estimate and
exact geometric projection are explicit load-bearing review targets.
No empirical comparator, assigned mode, cell wall or fitted
coefficient enters it.

## Exact execution record and scoped validation

Both scripts use the existing `substrate_framework.verification`
CheckLedger. The strongest oracle is the analytic operator/form/
kernel argument. SymPy verifies its exact full rows, pressure
eliminations, moment identities and action derivatives; it does not
claim to compute a numerical spectrum.

- `verify_metric_sturm.py`: first execution passed five predicates,
  then failed its axial pressure predicate because that predicate
  omitted the two background curl substitutions already used in
  the displayed physical velocity. Original failure and traceback
  are preserved in `first-metric.stdout`. The one-line predicate
  correction is recorded in `first-metric-repair.md`. No scientific
  equation or target changed. Corrected execution passes 12/12,
  process zero, in `repaired-metric.stdout`.
- `verify_pressure_feedback.py`: first execution passes 13/13,
  process zero, in `first-feedback.stdout`.
- Both final scripts have respectively twelve and thirteen lexical
  check calls and runtime executions. No loop count is conflated.
- Scoped Ruff passes for both scripts. No old scientific verifier,
  full suite, registry, release or shared API was rerun or edited.
  An untracked-file `git diff --no-index --check` probe on the main
  proof emitted no whitespace diagnostics (status one is its normal
  new-file difference status). Integration owns its final staged
  whitespace check; no broad check is claimed from untracked status.

Commands used the existing interpreter
`/home/dan/substrate-framework/.venv/bin/python` with `PYTHONPATH=src`,
pipefail and first-run tee capture. There is no numerical soft-mode
design, quadrature floor, discretized eigenvalue or convergence gate.

The proof now goes to one independent review before any canonical
claim promotion. Reusable metric/pressure/Schur definitions, if
accepted for promotion, still earn their extraction and direct tests
under the parent's implementation boundary. No canonical authority
or parent completion follows from this construction receipt alone.

## Frozen SHA256 boundary

- README.md: `170764bc18fd3751febf616fc9111e63f84c9696d9d19ac557136619672d0172`
- actual-operator.md: `9c010c7f45a547aea73aeb6826f807a80f8d623cf9fe899214918699ba35f641`
- mean-metric-response.md: `b9c45edd06e2edf036b690394cead12db890ef13545d9ab936dedccf68fb1eaa`
- pressure-gauge-and-feedback.md: `de1e37e1c71c92804a28d08202c49b7cb8f88943d6c67a41933b0a6bcb182245`
- verify_metric_sturm.py: `b4512d7b6d44dc60013c3a753c07cd7802957d9105010a09b0fcbc0f62e44246`
- verify_pressure_feedback.py: `66641e5b1d80aef63f7e8bd81b443c1a06856523edc634b743e20264aab16afe`
- first-metric.stdout: `fe3535db717554d7a3b567a078b74ba47168c998bb5b4d67beb811b90b2ab3b4`
- first-metric-repair.md: `9a246f5c274cd79f7e23b095bddfb34afb28011e9c2f3d0d03010a5eb9ae830b`
- repaired-metric.stdout: `0e3aaa5dc710bec438db7a00b611b765454e7b32f541580741ba830a6ced8735`
- first-feedback.stdout: `1fcde794407324fe561dbe0352c0f134b85a90a28befc55e36b869560d9c16d7`

Unchanged direct inputs:

- 0213 receipt.md: `6c466193740693088b674d11ee65f4cf23970d25ce480191a023f80e46cfcb12`
- 0211 affine-inner-seed.md: `a477727fa9a5d90becf2db8dea7b62882eea19ac6ca56f5ef421b50d3ed3dc52`
- 0211 same-ring-beltrami-torus.md: `8c2fa0c980c740e4d6acd55da324b66ac4a93c9a2862f0c9b2a816cf553291a6`
- 0206 low-harmonic-action.md: `89adfda0260d49b1c4ad33c7b0e1c2d0d534e7d97bb0793bb0a1553c82358184`
- 0206 parity-threshold-response.md: `4d08ff90ca7b42ac22d435596b49a012d28b1e4f4df8a82f972389bbcce954f9`
