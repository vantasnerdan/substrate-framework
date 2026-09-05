# 0130 — actual finite-window spin/centroid join

Owner /root. Baseline dbc2bf3956eec22da0c64edfa42a26f04d2f1cd6,
release v0.175.0; accepted claims remain unchanged. The stronger active
object is actual same-tube Euler angular/translational response, with
the full spatial constitutive closure still pursued in0125/0128/0129.

Positive target: combine0124's actual finite-time spin normalization,
0126's paired preparation/physical angle reconstruction and0117's full
hybrid observation into a quantitatively nonzero centroid-response window.
Use a quasiperiodic four-torus lift to realize the small-amplitude
nonlinear perturbation for every rotated periodic background and fixed
slow laboratory wave vector. Compare the commensurate-supercell route
and the lift by actual Euler constraints, complete pressure/ambient,
uniform finite-time estimates and measured nonzero signal.

The exact periodic-background/Euler existence inputs and finite derivative
bounds are fixed before packet and perturbation amplitudes. The oracle is
analytic Sobolev/Galerkin existence plus exact projector/observation
identities. There is no numerical spectrum, stability edge, tolerance
fit or empirical comparator. A first-jet response is not a positive
second-gradient Cosserat pencil or a nonzero homogenized microinertia in
a vanishing-cell limit. Those exclusions retain the actual parent scope.

## Checkpoint receipt inherited without rerunning

dbc2bf3 is pushed to origin/research/pr199-completion.0127/full.stdout
passes2580 tests in361.42s and all repository checks. The postcommit
command `python scripts/validate_changed.py --base 959db88 --head dbc2bf3
--print-only` succeeds and selects FULL because euler_phase.py changed;
this confirms, rather than repeats, that receipt after the failed tree-
argument invocation. GitNexus was refreshed atdbc2bf3:45710 nodes,
71489 edges,756 clusters,2 flows; explicit rg still supplies omitted
test/attempt call edges. The scientific work below changes attempt files
only; it does not stale the frozen canonical-source/test validation.

## Frozen result

`finite-window-join.md` derives an actual nonempty finite-wave-number
spin/centroid response window, including the full ambient mass and a
small-amplitude nonlinear Euler realization on the quasiperiodic lift.
The linear normalized signal is at least 3 k j Omega / 8; the declared
nonlinear amplitude window retains at least 5 k j Omega / 16.
The geometric and finite-time constants are fixed before choosing k
and amplitude; no uniform nonzero modulus in the vanishing-cell limit
or second-gradient constitutive closure is inferred.

First execution `PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python
proposals/P251-cosserat-from-vortex-euler/attempts/0130/verify.py` passes
12 exact checks, exit0, in1.62 seconds, captured in `first.stdout`.
Ruff passes. The symbolic checks corroborate the projector, lift and
observation algebra; the Sobolev existence and error estimates are
analytic arguments in the proof, not a claim of formalized PDE existence.

SHA256 proof: `6cb1f1442ec76bc4e6ab7c9368869910f3dab6f6ab82cedf83b0b209e719b779`.
SHA256 verifier: `e4324e6ec4dc3619e541df2ba092a772789cdf944bbf46b07439b1a4e8358922`.
Route verdict: established as stated for the finite-window construction.
Evidence scope: exact algebra and analytic finite-time Euler existence
with controlled physical first-jet observations. The active next
obligation is the full coupled spatial constitutive response.
