# 0082 — invariant tube centroids with a continuous material ambient

Parent P251 / issue #198; owner `/root/construction_review`; writes in this
directory only. Frozen exact theorem task: start from ONE material Euler
fluid, partition its reference labels into bounded invariant tube domains
D_a and their continuous ambient complement, and derive the hybrid momentum
and same-action cotangent/kinetic reduction. Tube centers are actual mass
centroids; ambient fluid is not divided into fictitious finite parcels.

Define hybrid momentum as tube-centroid momentum distributions plus ambient
point momentum. Retain the tube spin/shape multipoles and all tube/ambient
relative motions. Under an explicitly stated slowly varying common-translation
ensemble, calculate the difference between pointwise mean centering and
hybrid centering and the resulting physical action/field map. Determine,
rather than assume, the role of actual tube spin
A_D(S)=rho integral_D (x-X_D) cross P(S cross omega).

Inputs: 0072 exact point-mean cancellation, 0073 parcel moment identity,
0078 centroid/cotangent and Kelvin reduction. 0080 independently supplies
actual invariant-domain moment-response rank; it is not presumed proved here.
Original smooth-Euler, full pressure/ambient reaction, affine conditional
scope, and fixed circulation conventions remain unchanged.

Oracle: exact material calculus, distributional Taylor expansion with its
remainder, symplectic/metric block projection and full operator Schur algebra.
No isolated-cell inverse, added rotor mass, all-k invariant ansatz, fitted
coefficient or empirical comparator. Shared-face continuity and pressure
tractions stay in the single action; neither a tube-centroid field nor the
hybrid field is silently identified with the whole-fluid Eulerian point mean.

## Executed result

`hybrid-action.md` derives the exact hybrid first-multipole map, same-action
cotangent/kinetic centering and physical fixed-angle affine lift. Actual tube
spin, not global source angular impulse, controls the difference from the
point-mean cancellation. Under an admissible material lift and full reaction
reduction, hybrid centering leaves the absolute tube-angle source; all shape,
ambient relative motion and full reaction inverses remain explicit.

The known frozen-tag defect is repaired at the observable level by solving
the actual surface transport equation. The derivation supplies its exact
retarded boundary-spin history and the necessary-and-sufficient instantaneous
collar criterion `f0=L f1` for the two velocity columns. A missing material
commutator column produces a frequency-dependent boundary spin; including
that actual column removes the transport poles exactly. A compatible
material collar displacement supplies the column, and the coordinated 0084
calculation owns its full Jacobi and Kelvin-reconstruction action. This
continuation is executed, not replaced by a generic boundary caveat.

Route verdict: established for the hybrid identities and complete tag
transport/locality criterion. Evidence scope: exact conditional action and
observable bridge; full positive two-field EPS assembly depends on the
computed admissible collar/reaction block rather than moment rows alone.
The parent remains active; no completion or exhaustion is inferred here.

First verifier execution, preserved in `stdout.txt`, passes 20/20 exact
checks, exit zero, 1.763 seconds. It imports canonical Schur-jet/ledger
helpers; Ruff and scoped diff checks pass. Both helicities, a wrong point-
mean substitution, a missing commutator column and resonant transport are
exposing tests. No numerical small-ratio design or empirical value is used.
