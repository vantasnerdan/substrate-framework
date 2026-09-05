# 0075 — actual parcel centroids and compact induced-velocity candidate

Parent P251 / issue #198; owner `/root/orientation_construction`; this
directory only. Parent objective remains the full physical U,Phi Euler
construction, with the original slow-affine Cauchy--Born closure.

Frozen failure-generated obligation: after the 0072 mean-current audit,
seek a materially distinct physical coordinate construction. Use actual
material-parcel mass centroids and retain intrinsic angular momentum,
rather than renaming the Eulerian point mean. Derive any one-fluid mass,
internal kinetic, symplectic and pressure/current balance from the same
material action. No assigned rotor mass or desired kinetic energy.

Candidate A: exact finite material partition, with compact isovortical
internal tangents centered in each parcel and all ambient pieces counted.
Candidate B: construct actual compact induced velocities
`delta u=P(xi cross omega)` with nonzero angular moment and a physical
core-angle jet, thereby eliminating induced exterior parcel means.
Candidate C: if direct compact inversion fails, use a local covariant
potential or constrained displacement representation with explicit
pressure reaction and a computed centroid-to-point-mean multipole map.

Selection: exact Euler/Kelvin admissibility, real parcel observables,
nonzero retained physical contrast, finite support or explicitly retained
exterior motion, and full kinetic/canonical bookkeeping. The 0072 failure
is a route mechanism, not an obligation-level no-go. Analytic identities
and exact symbolic checks precede numerical work. Status: frozen.

Result: `material-partition.md` proves the complete reference-parcel
mass/cotangent/pressure identities and full-operator estimates
`P-P_centred=O(R^-7)`, `S_parcel-S_global=O(R^-1)` for the actual compact
sources with their retained Leray tails. This supplies a nonzero physical
parcel-spin row and a conditional mass-spin normal form retaining affine
inertia. It does not yet construct a time-homogeneous material shape
ensemble or the complete shape-orbit reconstruction. `compact-velocity.md`
proves the constant-vorticity compact-velocity route has zero angular
moment despite possible nontrivial core rotation jets. Both are scoped
results, not an obligation-level verdict.

The first exact verifier run passed 16/16 checks; its unchanged output is
in `first.stdout`. The checks concern the derived identities and an exact
rational full-operator oracle, not numerical Euler existence or stability.
Continuation: the actual invariant EPS tube/ambient material-phase route
is registered separately in 0080. Earlier frozen records remain untouched.
