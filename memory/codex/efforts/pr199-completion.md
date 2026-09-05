---
description: Complete the Euler to Cosserat construction and repairs for PR 199
author: codex
created: '2026-09-04T18:00:00+02:00'
updated: '2026-09-05T03:00:00+02:00'
tags: [substrate-framework, p251, pr-199]
category: efforts
confidence: working
status: active
---

## Goal and Success Contract

Meet the full positive objective of issue #198 and PR #199, preserving the
seven-node objective and explicit Euler, isotropy, and ensemble premises in
P251/proposal.yaml. Completion includes the microscopic action construction,
individual claim review, reusable implementation, and release promotion.

## Accepted Baseline

Frozen review head: 10a0f31bf0122f4c744a81dfe1cfa66506d57ad7. Campaign base:
v0.171.0 at e24497b0272f8582a71d4189581db818a0e276aa. The separate main
worktree pins v0.174.0 and has a user-owned AGENTS.md change; work here uses
research/pr199-completion in an isolated worktree. P251 has no accepted claim.
P242 homogenization is conditional imported infrastructure. Its sphere moments,
line tension integral and affine matching were read at source.

## Plan and Obligation Graph

Inherit N1–N7 and L1–L5 from the validated proposal, without reducing the
objective. N2/N3 are active. First repair the Euler momentum correspondence;
then construct angle-dependent interaction and its reduced action; replay
N3–N7, review each proposed claim and promote supported full results.

## Candidate Set and Analytic Receipt

Retain filament triads, core polarization and ergodic contrast. Compare exact
relative-angle variation, dimensions, Euler dynamics and assumption cost
before values. No empirical comparator. Audit the latest rate-quadratic
polarization candidate; independently explore nonaxisymmetric patches in
strain and finite vortex pairs. Exact algebra is the oracle for this stage;
there is no numerical remainder yet. A frequency squared alone supplies
neither an elastic potential nor its canonical coordinate/kinetic metric.

## Ownership and Impact

Codex owns repair implementation and attempt 0029+. An independent analytic
worker explores orientation interactions read-only, following AGENTS.md's
parallel replacement-candidate instruction. Existing attempts stay immutable.
GitNexus does not index P251 in this worktree; direct consumers are the seven
proposal verifiers and their new attempt receipts. No accepted API changes yet.

## Attempts and Validation

Before source-level execution, validate_repository.py passes: 252 accepted
claims, 12 proposals. Physics preflight passes. Previous valid receipts remain
evidence for unchanged scopes. New derivations will carry source and first-run
stdout. Repaired statements receive exact residual, dimension and mutation
checks; numerical design will follow its own analytic receipt if needed.

## Results and Remaining Achievements

Active investigation: attempt 0028 varies angular rates, while N4 requires
angles. Its coefficient map must preserve dimensions and the variational
equations. Independently audit attempt 0019's momentum residual operator before
using its allegedly corrected fields. The source PR remains active refactor.

Established in this run: 0029 derives the Cartesian Euler operator and restores
the original N2 velocity formulas; 0031 repairs the K1 derivative and derives
the m=1 constant 1/4-EulerGamma, corroborated by 40/60-digit signed roots.
0028's alpha_E has inertia units and cannot supply N4's elastic angle torque.
0033 repairs the mutual kinetic-energy kernel and derives a positive local
opposite-circulation tilt interaction, reopening the 0026 route.

Executed continuation: 0030 derives a genuine angle oscillator from an exact
elliptic Euler patch with prescribed ambient flow. 0032 supplies a closed
six-point-vortex relative-angle action, positive stiffness and positive
inertia, corroborated by a separate full Cartesian Jacobian. 0034 derives a
covariant collective-angle field map, one consistent alpha/j pair, conditional
helical-gradient wryness, and equality of the generalized and standard
transverse dispersion through order k^2. Cell density is L_cell=L_v/6.

The exact smooth-Euler continuum objective remains active. Finite-core
continuation, the reaction-bearing affine cage/translation closure, and
stationary EPS compatibility are not supplied by a point-vortex oscillator.
The user reaffirmed the exact smooth-Euler objective and explicitly directed
continued execution beyond checkpoint commits. No asymptotic scope reduction
is authorized. Continuation contract 0035 freezes finite-core/spectral,
material-action and stationary/EPS routes; worker 0036 and main 0037 execute
in parallel. The original objective remains the promotion boundary.
The new action gives a concrete candidate to review, not campaign success.

Continuation 0035–0043 (2026-09-05): 0036 constructs an actual C-infinity
finite-core Euler polygon by a bordered implicit-function argument, and its
exact two-triangle affine KKS angle action. Independent review establishes
that stated construction; it does not supply a stationary EPS embedding.
0037 derives the complete Euler displacement/Jacobi action with transport,
pressure Hessian and boundary reaction. 0039 derives material parcel mass,
locked inertia and a separate exact screw/Galilean stationarization identity.
These definitions are now conditional unpromoted euler_displacement.py APIs.

0038 tests actual stationary tube candidates. Its static Jacobi angle has
negative stiffness, motivating the representation change 0040 rather than
closing N3. 0040 derives a positive two-generator isovortical angle/cage
action on an exact stationary nonaxisymmetric Beltrami tube. Its physical
angle is the rotation of the quadratic core jet, not a rigid rotation of the
whole finite cross-section. Independent review establishes this scope and
adds a flat-torus Leray-transfer proof. The prescribed-knotted-tube bridge
is not inferred from mere local approximation.

0041 derives the full finite-core 3D Biot–Savart gradient matrix, including
the second vorticity variation. Its transverse-projector Gram form proves
positive twist and radial-momentum coefficients. Keeping the conjugate
radius gradient changes the optical k² coefficient: for the point-limit
affine-cage map the leading coefficient is 1129/220 times the old twist-only
value. Earlier 0034 twist-only numbers are superseded by this correction;
finite-core values are exact integral matrices, not this imposed ratio.

0043 independently derives positive affine shear energy from the same
stationary Beltrami fluid and its deformed periodic lattice: mu=4E0/15,
where E0 is the actual unstrained fluid energy density. Its new importable
euler_affine.py API inverts curl after pushing vorticity and wave covectors
together. No additional filament spring is postulated. This is an energy
coefficient, not yet the full joint kinetic reduction.

Active executable routes: 0042 derives common-angle material/orbit reaction
and repairs the distinction between a periodic spatial RVE and an invariant
material parcel. 0044 computes the stationary Beltrami sector's full Bloch
action, including its momentum/kinetic gradient terms before isotropic
averaging. 0045 constructs a physical core-tilt/cage sector directly inside
an actual EPS tube using compact negative-helicity generators. These are
campaign work on research/pr199-completion, not a terminal PR refresh.
The full common-frame/translational action, all N1–N7 closure, individual
claim promotion and release still determine completion.

Further completed constructions: 0045 is now independently reviewed as an
actual-EPS positive physical core-tilt/cage theorem, with explicit finite-k
norm bounds and exact integral coefficients. 0046 establishes the common
rotor intrinsically on one Euler orbit in an exact nonaxisymmetric Bessel
cylinder; review supplies the finite-domain Neumann estimates and requires
momentum-partner support away from the physical core. Both inertias derive
from KKS/Routh, avoiding 0042's Kelvin-inconsistent added-body-inertia route.
Next: 0048 extends this common rotor to the actual EPS relative-energy orbit.
0047 executes the same-EPS spatial lift; its initial carrier-only curvature
is negative, so explicit geometry-tied first-gradient cage return is being
constructed rather than averaging frequency squares into a false positive.

Checkpoint validation: fixed repository checks and the full 2501-test suite
pass (351.84 seconds), after the changed-scope selector's conservative
multiple-new-modules rule required full replay. The initial targeted nine
API tests also passed. Source PDFs are marked binary in their source folder;
this preserves byte-exact digests and avoids interpreting PDF xref padding
as source-code whitespace. No accepted registry or release changed.

## Canonicalization and Completion

Continuation 0047–0058 is now substantially constructed and independently
reviewed. 0048 gives a positive full common/internal action on the actual EPS
relative-energy orbit, with mixed Hessian terms and independently varied
time-reversed fluid reactions. 0049 extracts its reduction into euler_orbit.py.
0047/0050 keep finite-coherence gradient returns, full momentum reduction and
exterior kinetic energy. 0051/0052 derive exact material force/couple balances,
centroid kinetic bookkeeping and the physical-to-canonical boundary-current
map. 0055 proves a nondegenerate exterior response Gram and exact full KKS
centering with four disjoint response supports; its ambient impulse is retained.
0054 constructs positive shear from that same EPS field rather than copying
the periodic 0043 coefficient. 0053 derives the full simultaneous-rotation
average and same-map kinetic/potential gradient normal form.

The full-fluid assembly is still active, not inferred from a finite tube's
mass or from a whole-space relative energy. 0057 now executes both exact
material-partition accounting and a source-backed stationary Gaussian
Beltrami ensemble (2006.15033), whose finite-variance volume measure and
positive-density knotted invariant tori can supply the missing normalization.
0059 independently investigates its origin-independent common-rotation KKS
density. This new route reconstructs local coefficients and symmetry action;
single-decaying-field coefficients are not silently transplanted.

0056 found and repaired an additional old operator defect: longitudinal spin
has coefficient 2(c_s+c_tr), not 2(c_a-c_tr). The new unpromoted micropolar.py
action and full Fourier matrix are tested by direct real-field variation.
CST004/005 now use the same-orbit inertia/locking map, with 8/17 checks passing,
including DOP853 refinement. 0058 strengthens N6 through a full nonlinear
product-Haar energy argument: conservative coherent response vanishes while
positive fluctuations remain; correlated uniform marginals expose why
one-point isotropy alone is insufficient. A zero angular current is still
a separately named dynamical closure premise.

No claim promoted; no release change. Full completion and scientific exhaustion
are both unearned. Review and promotion follow construction, not a green tally.

## Cross-References

Authoritative goal and reproducible artifacts:

- proposals/P251-cosserat-from-vortex-euler/proposal.yaml
- memory/codex/efforts/pr199-review.md
- src/substrate_framework/homogenization.py
- https://github.com/vantasnerdan/substrate-framework/issues/198
- https://github.com/vantasnerdan/substrate-framework/pull/199
