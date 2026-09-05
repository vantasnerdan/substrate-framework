# 0087 — an actual generalized-Lagrangian mean and its pseudomomentum

Parent P251 / issue #198; owner `/root/construction_review`; this directory
only. Frozen alternative to finite-parcel collar and compact-velocity
construction: define a genuine volume-preserving flow factorization
g=Xi composed with g_mean and an explicit mean/projection at fixed material
labels. Determine whether its physical mean coordinate realizes the coupled
canonical displacement suggested by 0066, with the actual core angle varied
and full pseudomomentum, mass density and Kelvin constraints retained.

Primary import: Holm nlin/0103043, especially Theorem 3.3 and sections 3–4,
as a Lagrangian-averaged Euler–Poincare framework, not a closure theorem.
The warning about prescribing fluctuations rather than varying Xi is part
of the construction. Internal inputs: 0059 full relative-angle action,
0072 exact point-mean cancellation, 0082 material tag transport.

Candidate routes: (A) a mass-preserving near-identity GLM/Karcher or explicit
projection mean, with a jointly varied physical angle and derived momentum
map; (B) a volume-preserving coordinate factorization constrained by a
specified zero-mean fluctuation gauge. Selection is an actual averaging
operation, correct physical observable, retained same-action reaction and
no supplied inertia or coupling coefficient. No empirical comparator,
all-wavelength invariant-manifold condition or law asserted as closure.

Oracle: primary theorem hypotheses, exact factorization/variation identities,
second-order material expansions and the complete cotangent/kinetic map.
The mean definition must select the proposed coordinate independently of
the desired Cosserat equation; a formal variable rename alone does not meet
this route's positive intent.

## Completed route receipt

`volume-preserving-glm.md` constructs an actual mean, not a variable rename:
the arithmetic material mean is corrected by an explicit physical-space
Moser density flow. The result is idempotent, volume preserving and right-
equivariant under common volume-preserving relabeling. Its exact factorization
retains pressure, the mean gauge, all fluctuation variations and individual
Kelvin data. The full kinetic Gram gives both mean pseudomomentum and the
independent fluctuation equations.

Direct calculation derives the leading observable map
`u_E-u_mean=curl(S)/(2rho)-P div(Cov_t)/2`, with amplitude/spatial remainders
explicit. Thus a calculated physical fluctuation spin law can realize the
normal-form displacement through an actual averaging operation. No old
coadjoint inertia is inserted into that spin law. The new genuine positive
material-Jacobi mode is the next input; 0091's projected-rotation material
family supplies an exact spin/metric mechanism while retaining its physical
core registration, covariance evolution and Kelvin reduction.

Primary source read: Holm nlin/0103043v1, full sections 3–4 (including the
delta-Xi omissions, density and pseudomomentum formulas, and unclosed mean
equations), at https://arxiv.org/html/nlin/0103043v1. It is imported as a
variational framework, never as the missing fluctuation closure.

Route verdict: established for the explicit mean and full action/observable
identities. Evidence scope: exact conditional GLM averaging license and
computed retained-order map, not completed positive EPS material closure.
`verify.py` first execution (`stdout.txt`) passes 17/17 exact checks, exit
zero, 1.698 seconds; Ruff and scoped diff checks pass. Tests expose a wrong
spin-current sign, omitted pseudomomentum, and mistaking a positive scalar
fluctuation mode for nonzero fluctuation spin. No numerical small-ratio or
empirical comparator is used.
