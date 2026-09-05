# 0109 — a true Euler optical parcel with nonzero spin and translation

Owner `/root/smooth_core_review`, this directory only. Parent issue #198 /
P251, accepted release v0.171.0, and its complete objective are unchanged.
Uses physics-erdos-loop and the already-read material/Kelvin/source scope;
no numerical soft eigenvalue or empirical comparator is involved.

Frozen positive target: an actual transported finite core parcel in one
smooth stationary EPS-compatible Euler field, with nonzero physical spin
and translational momentum during a finite-frequency optical evolution,
controlled through one period by the full Euler equations. Preserve 0101
immutably and distinguish this parcel from the distant EPS invariant knot.

Registered candidates: (A) change the compact envelope to leave nonzero
complete packet angular impulse; (B) retain the necessary compact returns
but construct and transport an explicit parcel carrying nonzero momentum,
with the outside reaction retained; (C) spatially separated complementary
return circulation. Criterion: actual Kelvin preparation, actual material
moments, finite-time full Euler control, and explicit conservation balance.

Analytic observation motivating B: a compact divergence-free velocity can
have nonzero angular impulse in general, but the compact axial-primitive
class used for 0101 has zero impulse because its potential has zero axial
integral. This is a class-specific cancellation, not an Euler no-go. The
same cancellation does not apply to a material subparcel. The candidate
executed here uses an exactly polynomial inner envelope and the complete
moving-boundary spin; it supplies a positive matched physical inertia
without declaring the entire parcel rigid or changing Kelvin data.

Oracle: exact curl/Lin/material-moment/pressure algebra, explicit separated
norms, and a full Euler energy estimate. The naive 0101 L² bound loses too
much accuracy for the small parcel moment. An exact pressure-gradient
subtraction improves the return residual from `(ka)^-1` to `(ka)^-2` and
repairs that failure before any numerical design. All resulting scales
are finite and prescribed by the analytic error. No solver is needed.

Closure boundary: a two-dimensional optical family can have both physical
translation and spin while those responses remain linked, not independent
Cosserat fields. Global compensating returns, ambient pressure transfer,
and the distinction between transported parcel and persistent invariant
vortex structure remain explicit. This attachment claims no parent
completion or claim promotion.

Result: `material-optical-moments.md` establishes the scoped construction.
`verify.py` derives the full pressure and Lin residuals and actual material
moment/traction identities; all 21 checks passed on their first execution,
recorded in `first-run.txt`. Ruff passes. The completed analytic estimate
is `O(delta^(1/3))` for the physical parcel moments through one period,
including their inverse aspect-ratio sensitivity. No numerical solver or
unmeasured small-eigenvalue assertion was used.
