# 0123 — actual rotating-profile mechanical/canonical moment matching

Owner `/root/construction_review`, this new directory only. Base release
v0.175.0; accepted conditional claims unchanged. Parent objective remains
the stronger actual-Euler finite-time translational/angular realization.

Frozen positive theorem target: derive the mechanical tagged-spin row,
physical core-angle map, and canonical action of an axial-profile
displacement on exact uniform rotation; construct a smooth axial tail
whose derived profile integrals match the mechanical and canonical
inertias. Then derive the actual compact divergence/pressure completion
error and its transfer to the local gyro-dominant torus of 0120.

Route A is a nearly linear core plus a disjoint scaled odd axial tail,
using the exact quadratic moment equation and its simple small root.
Route B, activated by any completion error, uses finite-parameter
implicit continuation of the complete measured action/observable row.
Selection criteria are exact Euler/Lin compatibility, complete material
spin including shape/boundary terms, physical core normalization,
positivity, and controlled finite-time localization errors. A chosen
profile may be a constructive input; no measured or fitted frequency is.

Requires: original Euler material/Kelvin identities, 0119's physical spin
diagnosis, and 0120's fixed-curl torus approximation. Pass licenses:
profile matching and whatever explicitly bounded continuation follows.
It does not license all-time equality on a perturbed torus, an invariant
finite-dimensional Euler ansatz, or an autonomous coarse constitutive law.
Success activates the parent's same-action physical-observable join.

Oracle: analytic calculus and exact SymPy coefficient/integral identities,
with sign, factor-two, and omitted-completion mutations. No empirical
comparator, floating-point spectrum, numerical discriminant threshold,
or small-ratio solver is used. The profile root and hierarchy are chosen
from exact inequalities before any illustrative rational integrals.
Irreducible analytic question: whether a single simple-root parameter
controls the entire finite-time spin/action discrepancy, or only its
reference scalar row with a uniformly bounded temporal remainder.

## Frozen result and oracle receipt

`profile-match.md` establishes the smooth core/tail root, its exact
compact pressure-completed reference-time mechanical/canonical match,
and a uniform finite-time RELATIVE spin/shape error for actual linearized
Euler on uniform rotation. The proof uses the off-diagonal kernel of the
complete rotating Euler propagator, not an L2 estimate multiplied by an
uncontrolled radial lever arm. The remaining scalar spin connection and
time-dependent inertia mismatch are explicit. Global torus/EPS transfer
is delegated to the new disjoint attempt 0124.

`verify.py`: 24/24 exact checks, process exit 0, 2.397 seconds; captured
on its first scientific execution in `first-run.txt`. It derives the
divergence/curl completion, signed Kelvin residual, full KKS normalization,
material spin and nonzero shape moment, profile integrals, simple root,
physical moment-map determinant and full projected inertial-wave symbol.
Wrong source sign, missing divergence completion and missing factor two
are exposed. No numerical small-ratio test or empirical comparison occurs.

The integral fixture uses piecewise-polynomial C2 profiles only to check
exact integrals; the theorem uses smooth bumps with the same strict
margins. Smoothing and the root continuation are argued in the proof,
not inferred from a discretized eigenvalue.

The first Ruff pass identified two assigned-lambda style violations.
They were replaced by identical local function definitions, with no
scientific expression or predicate change. The corrected Ruff pass is
recorded below. Whitespace check reports no errors; the no-index diff
exit 1 denotes new files. No canonical API, accepted claim, or shared
parent file was changed. The concurrent 959db88 full validation belongs
to the parent and is not rerun or represented as this attempt's oracle.

Route verdict: established at the explicit uniform-rotation,
reference-match plus finite-time-error scope. Evidence scope:
`EXACT_PROFILE_AND_REFERENCE_MOMENT_MATCH_WITH_CONTROLLED_FINITE_TIME_
UNIFORM_ROTATION_EULER_ERROR`. This is not all-time moment equality.
