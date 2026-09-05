# 0137 — a smooth compact-vorticity column with the actual surface pole

Owner `/root/construction_review`; this new directory only. Base accepted
release v0.175.0 is unchanged. This failure-generated continuation uses
0135's fixed nonzero carrier before smoothing, and does not identify a
patch limit with a smooth Euler mode.

Frozen positive target: a genuine smooth Euler column mode continuing
the Rankine surface branch, with positive intrinsic material-angle
action/curvature and the same nonnegative-tag mechanical/canonical
moment match. Preserve its nonzero two-radial-moment gap, full pressure,
material transport, exterior reaction and physical time normalization.
The ordinary-column result is distinct from the generalized force-free
and local-torus transfer under construction by the parent in 0136.

Candidate A: smooth the Rankine vorticity across a thin annulus, keeping
it exactly constant in the interior and zero outside a finite radius.
At the fixed carrier, its optical root is nonresonant throughout the
vorticity support. Derive a radial Evans/matching continuation on the
isovortical compact-vorticity class, with irrotational exterior response.
Candidate B, if a true-pole transfer is not established: construct an
actual finite-time transported perturbation with a controlled Euler
residual and physical moment errors, retaining rather than discarding
the resulting time horizon. This is an alternative result, not a
substitute for a spectral theorem.

Selection is structural: an actual Euler solution, controlled complete
pressure response, nonsingular fixed-carrier observation, positive
action/curvature and physical spin equality. No empirical comparator,
mode fitting or numerical design is involved. Exact calculus and radial
ODE transfer come first. Any later soft numerical oracle first reads
the small-ratio skill and records a complete error budget.

Permitted imports: frozen 0135 mode/action/moment formulas; canonical
rankine_modes equations; the parent's actual 0136 background construction
only at its proved scope; primary Gallay--Smets source after its stated
operator and scope are read. Central expansion and repository-validator
receipt precede opening that new source body or computation.

The central expansion proposed for validation is: at fixed sufficiently
small ka>0, the optical Doppler sigma is near -Omega, hence separated
from zero across a sufficiently thin smoothing annulus. The exterior
critical particle radius lies beyond the vorticity support. On the
isovortical compact-vorticity class it cannot inject an omitted exterior
vorticity mode; the full irrotational pressure/velocity matching still
remains. The exact radial transfer matrix, including vorticity-gradient
terms across the annulus, is the next derivation.

Maximum license: the proved smooth-column mode or stated controlled
finite-time response. This does not by itself license arbitrary-knot
EPS, a whole-fluid autonomous two-field Bloch band, or parent completion.
Every candidate earns its own verdict; the stronger parent remains active.

## Frozen result

Candidate A is **established as stated** in `smooth-column.md`: the
bounded first-order displacement/pressure system gives an actual
smooth-column mode by a thin-annulus matching/implicit-function proof.
It preserves the positive material-angle action and intrinsic frequency
curvature, and the fixed-carrier nonnegative-tag spin match. It includes
the full exterior potential response and the actual material relations.
The controlled-time fallback was not needed for this ordinary-column
pole. Coherent sidebands, translation, and global EPS geometry remain
separate parent constructions, not exclusions from this theorem.

The same first-order equations also include arbitrary axial flow W(r).
Their W' terms cancel exactly. The proof identifies rather than discards
the finite-tag differential axial phase in a force-free transfer.

The exact verifier passes 13/13. Its original structural comparison
failed because SymPy `factor` returned unevaluated `0*I`; `simplify`
returned exact zero without changing any equation. Both executions are
preserved. The evidence is analytic and symbolic, with no numerical
eigenvalue computation, fitted quantity or canonical/release change.
