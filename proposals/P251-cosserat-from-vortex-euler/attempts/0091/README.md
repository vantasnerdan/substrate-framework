# 0091 — positive material core angle with exact complete-fluid spin metric

Owner `/root/smooth_core_review`; constructive P251 / issue #198 subtask,
this directory only. Parent objective and conditional slow-affine scope are
unchanged. Accepted base v0.171.0; campaign inputs 0084 full material
transport/Jacobi/spin, 0089 exact positive Fourier trial, 0090 compact
positive cage on one EPS-containing Euler field, and 0087's defined
volume-preserving material mean. No old coadjoint inertia is imported.

Positive target: a same-field material lift with actual core rotation,
positive stiffness, and equality of its complete physical spin-rate row
and material metric; derive its full common/relative-angle kinetic and
potential blocks before isotropic averaging, with ambient and Kelvin
reaction terms explicit. A tube-only angular impulse is not substituted
for complete fluid spin. No fitted elastic or inertia constant is allowed.

Candidate A: construct a radial rotational core lift and use a disjoint
positive material cage to satisfy the spin/metric moment equation exactly.
Candidate B: projected rotation of a divergence-free fluctuation covariance
for 0087's GLM registration. Candidate A is selected first because its actual
core jet and complete ambient spin are explicitly controlled. Candidate B
remains a possible mean-observation join, not a presumed identification.
Selection uses exact geometry, same-action coefficients and pressure/Kelvin
bookkeeping, with no empirical comparator or numerical eigensolver.

Oracle: analytic scaling bounds and exact polynomial/Schur identities.
The constructed mass equals an actual velocity-space moment, not an
assigned rigid rotor. Positive restricted material action alone is not
asserted to preserve each Kelvin leaf or form an invariant nonlinear Euler
ansatz. The exact circulation equation is retained and its status must be
decided explicitly in the completed statement.

## Completed construction

`joint-material-action.md` constructs `Z=B_L+tX`, with the explicit positive
root of its spin/metric equation. It proves `A(Z)=M(Z)=j>0` and `K(Z)>0`
at a finite coherence scale while preserving the core rotation jet. The
complete-fluid physical angular momentum, common/relative angle and all
affine kinetic terms are derived from the same material action. The actual
mean material shear is retained with its sign; 0094 supplies its independent
positive STF population, and 0093 supplies positive gradient stiffness with
the complete gradient mass. Neither old coadjoint inertia nor old affine
shear is copied.

`verify.py` passes 14/14 exact checks on its first run; `stdout.txt` preserves
the execution. The checker corroborates the exact moment, scaling,
kinetic-completion and same-field stiffness identities; the prose supplies
the field construction and its finite bounds. No numerical sign inference.

`route_verdict: established`; `evidence_scope: REPRESENTATION_SCOPED` for
the constructed positive moment-normalized material action. The actual
registered GLM observation and full individual Kelvin/reaction reduction
are explicitly distinguished from this conditional material pullback.
Their completion remains an active parent construction rather than an
assertion earned from a positive metric alone.
