# 0049 — reusable joint Euler-orbit reduction and ensemble averaging

Main owns this attempt and the new euler_orbit.py API/tests. This extracts
the same-action algebra required by 0046/0048, not another microscopic
mechanism. Frozen inputs: the complete positive compact Hessian in variables
(r,q,s), nonzero body/internal KKS pairings, and a common symmetry coordinate
B with zero energy row. The physical angles are B and q; both reaction
momenta remain independent until variation. No prescribed kinetic mass is
introduced.

Derive the exact Schur reduction including every mixed Hessian entry and
the resulting gyroscopic q Bdot term. For the time-reversal-paired ensemble,
tie only physical B,q, not the conjugate reactions (r_plus,s_plus) and
(r_minus,s_minus). Average the full unreduced actions, vary all four momenta,
and verify that the positive kinetic term survives while the odd gyroscopic
term cancels. Prematurely tying the momenta is an explicit mutation: its KKS
term cancels and it does not give the claimed inertia.

The final covariant affine-cage map uses the complete reduced kinetic
matrix, not the diagonal special case. Its nonzero denominator is stated
and checked; positivity follows from the underlying Schur complements.
Exact algebra, variation and source/test comparison are the oracle. There
is no numerical remainder, empirical comparator, or spectral tolerance.

API impact: new symbols were absent from direct source search and the
GitNexus index. Existing collective_coordinates.py only supplies scalar
pullback/reparameterization and is unchanged. New APIs are conditional,
unpromoted infrastructure; they do not themselves construct EPS data or
assert that a finite-dimensional constrained orbit is PDE invariant.

The importable reduction is now implemented in euler_orbit.py. Four tests
pass in pytest-first.stdout.txt, including independent variation of the full
mixed Hessian and of the unreduced paired ensemble. The wrong tied-momentum
ensemble is explicitly shown to lose every time-kinetic term. Ruff and
diff checks pass. This API supplies algebra for 0048 and the continuing
joint gradient reduction; it does not replace their microscopic proofs.
