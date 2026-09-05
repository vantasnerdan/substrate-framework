# 0055 — exterior-fluid centering and exact KKS restoration

Parent: P251 / issue #198. Owner: `/root/construction_review`; this directory
only. Attempt 0052 is preserved as the completed material/current bridge
with its then-unresolved centered common-response identity.

Frozen continuation: expand the finite coherence region to include compact
response generators in the surrounding fluid OUTSIDE the invariant EPS
tube D. All fluid is the same actual smooth Euler/Beltrami field; no
external angular reservoir or fitted coefficient is introduced. Prove that
the common angular moment is independent of the three physical cell-mean
responses using their different far-field orders. Construct fixed compact
response generators, project off those moments, and restore the complete
four-dimensional KKS structure without changing physical core jets.

Analytic oracle: the finite spherical-harmonic EPS far-field expansion,
the explicit Leray dipole of 1_D e_i, analytic continuation in the connected
exterior, finite positive Gram inverses, disjoint-support KKS identities,
and the finite high-frequency Hessian bounds already established in 0045.
No numerical eigenvalue, target modulus, or empirical comparator is used.

The physical cell is still centered at its mass centroid. Surrounding-fluid
angular momentum and its current improvement are retained when comparing
the full coherence-region canonical moment with the physical spin of D.
This is an append-only geometric/response expansion, not a change to the
parent smooth-Euler or conditional affine objective.

Status: active. Frozen before the exterior response proof.

The completed theorem and finite-support construction are in
`exterior-centering.md`. It proves the full four-response Gram is positive
on every exterior ball, builds four mutually KKS-isotropic dual responses
in separate balls, and preserves the full centered rotor KKS matrix exactly.
The explicit current map retains the surrounding angular impulse.

`verify.py` derives the dipole and checks the full moment projection, KKS
matrix, physical jets, finite positivity threshold, and ambient-current
accounting. Its first execution is preserved in `stdout.txt`; Ruff and
scoped diff checks pass. The positivity oracle is the analytic far-field
and Gram proof, not numerical eigenvalues of an unspecified EPS field.

Route verdict: established as stated; independent review pending. The
parent's material/coherence assembly remains a distinct active construction.
