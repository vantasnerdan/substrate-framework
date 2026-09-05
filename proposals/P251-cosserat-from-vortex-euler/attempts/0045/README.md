# 0045 — positive physical core-angle action directly in an EPS tube

Parent: P251 / issue #198, original conditional affine smooth-Euler objective.
Owner: Codex `/root/construction_review`; write surface: this directory only.
0040 is frozen. This attempt uses an actual EPS stationary Beltrami tube,
rather than transferring the periodic prototype's complete geometry.

## Frozen object and positive deliverable

Take a smooth stationary Beltrami field u0 on R³ supplied by EPS 1210.6271,
with curl u0=lambda u0, lambda nonzero, and an invariant material solid torus
D. Select an interior point x0 with omega0(x0) nonzero and a small coordinate
ball compactly inside D where omega0·e_z has one strict sign. Orient e_z so
that sign is positive. Constant density rho>0 and the ordinary Euclidean
Leray projection on R³ are retained. The background may have infinite total
energy; the compact isovortical tangent energy Hessian is finite.

Construct a compact physical core-rotation generator in an inner ball and
a disjoint compact negative-helicity cage pair in an annulus. Derive the
FULL two-coordinate Euler orbit Hessian and KKS form, including all cutoff
and projection terms. A finite-wave-number bound must prove strict
positivity/nondegeneracy before eliminating the conjugate cage coordinate.
The desired result is a positive ANGLE action carried by this actual smooth
stationary EPS tube, with exact finite-parameter integral coefficients and
no fitted modulus or imposed external reservoir.

Candidate: compact curl-generated circular polarizations with carrier signed
wave number k having the same sign as lambda. The microscopic core angle is
attached through a separate rigid-rotation jet. Structural selection uses
smooth Euler, exact volume preservation, physical angle observability,
same-tube reaction, full positive Hessian and nondegenerate symplectic form.
No empirical comparator or macroscopic target number enters the construction.

## Analytic design and scope

Use delta u=P(xi cross omega0), the physical orbit energy Hessian
rho integral[|delta u|²−delta u.curl(delta u)/lambda], and the KKS sign
derived in 0040. The irreducible operation is an exact Leray projection,
not a discretized PDE. Fourier-symbol inequalities bound its oscillatory
remainder in L². All positive thresholds are derived from norms of the
declared smooth background and cutoffs before any numerical value is used.
There is no production numerical spectrum or stability calculation.

The two-generator ensemble is an explicit constrained microscopic family,
consistent with the original conditional coarse-graining scope. The result
does not demand that every microscopic Euler solution remain in that family.
Core/cage forces are the derivatives of the same Euler orbit energy; no
independently prescribed strain source is admitted. Common-angle/macro
kinetic coupling and spatial homogenization remain parent obligations.

The construction and exact finite-norm theorem are in `proof.md`.
`verify_direct_eps.py` independently derives its local vector identities,
carrier helicity, integration-by-parts sign, observable core tilt, and full
off-diagonal canonical elimination. The analytic inequalities, not a
numerical eigenvalue tally, are the positivity oracle.

Route verdict: established as stated (individual review pending). The parent
campaign remains active for common-angle/translation and spatial closure.

Verification: `stdout.txt` records 17/17 exact checks; Ruff passes. The
first execution's identities all passed but its ledger completion call used
the wrong API. `first-run.txt` preserves that implementation failure; the
repair uses the repository's `CheckLedger.finish()` without changing any
equation. There is no numerical stability window or unmeasured small ratio.
