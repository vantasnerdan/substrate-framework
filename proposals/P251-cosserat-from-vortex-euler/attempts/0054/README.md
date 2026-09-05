# 0054 — same-EPS finite-coherence affine shear

Parent P251 / issue #198; owner `/root/orientation_construction`; this
directory only. Base release v0.171.0. The frozen parent smooth-Euler,
stationary-EPS and conditional slow-affine continuum objective is unchanged.

Frozen child obligation: construct the five symmetric tracefree physical
core-deformation jets on the SAME actual EPS field used by 0045/0048/0050.
Use genuine compact volume-preserving displacements and full second-order
transported-vorticity energy, including whole-space Leray and exterior
contributions. Obtain an exact finite shear matrix and, if required, a
same-field affine-bond cage construction making it strictly positive.
Coordinate the common/translation moment projections with 0052; preserve
physical core jets and retain all fixed finite-dimensional corrections.
Derive the isotropic shear scalar and the vanishing of isotropic
symmetric-strain/rotation cross terms without assuming a desired modulus.

Candidate A: fixed affine core jets with explicit compact outer returns,
then a complete finite-dimensional momentum Schur reduction. Candidate B:
negative-helicity cages driven by a full-rank set of affine bond extensions,
using fixed off-core moment corrections. Selection criteria: exact Euler
transport, physical affine jets, finite exterior-inclusive integrals,
proved rank and jet preservation, full cross-term retention, positive
derived coefficient, and explicit falsifiable kinematics. No empirical
comparator, soft numerical eigenvalue, assigned mass, or target modulus.

Analytic oracle: curl identities, finite-dimensional moment/Schur algebra,
the actual-EPS carrier estimates in 0045/0048, and the exact isotropic
representation average. Read the 0052 joining record before selecting a
moment projection. Status: frozen before substantive construction.

Result: `affine-shear.md` constructs the exact physical affine jets, a
nine-bond full-rank cage attachment, the actual-EPS energy/exterior integral,
the rank-aware fixed-response projection, and a strictly positive complete
shear Schur matrix at an explicit finite carrier threshold. The exact
isotropic scalar is `mu=n_cell tr(C)/10`; local spin/shear coupling vanishes
by the verified representation algebra. No periodic surrogate modulus or
mass is imported. The joint positive momentum block remains an explicit
parent dependency, not a coefficient supplied by this child.

During coordination with 0052 we identified and repaired one topology
precision: the Gram range argument applies directly to our actual curls
of compact potentials inside the solid torus. Extending it to arbitrary
compact solenoidal fields would additionally require its topological flux
sector. This correction changes none of the constructed fields or jets.

`verify.py` passes 18/18 exact checks on its first execution, preserved in
`stdout.txt` and `stderr.txt`; Ruff and scoped diff checks pass. Route
verdict: established as stated, individual review pending. The parent owns
the final joint kinetic reduction, continuum joining and claim promotion.
