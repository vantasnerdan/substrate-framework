# 0093 — material gradient cages with their full kinetic normalization

Parent P251 / issue #198; owner `/root/construction_review`; this directory
only. Frozen positive task: using the ACTUAL material Jacobi action and
0091's physical-spin-equals-mass base, construct gradient-only attachments
whose stiffness gain remains positive after their added gradient inertia
and all retained mean/field-map and Schur corrections are included.

Inputs to verify at source: compact positive material cage X from 0090,
rotational profile B_L=chi(r/L)(e cross r), A(B_L)-M(B_L)=d L^5>0,
and spin/mass completion Xi=B_L+tX with
M_X t²-A_X t=d L^5. No modulus or target frequency is fitted: t enforces a
computed physical spin/kinetic identity and the cage ratio is measured from
the same Euler material Hessian and kinetic norm.

Candidate mechanism: the base optical ratio tends
(K_X/M_X)*(d/a), with 0<d/a<1, so the SAME witness cage has a strictly larger
stiffness/mass ratio than the base gap. Attach disjoint finite-difference
gradient cages and retain BOTH stiffness and kinetic jets. Positive residual
requires K_grad-(kappa/j)M_grad>0; K_grad>0 alone is insufficient.

Oracle: exact scaling, finite inequalities, physical moment/current map,
full quadratic-form/Schur algebra and isotropic tensor contractions. No
empirical comparator, soft eigenvalue numerics, all-wavelength invariance,
isolated unlicensed energy sum or omitted kinetic correction. Source/frame
and boundary assumptions remain those of the actual material construction.

Completed route: `material-gradient.md` proves a finite same-cage optical
ratio margin by the exact identity
`(K_X/M_X)j-kappa=(K_X/M_X)M_B-K_B`. Paired translated gradient cages have
zero zeroth AND first material moments, so their coherent mean starts at
order k³ after the gradient attachment. The proof retains their added
gradient inertia and the full physical/normal field map, giving a strictly
positive gain `K_grad-(kappa/j)M_grad` in both curvature eigenvalues.
Existing finite reaction profiles are kept support-disjoint; their complete
Schur reduction commutes with this new jet. This is the declared material
Cauchy--Born construction, not freely relaxed nonaffine Euler shapes.

`verify.py` first execution is preserved in `stdout.txt`: 23/23 exact checks,
exit zero. Mutations expose omitted gradient mass and an overlapping free
reaction falsely treated as support-orthogonal. No soft numerical oracle
or empirical comparator is involved. Established conditional gradient
construction;0091's actual mean/core/Kelvin action license and0094's actual
shear are separate parent inputs, not claimed to follow from this proof.
