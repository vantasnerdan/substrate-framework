# 0078 — finite material Euler parcel: centroid, frame and Kelvin reduction

Parent P251 / issue #198; owner `/root/construction_review`; this directory
only. Frozen task: derive the exact finite-parcel mass-centroid and intrinsic
angular momentum decomposition from the full material Euler action, with an
explicit equivariant frame/shape chart and retained circulation constraints.
Establish when a positive absolute-spin kinetic term survives fixed-Kelvin
reduction, distinguishing geometric locked inertia from an Euler-orbit
effective inertia. No rotor mass or desired Cosserat law is supplied.

Fixed theorem routes: exact material kinetic/cotangent decomposition; polar
or Eckart frame chart; mechanical connection and Schur complement; commuting
spatial/relabeling momentum maps and a positive physical-angle orbit sector.
Sources inside the campaign: 0037 material Jacobi action, 0042 Kelvin mismatch,
0052 finite centroid/currents and 0073 parcel multipole map. This is a general
action license, not 0075's constructive partition or compact-velocity work.

Oracle: exact variations and constraint reduction, finite-dimensional block
identities where they expose the infinite-dimensional theorem, and explicit
frame/gauge counterexamples. No numerical eigenvalue, empirical comparator,
fitted coefficient, new all-k gate or parent completion claim.

## Completed local receipt

`finite-parcel-reduction.md` gives the exact centroid/cotangent split,
equivariant polar frame, mechanical connection and shape-momentum Schur
formula. It distinguishes a label-dependent geometric locked inertia from
a physical Eulerian vorticity-angle observable on a fixed Kelvin orbit.
Commuting spatial/relabeling actions preserve the physical angular momentum.

A positive constructive fixed-Kelvin theorem is supplied: a physical rotation
angle and same-leaf momentum-changing direction with B=dJ!=0, h=H_etaeta>0
give I_red=B²/h. A smooth Beltrami field in a finite spherical slip/material
domain realizes the hypotheses; compact same-leaf responses and finite
negative-helicity cages supply the positive h. This example explicitly
retains its spherical shape/pressure boundary constraint. It is NOT silently
identified with the parent's EPS/ambient partition, whose compatibility is
0075's separate construction.

Route verdict: established as stated. Evidence scope: exact finite-parcel
action and positive fixed-Kelvin spin license, not parent Cosserat closure.
The locked geometric inertia is never added to an already reduced orbit
inertia. All new coefficients are calculated mass moments or Euler KKS/
Hessian expressions, with no fitted parameter or numerical stability gate.

Verification: `verify.py` imports the canonical Schur-jet and ledger helpers;
17/17 exact checks pass, exit zero, 1.866 seconds. The first execution caught
a checker implementation error: substituting z² in an expanded radial
constraint left its z³ term unreduced. `first-run.txt` preserves that abort;
polynomial remainder modulo x²+y²+z²-r² repairs the checker without changing
the field equation. `stdout.txt` preserves the repaired execution. Ruff and
scoped diff whitespace checks pass.
