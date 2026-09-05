# 0096 — compact coadjoint shear and gradients with the full reaction

Parent P251 / issue198; owner `/root/construction_review`, this directory
only. Frozen constructive continuation of0085's positive compact
angle/reaction sector. Use the actual same-isovortical Euler H and KKS,
with BOTH generator and induced velocity compact inside the invariant EPS
tube. Construct macro STF and spin-gradient attachments, eliminate the
complete independent fluid momenta, and retain their full gradient mass,
mean identification and physical current. No coefficient from0093's
different material-Jacobi functional is imported.

Candidate A: zero-mechanical-spin compact quadrature pairs, with a large
finite carrier giving a strict reduced stiffness/inertia ratio margin;
attach their coordinates to STF strain or neighboring relative angles.
Candidate B: paired nonzero-spin cells with explicitly matched opposite
current amplitudes. Prefer A if its full H/KKS and finite mean moments
close, because its zero current is an exact Euler moment condition.

Oracle:0085's finite analytic operator/WKB bounds, full block Routh/Dirac
algebra, material centroid cotangent identity and exact moment expansion.
Mean mass comes only from the joint Euler material action, not from an
added isolated rotor or a whole-space Galilean relative-energy integral.
Original finite conditional Cauchy--Born geometry remains permitted;
unrestricted all-wave-number invariant manifolds are not the target.
No empirical comparator, fitted modulus or soft eigenvalue solver.

Completed route: `compact-action.md` constructs compact zero-spin canonical
cages directly from0085 and proves their full reduced ratio
`K_g/J_g=det(H_g)/B_g²` grows at least as a positive constant times the
carrier squared. This supplies a finite strict gradient stiffness/inertia
margin. Exact zero zeroth/first velocity moments remove the added mean
return through degree two. The full reaction square, including macro
forcing, yields positive STF shear and retains its gradient inertia.
Both normalized spin curvatures can then be made positive at finite
attachment strength; no material-Jacobi coefficient is transferred.

`verify.py` first execution is preserved in `stdout.txt`:22/22 exact checks,
exit zero,0.554321772 seconds. Ruff and scoped diff checks pass. Mutations
expose averaged inverses, omitted reaction crosses and omitted gradient
mass. The exact hybrid mean/current identification is part of this route;
0097's shared material/Kelvin affine lift remains the explicit joining
input, not a consequence of an isolated coadjoint pair. No parent claim is
promoted by this attempt.
