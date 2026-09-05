# 0048 — direct EPS relative-orbit common rotor

Parent: P251 / issue #198, same smooth-Euler conditional affine objective.
Owner: Codex `/root/construction_review`; this directory only.

Frozen positive construction: use the actual R³ EPS stationary Beltrami
field from 1210.6271 and its invariant vortex tube, not the finite cylinder
as a surrogate. Enlarge its compact isovortical perturbations by one actual
global Euclidean rotation about a chosen axis through the core point.
Derive a finite relative quadratic energy and KKS action, a physical common
angle and positive conjugate response, and combine these with a compact
physical core-angle sector on the SAME Euler orbit. Obtain both inertias by
eliminating the same action's conjugates; supply no body mass or desired
modulus. Preserve source reaction and physical angle observability.

Candidates: (A) renormalized energy over rotation-invariant balls, exploiting
the exact global symmetry and its zero pointwise Beltrami Hessian; enforce
the finite KKS cross constraints on the compact generators and retain the
full positive momentum block. (B) if that relative action fails, a compact
rotation with an explicitly computed interface/ambient reaction. Structural
selection is Euler/Kelvin consistency, finite well-defined coefficients,
nonzero physical moment map, full positive Schur complement, and actual EPS
ambient compatibility. No numerical spectrum or empirical comparator.

Analytic oracle: exact symmetry and integration identities, source-supported
far-field bounds, finite-dimensional symplectic linear algebra, and the
explicit finite-wave-number estimates of 0045. The relative orbit and
renormalization are part of the declared construction, not implicit
permission to subtract an arbitrary counterterm. Verify source asymptotics
before relying on them. The parent spatial homogenization remains distinct.

Status: active. Contract frozen before the derivation.

Append-only candidate expansion: the full four-coordinate reduction retains
an odd-in-circulation gyroscopic term. Use the time-reversal-paired ensemble
{u0,-u0} with equal weights and coherently tied physical B/q displacements.
Both realizations have identical EPS tube geometry and pressure. The
zero-handedness/equal-circulation-sign premise is independently falsifiable;
it selects neither a modulus nor a target dispersion. Derive the action for
each realization first, then average: H is unchanged and Omega changes sign,
so even kinetic/restoring terms remain and odd gyroscopic terms cancel.

The conjugate fluid shapes are independent across the two signs; only B,q
are coherently tied. This is essential to retaining the physical inertia.
The finished construction is `direct-rotor.md`, including the exact relative
orbit, EPS far-field source proof, fixed-eta0 moment projection, full finite
positive Hessian estimate, independent-momentum paired action, and general
physical field map. `verify.py` gives 20/20 exact checks in `stdout.txt`;
`first-run.txt` preserves the earlier 18-check execution before the two
explicit ensemble-order checks were added. Ruff and scoped diff checks pass.

Route verdict: established as stated, individual review pending. The parent
continues spatial affine/gradient joining on this same EPS action.
