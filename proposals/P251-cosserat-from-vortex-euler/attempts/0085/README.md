# 0085 — compact induced velocity in a genuinely varying Beltrami field

Owner `/root/orientation_construction`, this directory only. Continuation
of 0080's physical-tag reconstruction obligation. Parent objective and
permitted slow-affine material closure remain unchanged.

Frozen candidate: construct compact divergence-free xi in the actual
smooth EPS tube with div(xi cross omega)=0, so the full Leray induced
velocity is itself compact and the material boundary has no omitted
normal motion. Require a nonzero actual angular response and retain the
physical core-angle/KKS construction. The constant-vorticity no-spin
result in 0075 motivates using genuine spatial Beltrami variation, not
reusing that local jet as an existence proof.

Candidates: A exact variable-coefficient potential/flow-box construction;
B an exact collar correction cancelling the harmonic pressure trace,
with all corrections and norms retained. Selection is by exact compact
support, Euler/Kelvin admissibility, nonzero physical spin, and analytic
control of corrections. Finite moment prescriptions alone do not solve
the boundary condition. No soft numerical oracle is designed.

Result: `compact-spin-theorem.md` gives an exact order-six differential
operator construction of compact induced velocities with all three
physical spin components on an actual EPS field. Its exposing finite
jet ranks are235/236/238 over both a finite field and exact rationals;
the analytic universal rank upper bound supplies the proof bridge.
`positive-compact-pair.md` then retains the full quadrature return to
prove positive compact H and nonzero KKS, attaches an actual core angle,
and matches physical spin by disjoint compact interior responses.
Both coordinate and reaction velocities have zero material boundary
normal trace. Complete time-reversal pairing retains the physical
current after momentum elimination.

The scripts and first outputs distinguish the initially incomplete
modular certificate from its analytic rank repair and independent exact
integer replay. `wkb_pair_verify.py` is a separate exact nine-check
oracle of the corrected symbol, not a soft eigenvalue computation.
The parent retains continuum gradient, macro action and governance
obligations. No nonlinear unrestricted invariant-manifold claim is made.
