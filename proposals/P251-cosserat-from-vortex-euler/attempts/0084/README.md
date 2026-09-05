# 0084 — constructive material transport and collar repair

Owner `/root/smooth_core_review`; authorized constructive subtask of P251 /
issue #198, not an independent review of this author's result. Write surface:
this directory only. Accepted base `v0.171.0`; campaign imports 0037 exact
Euler Jacobi action, 0078 material centroid/Kelvin distinction, and the
stationary invariant EPS tube selected in 0080. No registry delta is made.
The parent objective and its conditional slow-affine scope remain unchanged.

Positive deliverable: reconstruct an actual material tag from a velocity
perturbation; construct a divergence-free collar lift and its full-action
local spin map; state an executable condition for eliminating tag history
into angle/rate variables without replacing the actual centroid or adding
fixed-Kelvin rotor mass.

The observed failure motivating the route is precise: an instantaneous
coadjoint generator compact away from a tube boundary leaves that boundary
unchanged, but its induced Leray velocity can cross it. Equating the two
deformations loses actual material transport. Repair candidates are (A) exact
retarded material displacement for an existing Eulerian velocity sector,
retaining the resulting tag state, and (B) a material displacement/collar
construction first, whose velocity and full Jacobi action are then derived.
Both are pursued analytically. The structural selection criterion is actual
transport with complete pressure/ambient/Kelvin bookkeeping, not agreement
with a desired inertia or optical coefficient.

Sibling 0082 owns invariant-torus Fourier/retarded locality analysis; 0080
owns instantaneous response rank. This attempt supplies the complementary
material lift and action, without duplicating their constructions.

Oracle: direct material differentiation, divergence theorem, explicit curl
extension and exact symbolic checks using the existing importable
`euler_displacement` API. No numerical remainder, soft eigenvalue estimate,
new external source, empirical comparator, all-k closure or full nonlinear
invariant ansatz is assumed. A positive kinematic metric is not a proof of
positive fixed-Kelvin reduced inertia or parent closure.

Next dependencies unlocked: apply the constructed lift to a selected EPS
sector, then carry its retained circulation constraints and full reaction
operator through the same-action reduction. A fixed-sector locality claim
requires its displayed compatibility identity, not a new parent criterion.

## Completed receipt

`material-lift.md` constructs the retarded material displacement and actual
tag, a divergence-free collar with any volume-compatible normal trace, and
the local full-spin formula `delta S=A zdot+C z`. Its velocity contains the
required commutator, and the full Jacobi mass/gyro/stiffness and Kelvin
one-form are retained. The surface transport operator includes the actual
surface divergence; no integrable collar foliation is presumed.

The exact Beltrami stiffness identity is
`K/rho=lambda integral F.curl F-integral|curl F|²`, `F=Xi cross u`.
An explicit smooth single-wave background exposes negative compact material
stiffness; the same identity supplies the distinct positive same-helicity
low-curl-band candidate pursued by main in 0088. Neither result borrows the
old orbit cage's sign or assigns fixed-Kelvin rotor mass.

`route_verdict: established`; `evidence_scope: REPRESENTATION_SCOPED` for the
construction and its stated conditions. The old instantaneous EPS sector's
complete constrained material reduction is a remaining parent construction,
not an established consequence of this receipt.

`verify.py` imports the existing displacement/Jacobi API and passes 15/15
exact checks. The preserved first run exposed an incorrectly signed
solid-rotation test pressure; the fixture was repaired and a stationary
Euler residual added before testing the unchanged identity. `first-run.txt`
and `stdout.txt` retain both executions. No numerical stability inference
or new canonical symbol was introduced. Direct consumers: 0080/0082 material
joining and main's 0088 sign construction; independent review is requested
at this frozen boundary.
