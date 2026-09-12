# Interrupted author state — user-requested pause, 2026-09-07

This file preserves work and conversation-level findings without completing
or adjudicating the interrupted attempt. The activated README is SHA-256
`59cfe8d5ec379fdded8f34ecdf951ae9a66c6b910603f8dbef1fc1e94f91f345`.
Its original activation and corrected momentum-domain replay quartets remain
append-only. `pause-manifest.sha256` identifies the scientific/activation
bytes at the pause, before recovery documentation was added.

## Route A: exact relative-momentum leaf

The draft derives forced hydrodynamic impulse and Maxwell momentum exchange:

    d I_h/dt = integral f,   d P_EM/dt = -integral f,
    J_ren = Delta I_h + Delta P_EM.

The earlier velocity-integral momentum domain was rejected: divergence-free
L1 velocity has zero integral, while Lorentz exchange can generate a nonzero
impulse/r^-3 tail. Hydrodynamic impulse and tag center are auxiliary initial
slices, not additional conserved coupled rows.

For a nonzero localized divergence-free L2 field F, the draft proves
injectivity of y -> P_L(F cross y) and builds compact smooth witnesses using
vector-potential cutoffs followed by mollification. Applied to omega_g and
B_g, this gives impulse and electromagnetic matrices H and G. The actual
joint row derivative is lower triangular with K_I and K_C retained. A closed
tagged streamline circulation supplies the carrier-specific B_g != 0 argument.
High toroidal character cancels the physical linear vector rows. A material
push-forward plus electric transverse correction is proposed to give the
exact fixed-J_ren initial curve.

The local Maxwell oracle's first failure is preserved. The corrected identity
exposes the unconstrained defect `-(B/mu_EM) div B` and recovers the physical
momentum equation on `B=curl A`. The recorded corrected run reports ten
algebra checks and exit 0; it does not prove the witness/IFT/PDE bridges.
The verifier was being edited during the interrupted author pass. Its final
pause hash and historical run receipt must be matched on resumption before
calling that run validation of the final source.

Known open integration point: Section 9 calls the nonlinear obstruction
established while the route ledger retains a prospective status. 0104 supports
a conditional same-target pathwise differentiation implication; it does not
independently establish all antecedents for the newly corrected curve. Finish
the fixed-(j,N) finite-time energy argument and check carrier regularity,
exact rows, and the input/output norm before choosing the final verdict.
Order of limits is fixed j, fixed N, amplitude -> 0, then N -> infinity.
No uniform N/j differentiability or full operator-norm consequence is needed.

## Route B: distinguish angle and compact-edge limits

For C=x[[0,-Z],[R,0]], H0=diag(R,Z) is a positive order-zero symmetrizer
independent of x whenever R,Z have common positive bounds. The KKS Hessian
|x|H0 degenerates as x -> 0, which prevents that energy from directly
controlling the physical norm; this alone is not linear instability.

At the compact edge Z -> 0 with R -> R_edge > 0, the frozen-column propagator
has a lower entry of magnitude sqrt(R/Z) at
t=pi/(2|x|sqrt(RZ)); at Z=0 it is a nilpotent shear. The finite-Cao accessible
packet/interface transfer is still open, so this is a column mechanism, not
a completed global Cao conclusion. The separate reviewed center hyperbolic
tongue remains valid at its fixed-delta small-charge scope.

## Route C: recover the unwritten Hill result

p4 was stopped before creating `hill-route-c-hessian.md`. The following is a
record of its preliminary author finding and coordinator/supervisor algebra,
not a substituted final proof or independent review verdict.

For a compact-interior, divergence-free displacement xi of an uncharged
relative steady Euler state, the proposed exact material Hessian is

    Q/rho = ||P_L(xi cross omega)||_2^2
            + integral W dot [xi cross curl(xi cross omega)].

Use a constant unit covector k and circular polarizations a_sigma on k-perp.
An exact divergence-free compact packet is obtained as a curl of a compact
oscillatory vector potential divided by N. The first term is O(1). The
second has opposite leading signs proportional to
N (k dot omega)(W dot k), with the real-packet factor fixed by the chosen
polarization normalization (p4 reported N/2 times the cutoff-square integral).
Cutoff corrections are O(1), and compact interior support eliminates the
interface displacement. The whole-space Hodge operator retains exterior
velocity; it is not replaced by a boundary truncation.

For Hill radius a, relative vorticity lambda, and c=2 lambda a^2/15, at
(s,z)=(a/2,0): W_z=3c/4, omega_theta=lambda a/2. With the covector fixed at
that point in direction (e_z+e_theta)/sqrt(2), the product is
3 lambda c a/16 != 0 and retains its sign on a small cell. The reported
two-by-two Hermitian symbol has eigenvalues of opposite sign. On the reported
first-order DA graph normalization, amplitude 1/N gives Q of opposite
O(1/N) signs. Norm order, exact packet factors, translation/constraint
corrections, and the full written proof still need pinning.

This is promising evidence against one-sign full-3D Hill orbit-Hessian
coercivity. It does not establish spectral or nonlinear instability or close
P2. The symmetry axis has omega=0 and the stagnation circle has W=0, so neither
is the chosen exposing cell.

A proposed generalization uses any smooth steady cell with W and omega both
nonzero, choosing k outside their two orthogonal planes. This remains a
separate candidate theorem requiring its own contract/domain proof. The
suggested next strategy, if that theorem is earned, is a dynamical or invariant
graph-metric route rather than another local one-sign Hessian search. No new
attempt has been activated for that generalization during this pause.

## Resume and evidence boundary

Preserve draft bytes and historical failures. Read `../../PAUSED.md` and
0104/0106 final reviews before resuming. No result.yaml, final validation,
completion manifest, or independent 0107 review existed at interruption.
Author claims in derivation.md are therefore working statements. This pause
does not promote them or change any earlier review.
