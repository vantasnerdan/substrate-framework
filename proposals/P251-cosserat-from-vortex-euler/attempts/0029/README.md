# Attempt 0029: Euler correspondence and physical field types

Review boundary: PR #199 at 10a0f31, incorporating de980dc and attempts
0026–0028. Original P251 objective and obligations remain unchanged.
Method repair: derive the Euler operator in Cartesian components before
comparing any historical formula, then track dimensions and variational
coordinates before matching coefficients. Exact algebra, no numeric remainder.

The background v0=Omega*(-y,x,0) gives the polar acceleration components

    -i*wt*vr - 2*Omega*vtheta,
    -i*wt*vtheta + 2*Omega*vr,
    -i*wt*vz.

The original verify_cst002 velocity formulas solve these equations exactly.
Attempt 0019's alternative residual operator reverses both Coriolis signs.
Its supposed field correction therefore solves the wrong operator. This
corrects the previous review's reliance on 0019; it restores the original
field evidence rather than rejecting that supported part of N2. The energy
receipts depending on 0019's altered fields cannot supply the modulus/action
bridge. Attempt 0031 continues by repairing the separate Bessel derivative
and its subleading bending constant.

Attempt 0028's proposed alpha_E has dimensions mass/length (microinertia).
N4's angle stiffness has pressure dimensions, mass/(length*time^2). Equality
of those coefficients, or of their purported coincidence condition, is
dimensionally invalid. The Euler–Lagrange variation of J*q_dot^2/2 produces
J*q_ddot, whereas K*q^2/2 produces K*q. Integrating a rate in prose changes
neither the functional nor its units. The revised verify_cst003 preserves
exact conditional moment and angle-energy matching while removing the
unsupported map and the old truncated-rotation energy from its general probe.

The finite-k composition also needs repair: transforming to an ambient
rotating frame leaves the Coriolis term in three-dimensional Euler. A bulk
axial inertial-wave probe at the proposed contrast-only frequency yields
the nonzero determinant 4*Omega_o*(2*Omega_i-Omega_o). This refutes the
general replacement of rotation by contrast at finite k; the planar k=0
contour-advection identity is retained.

`route_verdict: established` for these correspondence identities and explicit
counterexamples. `evidence_scope: EXACT_OPERATOR_AND_FIELD_TYPE_AUDIT`.
N2/N3 stay active. The next executed constructions are 0030 (elliptic-patch
angle action), 0031 (corrected mode equation), 0032 (closed collective vortex
action), and 0033 (repaired mutual kinetic-energy kernel).

Verification: verify_correspondence.py exits zero, 10/10 checks, captured in
stdout.txt with empty stderr.txt. The repaired verify_cst003 exits zero,
23/23 checks, captured in cst003.stdout.txt and cst003.stderr.txt. These are
scientific predicates, not a campaign-completion tally.

No accepted claim or existing canonical API changes. GitNexus's index targets
the separate main worktree and lacks P251 symbols; its impact lookup reports
unknown and its compare report concerns P250, not this transaction. Direct
source search finds the changed verifier functions called only by their own
main entrypoints. New rankine_modes consumers are the N2 verifier, 0031
receipt and tests/test_rankine_modes.py. This bounds replay to those paths;
N1/N4–N7 retain their previously validated conditional evidence.
