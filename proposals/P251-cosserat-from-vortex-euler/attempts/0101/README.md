# 0101 — a Kelvin-prepared optical tilt packet with controlled Euler memory

Owner `/root/smooth_core_review`; P251 / issue #198, this directory only.
The parent conditional slow-affine objective and accepted base v0.171.0
are unchanged. This is a new constructive concept after 0095, not a review
of that complement theorem or an exhaustion claim.

Positive deliverable: an actual isovortical initial family whose unrestricted
linearized Euler evolution retains a finite-frequency physical core-angle
signal through an optical period, with an analytic remainder rather than
zero initial complement or an arbitrary positive Dirac pair.

Registered candidates: (A) the exact near-axis rotation symmetry of a
smooth Lundquist Beltrami vortex and a compact paraxial tilt packet;
(B) the nearly field-parallel compact WKB sector suggested by main, with
careful separation of Eulerian Doppler frequency from material rotation;
(C) symmetry/Jordan-chain closure. The ordinary exact rotation symmetry
alone gives zero-frequency directions, so it does not supply the requested
optical signal. Candidate A is executed first: it supplies an intrinsic
frequency in an axial Galilean frame, a direct PDE residual, and exact
Kelvin preparation. No MHD stability or dynamics theorem is imported.

Inputs: 0037/0084 Lin and Jacobi identities; 0095 exact residual/propagator
logic; smooth EPS existence and decay already pinned in the campaign;
elementary Bessel identities verified from their differential equations.
The older singular Rankine/Kelvin-wave route is related prior construction,
but its sheet matching and long-wave coefficients are not imported into
this smooth finite-core packet.

Oracle: exact vector calculus and separated norm bounds, then the full
linearized Euler L² energy estimate. All carrier/support choices are made
from these analytic bounds. No solver, eigenvalue discretization, empirical
comparator or observed-order numerical inference is needed.

Geometric boundary: adding a small same-eigenvalue Lundquist field to an
EPS seed preserves its knotted invariant domain. The new packet occupies
a distant near-Lundquist finite straight flow-box on that SAME stationary
field; its physical tilt is not claimed to be the angle of the distant
knotted torus. This distinction is retained while assessing the parent
consumer.

Result: `optical-packet.md` establishes the scoped finite-time theorem.
`verify.py` derives the compact curl, full residual, Kelvin preparation,
Bessel orders, physical tilt, and parameter/error scalings: 20/20 checks
on the first execution, archived without a repaired scientific rerun in
`first-run.txt`. Analytic L² evolution and actual transported positive
weight estimates supply the non-discretized PDE remainder. The leading
complete packet has zero angular impulse; its local tilt does not claim
a nonzero total-spin rotor or parent completion.
