# Independent scientific review: C-CST-012

Reviewer: /root,2026-09-05. Author of reviewed0161: construction_review.
Base v0.176.0/dbf0c04. This is one substantive review of the new
triangular-array theorem, not another review of accepted C-CST-008..011.
Reviewer authored earlier unchanged finite-phase algebra and registered
the candidate; reviewer did not derive0161's group proof or write its
verifier. Those prior roles are disclosed, not counted as independent
evidence for the new theorem.

## Claim and positive role

For Psi,lambda,rho>0, let psi=Psi sum_j cos(b_j.x), with three planar
wavevectors of length lambda summing to zero, on their triangular
periodic cell, and u=(J grad psi,0). The full zero-axial-wavevector
horizontal linear Euler group is uniformly bounded in L2 vorticity
on the real vector isotypic sector under sixty-degree rotations,
including its translation kernel. For axial Bloch k approaching zero,
the actual common-velocity and Kelvin-prepared displacement Euler
histories obey0161(11)-(14) on every fixed slow interval t<=T/|k|,
with O_T(|k|/lambda) normalized error and C_v=3 Psi² lambda² I/4.
The exact initial material phase pairing has mass rho and initial
displacement stiffness rho k² lambda² C_v/(lambda²+k²); the complete
moving action in the compensated-current slow chart has the stated
leading mass and stiffness with controlled normalized corrections.

X here is explicitly the time integral of the full Euler velocity mean,
not an unproved identification with a marked material centroid. The
compensated current p_c differs from that mean and is retained before
differentiating. Common-V and prepared-D have their specified different
initial circulation variations. Nonlinear amplitude is taken last on
each finite interval; no k-uniform nonlinear amplitude is asserted.

This supplies actual smooth-fluid acoustic-time dynamics rather than
a prescribed oscillator, fixed-time Taylor coefficient or assumed
full-cell stability. It advances the parent translational obligation.
The constant-curl lift, generic wavevectors, isotropy, EPS topology and
same-field optical closure are not in this claim and remain active.

## Frozen transaction and evidence

The reviewed primary proof is0161/planar-acoustic-array.md, SHA256
8238ac33e81e3750420234a8dc0e74ecba3b6af3a0fa6baf074e244eab6bcede.
Its verifier SHA256 is
a5a8ea8f7325e2dd80cb2abe3e68076a8a71892a105b57baca1c484b769e19f5;
the successful output is
cb5af12d8f22dee411357102ba3902f9dc4c8cd298d0d56b1c9eb71d9325630d.
0146/acoustic-normal-form.md sections1-3 supply the full-pressure
normal-form proof, read directly at source. These are new analytic
evidence under the standard incompressible Euler/Lin, periodic Fourier,
transport-group and finite-time smooth-solution imports, not authority
borrowed from an unaccepted campaign number.

The primary oracle is the complete analytic proof. The20 recorded exact
checks support its field, oblique lattice, boundary, representation,
forcing and pressure/action normalization. They are not20 independent
proofs of the PDE bound. In particular the isolated arithmetic check
1-1/3=2/3 is only a consistency anchor: the all-mode gap follows from
the positive integer quadratic form and its modulo-three exclusion,
not that literal check. No numerical spectrum supplies the result.
The original implementation failures and their repairs are preserved;
their corrected coefficients and source match the proof. Reuse their
unchanged passing execution rather than regenerate the same receipt.

## Oracle audit and independent reasoning

The reciprocal quadratic form n1²+n2²-n1n2 is positive for every nonzero
integer pair. Modulo three it equals (n1+n2)², excluding2; hence the
complement of its actual first shell has H>=2/3. Rotation by pi removes
the cosine shell in the excited vector sector. The actual three-sine
rotation matrix splits into a two-dimensional vector and a distinct
eigenvalue-minus-one component, leaving precisely the two translations.
No reflection symmetry of a signed vortex is needed.

The separatrix factorization partitions the universal cover into bounded
invariant polygons. Their centered coordinates are periodic and belong
to the transport domain because normal transport flux is zero on every
edge. Integration of partial_i(r_j psi) on each polygon retains the
common boundary level; after summing, the matrix is
(psi_s-mean psi)delta_ij=-Psi delta_ij. This also shows streamfunction
gauge invariance: setting psi_s to zero while forgetting the changed
mean would be an invalid mutation. The paired adjoint rows satisfy
A h_i=t_i and have a nonzero antisymmetric matrix on the kernel.
Their conserved values control both kernel coefficients after the
positive-complement energy controls f_perp. This is the essential
step that prevents secular kernel growth; energy positivity alone
would not prove the claim.

The C0 realization follows from unitary smooth transport plus bounded
lambda² A(-Delta)^-1. The adjoint rows have only permitted transport
jumps, so density extends the energy/row identities to the claimed
space. The independently expanded curl of F has only smooth zeta,
grad zeta, grad(-Delta)^-1 f and f: no derivative-loss assumption is
hidden in the Y mapping. The complete pressure maps are bounded in
the actual periodic mean-zero spaces and preserve the vector sector.

I checked0146's slow scaling directly. In original time each coupling
of y,Z,x=|k|X,m outside the bounded fast groups is O(|k|); pi_r is
bounded by y+|k|Z+|k|x. Gronwall therefore has exponent C T, not C T/|k|.
The compensated-current equations then have uniformly O(|k|) slow
remainders. Duhamel gives the actual physical mean estimate. In the
prepared displacement data, retaining both pressure-return components
gives the required small y0,Z0 and the two-column C1 slow-chart estimate.

The material pairing follows from eta_D(0)=X0 and
pi_V(0)=rho V0, while the prepared-D return has zero horizontal mean.
Its full axial Leray projection has squared norm lambda²/(lambda²+k²),
giving the exact initial stiffness, not the unprojected k² coefficient.
The pressure Hessian averages to zero on the same periodic cell.
At later times the slow chart uses p_c, whose derivative is controlled;
the proof explicitly avoids differentiating an O(k) estimate for m.
Thus the exact moving-action identity legitimately supplies the leading
slow action. It does not turn the initial phase Hamiltonian into an
exact autonomous finite-k Hamiltonian.

## Findings and compatibility

No load-bearing defect found within this claim. No scope reduction or
correction requested. The finite action, density and pressure conventions
are consistent; all inverse operators and conserved quantities come
from the same fluid cell. Declared circulation data are genuine initial
data, not an extra mechanical mass. The original compact-core candidate
is neither disproved nor replaced by this alternative smooth array.

The next construction consumers are0163 and0166; neither may inherit
the planar theorem without deriving its new axial terms. Canonical
extraction and tests remain part of the promotion transaction, with
targeted downstream validation under the user's explicit scope choice.
These implementation tasks do not reopen the scientific proof.

## Four-axis decision and promotion

Verification: symbolic_verified for the exact algebra, with the analytic
PDE/group proof explicitly recorded as the primary evidence.
Review: accepted.
Compatibility: compatible_extension.
Epistemic: active.
Relationship: additive actual smooth-fluid acoustic theorem.

Promote only the statement above with its exact domains and errors.
Materialize reusable field/projection definitions and their tests,
immutable adjudication, registry, pinned release, generated views and
memory. Reuse unchanged receipts and validate the actual delta; no full
suite is requested. Correction check: not needed. Parent objective
remains active;0163/0164/0166 are its immediate constructive work.
