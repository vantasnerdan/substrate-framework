# Independent individual review — continuum action and operator repair

Reviewer: Codex `/root/construction_review`, 2026-09-05.
Author/implementer of 0056: Codex `/root`. This reviewer did not implement
`micropolar.py`, its tests, or the changed CST004/CST005 consumers.
Previously authored microscopic dependencies are not independently reviewed
again here. Scope: this fixed action/operator reconciliation, not promotion
of the unresolved material assembly.

Read the frozen README, `src/substrate_framework/micropolar.py`,
`tests/test_micropolar.py`, both changed root verifiers, and all five captured
receipts. No test rerun was necessary. The first pytest failure compares
expanded and factored SymPy expressions structurally; its recorded repair
tests the expanded difference against zero. It changes no mathematical
equation or tolerance. The repaired receipt gives 4/4 tests; CST004 and
CST005 give 8/8 and 17/17 checks, respectively.

## Independent real-field derivation

Use `h_ij=partial_j U_i`, `G_ij=partial_j Phi_i`, and
`r=curl U-2 Phi`. Direct differentiation of the stated real energy gives

    m=2 c_tr tr(G) I+(c_s+c_a) G+(c_s-c_a) G^T,
    div m=(c_s+c_a) Delta Phi
            +(2 c_tr+c_s-c_a) grad div Phi,
    partial W/partial Phi=-2 alpha r.

For a longitudinal field `Phi=(0,0,f(z))`, the curvature energy itself is
`(c_s+c_tr) (f')²`: its skew gradient is zero. Variation gives the positive
Fourier stiffness `2(c_s+c_tr) k²`. Equivalently, longitudinal projection
ADDS the Laplacian and grad-div coefficients:

    (c_s+c_a)+(2 c_tr+c_s-c_a)=2(c_s+c_tr).

The retired `2(c_a-c_tr)` expression is therefore concretely wrong; it
subtracts terms having the same Fourier sign and incorrectly leaves the
skew-curvature coefficient in a purely symmetric longitudinal gradient.

With `ax(sigma)_i=epsilon_ijk sigma_jk`, the locking stress has
`ax(sigma)=-2 alpha r=partial W/partial Phi`. Thus the spin equation is
`j Phi_tt=div m-ax(sigma)`. CST004's explicitly indexed `eps_sig` is the
NEGATIVE of this ax convention, so its displayed `div m+eps:sigma` is the
same equation, not an opposing torque convention.

For exp(i k.x), put C=i[k cross]. C is Hermitian for real k, and
`C²=|k|² I-k k^T`. The independently varied stiffness blocks are

    K_UU=(mu+alpha)|k|² I+(lambda+mu-alpha) k k^T,
    K_U_Phi=K_Phi_U=-2 alpha C,
    K_Phi_Phi=[4 alpha+(c_s+c_a)|k|²] I
                 +(2 c_tr+c_s-c_a) k k^T.

These are exactly the new API. The real cosine/sine test uses
`Phi_amplitude=-i phi`, matching Phi=phi sin(k.x), and its factor one-half
in the amplitude Hessian is correct. It checks the full real-field action,
not only a determinant re-entered with the desired coefficient.

At k along z, `(1,i h,0)/sqrt(2)`, h=±1, is a curl eigenvector of eigenvalue
h k. Each transverse block therefore is

    [[(mu+alpha) k², -2 alpha h k],
     [-2 alpha h k, (c_s+c_a) k²+4 alpha]].

The coupling product is 4 alpha² k² for BOTH helicities. The mass pencil
used in CST005 is minus `(K-omega² M)` in this two-dimensional sector;
that changes no determinant or dispersion roots. Its optical root at zero
wave number is 4 alpha/j. Longitudinal displacement is excluded by k.U=0;
longitudinal spin is not excluded and retains
`j omega²=4 alpha+2(c_s+c_tr) k²`.

## Coefficient source and singular limit

CST004 now takes the complete mixed Euler-orbit H/KKS reduction through the
shared Schur and affine-cage APIs. Its numerical positive matrix is explicitly
an algebraic probe, not a fitted EPS field or a microscopic coefficient claim.
The isotropic factors

    j=n_cell J_Psi/3,  alpha=n_cell K_Psi/12

give `4 alpha/j=K_Psi/J_Psi` exactly. Neither the old free-director tension
locking substitution nor an independently appended rigid-body inertia is
used by the changed consumers. The singular structure-free check preserves
its declared common density scaling: at positive density the spin ratio can
remain finite, while at zero density its kinetic and potential weights both
vanish and the spin coordinate is removed. The retained DOP853 refinement
is a well-conditioned operator regression, not evidence for Euler material
assembly or a source of fitted microscopic constants.

## One prose correction

At review time the CST005 introductory docstring still described both
longitudinal formulas as outside incompressible Euler. Replace it with the
precise distinction already implemented in the new module: incompressibility
restricts displacement; longitudinal spin remains an allowed mode of the
conditional micropolar system. Only the compressible longitudinal
displacement branch is a formal extension. This is an exposition correction,
not a defect in the new operator or the longitudinal-spin test.

## Individual verdicts

* **Quadratic real-field action and full Fourier operator: established as
  stated.** Variation, Hermitian coupling sign, both helicities, angular
  reaction and rigid-frame covariance agree.
* **Longitudinal-spin repair: established as stated.** The coefficient is
  exactly 2(c_s+c_tr); the exposing real-field projection rejects the old
  subtraction with a named mechanism.
* **Same-orbit coefficient plumbing and limit/regression consumers:
  established in their stated conditional scope.** The substituted inertia
  and locking terms now come from the same full reduced action, with no
  inference that the material/coherence assembly has already been promoted.

The new APIs are correctly marked conditional and unpromoted. With the
single docstring clarification above, this fixed reconciliation is ready
for its intended repository use. No further scientific narrowing or new
acceptance requirement is introduced by this review.

## Correction check receipt

The CST005 docstring now explicitly restricts only the displacement sector
and retains longitudinal spin in the conditional incompressible system.
Confirmed by reading its corrected introductory lines; no equation changed
and no scientific test rerun was needed. The requested correction is closed.

SHA-256 at this correction check:

    verify_cst005.py
    e53f1173a321509e596b421496bb015f336a51489c6999e3210485c27c991ef6
    src/substrate_framework/micropolar.py
    7187523783e3dfe30ce062c5b547e5d0524631584dcde13326646d45547a2693
    tests/test_micropolar.py
    d78808625a9af88d6c8151e6dbba4f2f4543273e7e73dfbd1547f16b25664b09
