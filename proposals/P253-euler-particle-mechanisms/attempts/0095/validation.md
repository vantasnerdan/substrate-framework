# Validation and claim boundary

The analytic derivation is load bearing.  The exact verifier and focused API
tests expose algebraic signs, normalizations, and paired-mode structure; they
do not replace the compactness, pseudodifferential, or nonlinear path proofs.

The Recovery-7 Sol-High verifier checks fourteen exact predicates:

1. the comoving transverse-Maxwell characteristic polynomial;
2. the column BAS square and Rayleigh discriminant;
3. the positive column metric skew identity;
4. the exact streamfunction-normal identity return;
5. the universal metric-only first-curvature Melnikov row;
6. the Cao odd-cell affine solvability moment;
7. the full `r^2/2` speed polynomial's zero BAS strain;
8. the HF trace-free-system-to-Hill elimination;
9. the HF smooth-center forcing coefficient `15 Omega_c/128`;
10. the physical HF growth bracket `15/256`;
11. both reciprocal uncharged Hill couplings and their hyperbolic discriminant;
12. the exact-return detuning polynomial `lambda^2+d^2-b^2`;
13. the charged-column slow-block square `-x^2 Z_tau R_tau I`; and
14. the two-sign constant-lambda helicity symbol.

The focused public API tests independently exercise the HF elimination and
center coefficient, the detuned Hill tongue, and the charged-column block.
Both use the repository interpreter and import
`substrate_framework.euler_p2_principal`; the attempt verifier does not import
the compatibility wrapper.

Analytic checks outside the script establish the exact phase identity
`(partial_t+W dot grad)(F(P)+ell theta)=0`, return of its differential, the
axisymmetric-row Fourier orthogonality on the smooth compact packet core, the
exact raw coadjoint/tag/Gauss path, the `A_2` weight interval
`|sigma|<3/2`, the fixed
observation-operator quantifiers in (36bf)--(36bg), including the essential
norm of the fixed observed propagator `J_loc S_g` from the weighted global
complete-state input to the local same-`H^s` output, and the conditional
fixed-time same-target differentiation bridge once an exact fixed-row path
exists.  No continuous finite-row submersion on the completed weighted space
is claimed.  The source reconstruction
uses the cached HF PDF at SHA-256
`c6a35c44d55d26b8bb3e8fc9e55b29b638d0983b646b17eb522a4192b6444fb9`.

No production numerics are used.  The verifier does not prove a uniform
charged Cao-to-column branch, a delta-uniform monodromy constant, exact
charged reciprocity, a relative-total-momentum domain, the carrier-specific
finite-row witness determinant, the combined nonlinear Lipschitz refutation,
a full Hill hessian, a one-datum nonlinear escape theorem, or LP2.
