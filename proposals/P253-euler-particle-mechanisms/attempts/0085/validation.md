# Validation

## Exact algebra

`tests/test_euler_charged_hessian.py` and the attempt verifier derive:

1. equality of the physical `E_EM-c P_EM,z` density and its completed-square
   form;
2. the parallel and transverse limits of the anisotropic charge Schur
   multiplier;
3. the multiplier by explicit minimization over the sourced transverse field
   pair;
4. both regular-band material tag-locking coefficients;
5. the translating-frame Maxwell wave denominator and its steady subluminal
   limit;
6. the nonzero-frequency radial shell root, radial derivative, and coarea
   weight; and
7. invalid-domain rejection.

The early focused runs exposed three structural-equality assertions against
factored SymPy forms.  All failures are preserved.  The corrected assertions
use derived algebraic equality; the final focused receipt has seven passes.

## Analytic evidence

The analytic derivation, not the algebra tests, supplies:

- the exact instantaneous constraint interpretation;
- the two-sided homogeneous `H^-1` bound;
- the dynamically accessible high-frequency Cao packets and compact Green
  estimate producing strict negative fluid directions;
- the pure-field positive directions and infinite-index saddle conclusion;
- the scaled compact-curl cutoff plane-wave sequences at nonzero frequency and
  the joint `k_L->0` construction at zero, including graph domain, weak
  nullness, localized/tail decoupling, finite-row removal, and the left-
  Fredholm-regularizer contradiction; and
- the distinction between an embedded eigenmode, a dark source, and an
  outgoing resonance.

The charged-branch transfer uses the independently accepted P253/0080/0084
finite-window branch.  No test proves the limiting-absorption theorem, the on-shell
current value, spectral stability, nonlinear persistence, or particle
interpretation.
