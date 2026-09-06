# Validation

## Exact oracle

The repository interpreter executed:

    python proposals/P253-euler-particle-mechanisms/attempts/0058/verify_nonlinear_range.py

The exact command, standard output, empty standard error, and exit `0` are
materialized as `nonlinear-range.command.txt`, `nonlinear-range.stdout`,
`nonlinear-range.stderr`, and `nonlinear-range.exit`.  This was an exact SymPy
run, not production numerics.

The thirteen assertions expose:

1. both input vorticity amplitudes are transverse;
2. the output numerator is transverse;
3. both Biot--Savart velocity amplitudes are computed rather than prescribed;
4. the quadratic numerator is exactly `(1,-1,0)`;
5. the selected output transport divisor is exactly zero;
6. the numerator at that zero is nonzero;
7. the DA triple-product symbol is derived exactly;
8. the positive-core nonzero-harmonic centralizer determinant is
   `i ell^3 zeta^3`, hence nonzero;
9. the Hodge cross-product square on the transverse plane has the correct
   sign; and
10. the physical rotation convention gives
    `J'(a)=a ell sigma` and `j_2=ell sigma/2`.

Several lines above contain more than one primitive assertion; the verifier's
actual count is thirteen and is printed by the captured run.

## Strongest claims validated

- The Euler quadratic tensor is not universally divisible by the raw output
  transport factor.  The counterexample is exact and sign exposing.
- A positive toroidal Cao core has no nonzero azimuthal harmonic in the
  centralizer of its vorticity.  Therefore an arbitrary nonaxisymmetric
  same-leaf profile trace cannot be introduced through that centralizer.
- DA displacements have the stated transverse principal symbol, while the
  full frozen linearized operator retains a distinct Hodge term.
- The supervisor-corrected KKS/angular-momentum sign is preserved.

## Deliberate false-green exclusions

The oracle would fail if the output numerator were hard-coded as zero at
resonance, if `eta_2` lost transversality, if the centralizer determinant were
treated as singular for `ell zeta!=0`, if the Hodge cross-square sign were
reversed, or if `j_2` were assigned the old negative sign.

The oracle intentionally does **not** validate:

- the 0052/0054 finite cluster, graph domain, Riesz projector, commutator, or
  nonnormal resolvent;
- the location or adjoint polarization of an exact Cao critical layer;
- vanishing or nonvanishing of the source-specific Cao recursive trace;
- existence, convergence, pressure reconstruction, or stability of a
  rotating Cao branch; or
- any particle or quantum statement.

Those exclusions matter because a finite algebra oracle cannot substitute for
the full-Hodge critical-layer reduction and inner/outer nonlinear matching
identified in `derivation.md`.

## Post-execution Sol-High scope correction

The exact verifier checks the algebraic generic pendulum width/frequency and
volume-preserving carrier-map identities. It does not construct the full-Hodge
generalized adjoint mode, physical symplectic measure, scaled cutoff/Bogovskii
Sobolev ledger, or a centralizer field with pendulum critical topology.
Accordingly equations (33)--(40) are now classified as a conditional generic
normal form plus exact residual formulation, not an established Euler inner
continuation. No formula predicate changed and the oracle was not rerun.
