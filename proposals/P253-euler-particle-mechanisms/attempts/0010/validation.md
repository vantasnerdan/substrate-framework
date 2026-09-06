# Validation receipt

## Boundary

Only `proposals/P253-euler-particle-mechanisms/attempts/0010` was written. No source,
canonical governance, central proposal, skill, memory, or generated-document file was edited,
and no commit was created by this provenance repair. The historical attempt baseline is commit
`e38a8e95ed04269de43f6ce067d05b190628585c`.

## Strongest analytic oracle

The main new bridge was checked by independent elimination and square completion. For

    f(xi)=xi^2/(2r^2)+(alpha/e)xi-e zeta xi,

symbolic differentiation returned

    xi_* = r^2(-alpha+e^2 zeta)/e,
    f(xi_*) = -r^2(-alpha+e^2 zeta)^2/(2e^2).

The same symbolic check returned zero for the residual in

    chi^2/r^2-2e eta chi
      = (chi/r-e r eta)^2-e^2r^2eta^2.

This independently reproduces Cao--Zhan's
`J(r,s)=r^2(s-alpha)_+^2/2` after positivity restriction and verifies that its
reduced Hessian is the Schur complement of equation (17) in `derivation.md`.

The exact historical Python one-liner is now materialized source-equivalently as
`verify_variation.py`. Its original stdout is preserved in
`original.stdout.reconstructed-from-transcript.txt`; the filename deliberately records that
the bytes were recovered from the durable Codex session transcript rather than produced by a
new execution of the materialized file. In that output, symbols `a` and `z` denote `alpha` and
`zeta`, respectively.

The physical dimensions provide a second check:

    [H T] = [rho0 U^2 L^3][L/U] = [rho0 U L^4] = [j],

consistent with the KKS period `4pi j`. The result leaves both Euler scale parameters free;
no quantum action was fitted.

## Commands and exits

- Frozen activation receipt: `activation-schema.exit` contains `0`.
- The exact original composite command is preserved in
  `original-command.reconstructed-from-transcript.txt`.
- Codex session `01a07656-64cf-74e2-ac94-9ba86652daef`, command ordinal `965` and execution
  result ordinal `967`, records the displayed SymPy output, empty stderr, and composite exit
  `0`; `transcript-provenance.txt` records that lookup.
- The historical command joined the Python one-liner, `git diff --check`, and `git status` with
  semicolons. Therefore its process exit `0` is the exit of the final `git status`; separate
  Python and `git diff --check` exit codes were not captured and are not reconstructed here.
- `verify_variation.py` was not rerun for this receipt repair: the recovered transcript is the
  historical execution evidence, not a newly manufactured tally.
- No small-ratio numerical verifier was designed or run, so the small-ratio-numerics skill
  was not activated.

## Claim boundary

The oracle validates the algebraic correspondence and action scaling, not nonlinear Euler
stability. Source-level existence statements remain external imports in this active attempt.
The exact open proposition is coercivity (or a growing mode) of the joint Hessian on the
dynamically accessible Euler tangent space, modulo translation and then in non-axisymmetric
sectors.
