# P253/0031 validation receipt

## Activated and source boundary

The central proposal registers `0031` to `particle-balance-review` with scope
“P4 same-field solitary-wave relative moments, helicity and finite physical
action.” The recorded schema command is

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python scripts/validate_repository.py

and `activation-schema.exit` is zero. The frozen README hash is
`ba377906ac1f4e8b75444815c7389f6b86a92afa78bb16c16fc296a26bd7defe`.
The activation command, stdout, stderr, and exit hashes are respectively
`5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

The load-bearing supplier files are `0027/solitary-wave-construction.md` hash
`c18d394a757f9fd967975ba23a486e62087c285c8fda4f08d8aeaa0159c6f850`
and `0027/exterior-construction.md` hash
`0c2cda195dd5bea9163a7a6c582bfcff746a608a2e48d0a9cd7a007c1b260fc3`.
Their independently established scope is recorded in `0028/review.md` hash
`91fa6ffc55c43acec1d87f05b18100e6f0d30f0a18734f44aa6696399601ddc0`
and `0028/verdicts.yaml` hash
`cba564a7e38818d6fa1c156cfb08f0b06ada5d502236373bba8ac66c5076d3c9`.
Current `0027/result.yaml` and `validation.md` only record that earned review
state; the proof/API/test boundary used here is unchanged. No `0029` or `0030`
body was inspected.

## Strongest analytic oracle

The primary validation is the derivation itself:

- direct cylindrical curl gives the complete laboratory vorticity and raw
  helicity density;
- weighted Green integration retains the axis, axial-cap, and radial-infinity
  helicity surface terms and proves their limits;
- the direct axial angular momentum is absolutely convergent, while its
  vorticity rewrite retains the radial surface row;
- axial impulse is derived from the toroidal vorticity and the exterior
  zero-Fourier mode, which produces the nonzero infinity term;
- the full background/swirl/poloidal kinetic excess and finite-time material
  action retain `rho_m` and the cylindrical measure;
- the KKS sign is propagated through the Euler Hamiltonian and both spatial
  momentum maps on the fixed-background orbit; and
- weighted convergence to the exact homoclinic gives controlled little-`o`
  expansions with strictly nonzero coefficients, rather than replacing exact
  finite-`mu` integrals by leading formulas.

The analytic result also separates the finite local KKS/action construction
from the absent compact internal orbit. Axial `SO(2)` is a stabilizer, axial
translation is a noncompact `R` orbit, and tilted rotations leave the
finite-excess background class. This is a source/domain conclusion, not a
numerical no-go.

## Executable history and method repair

The first invocation used a nonexistent worktree-local interpreter,

    PYTHONPATH=src .venv/bin/python proposals/P253-euler-particle-mechanisms/attempts/0031/verify_relative_observables.py

and exited `127`. Its exact command, exit, and stderr are retained as
`initial-environment-run.*`. This was an environment-path error and made no
scientific assertion.

The next invocation used the repository environment and derived four identities
before SymPy declined to simplify the improper `sech^2` integral, exiting `1`.
The exact command/stdout/stderr/exit are retained as `initial-symbolic-run.*`.
The method repair replaced that unsupported direct integration with an exact
antiderivative check and its two endpoint limits; the integrand and claimed
factor were not changed.

At the final unchanged script boundary, the command

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P253-euler-particle-mechanisms/attempts/0031/verify_relative_observables.py

exited zero and printed seven exact checks plus four exposing mutations. The
final command, stdout, and exit hashes are
`cff115140be4b7c205a678537db007ec762643c2c6e7c951ebd31fdf7479cc00`,
`c5b42dce0f83519f22ee4b9978698d56f4525a788dd08d1bade8c9ed107ddea2`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
The mutations reverse the helicity boundary sign, remove half the azimuthal
factor, drop one impulse surface coefficient, and reverse the KKS sign; each is
algebraically exposed.

The script corroborates exact local identities and factors. It does not prove
absolute convergence, the exterior zero-mode estimate, global group
integration, or the nonlinear existence theorem; those roles remain with the
analytic proof and the independently reviewed supplier.

No source/API/test file changed, so the existing `0027` four-test receipt was
not rerun. No production numerical eigenvalue, force, energy splitting, or
stability edge is used, and the small-ratio-numerics prescriptions do not bind.

## Bounded Hodge-domain proof correction

A later bounded proof check identified one false intermediate premise in the
momentum-map bridge: a compact coadjoint generator makes
`delta omega=curl(eta cross omega)` compactly supported, but the recovered
velocity tangent `delta u=P(eta cross omega)` generally has a Hodge pressure
tail. The observable values, asymptotics, KKS sign, and route verdicts are
unchanged.

The corrected proof differentiates impulse through the compact field
`eta cross omega` itself and obtains

    delta I_z=rho_m integral (eta cross omega)_z
             =rho_m integral eta dot partial_z u
             =Omega(X_Z,X_eta).

For axial angular momentum, the Leray gradient has zero pairing because its
complete azimuthal-cylinder integral is `integral partial_theta p=0`. The
result now states only the finite local KKS form and weak Hamiltonian identity
on compact-generator tangents, together with the exact noncompact translation
orbit. A general noncompact Euler generator is not placed in `Diff_vol,c`, and
an arbitrary-nearby asymptotic group remains conditional on a separately
specified domain.

This correction changed no executable identity or source/API/test input, so
the unchanged symbolic oracle was intentionally not rerun. Exact pre/post
artifact hashes are recorded in `hodge-domain-correction.md`.

## Final artifact checks

At the final result boundary,

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python -c
      'import pathlib,yaml; ... yaml.safe_load(...); assert both route verdicts'

exited zero and printed `0031 result YAML parse PASS`. A scoped
`git diff --no-index --check /dev/null <file>` pass over every regular file in
`0031` also exited zero after interpreting the ordinary new-file status
separately from whitespace errors. The derivation, structured result, verifier,
and final symbolic stdout hashes at that boundary are respectively

    8a44074412f41b3fcee94dcd5154f87ab060adca51024af29fb38adc541002cd
    e2bde6a06bca0dd968db1032e0d708e6d97aa0184e790d8650ffd4e8c437298e
    08c6f7bd218d485cda047188648a11bb236ace2bdc1c81bdbbce53e0a1ba590a
    c5b42dce0f83519f22ee4b9978698d56f4525a788dd08d1bade8c9ed107ddea2.

No repository-wide test or claim-registry replay is indicated by this
append-only proposal artifact.
