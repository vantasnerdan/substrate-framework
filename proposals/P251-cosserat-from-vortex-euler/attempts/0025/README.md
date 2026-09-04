# Attempt 0025 — independent PR review repair: operator-derived dispersion, scoped ergodic closure, and reopened N2→N3 bridge

## Boundary

- Reviewed PR: #199 at `16dae73077d219cb8d864384fd1e07adae279967`
- Base issue: #198 (created before the PR)
- Review role: independent non-author correction pass
- Small-ratio prescriptions: the recorded `c4 = -0.253 +/- 0.015` estimate is
  not promoted here because the campaign does not supply its crossed refinement,
  evaluator-noise floor, or derived error budget. The repaired N5 IVP uses a
  well-conditioned mode and records float64, DOP853, tolerances, refinement,
  state error, and energy drift; no soft-eigenvalue sign claim is made.

## Strongest retained results

- N1: the Bishop/Frenet frame-transport identities remain exact and unchanged.
- N2: the Poincaré reduction, inertial-wave limit, and recorded asymptotic
  Rankine-vortex branch calculations remain useful conditional evidence.
- N3/N4: the sphere-moment contractions, coefficient matching, and
  Euler–Lagrange identities are exact *conditional on the supplied micropolar
  energy*.
- N5: the corrected conditional transverse determinant and branch limits are
  exact, with an independent SciPy IVP replay.
- N6: uniform phase decorrelation kills the coherent signed first moment; the
  Monte Carlo routine now actually executes and tests that same observable.
- N7: the hash-pinned EPS sources supply the declared vortex-tube existence
  input, while the simple periodic Beltrami field is retained only as an exact
  Beltrami-to-steady-Euler regression.

## Repairs

1. `verify_cst005.py` now builds its determinant from the N4 2×2 operator.
   The coupling product is `4 alpha^2 k^2`, not `8 alpha^2 k^2`. The previous
   doubled value made the N3 parameter specialization unstable at small `k`.
2. The `L_v -> 0` limit is evaluated with `mu, alpha, c_*, j` all proportional
   to `L_v`. It has one neutral displacement root and a finite spin ratio; at
   `L_v = 0` the zero-weight spin coordinate is removed. The old substitution
   `j = 1` tested a different model.
3. The promised DOP853 IVP replay was added with a three-level tolerance
   refinement and an energy-drift monitor.
4. `verify_cst006.py` now calls its Monte Carlo check and replaces literal
   `True` predicates with exact sphere/phase integrals. It separately verifies
   `<L> = 0` and `<L.T L> = I`, preventing zero mean from being mistaken for
   zero quadratic fluctuation energy.
5. `verify_cst007.py` checks all three EPS source digests, audits theorem
   statement locations, and labels the constant-speed plane-wave field as
   nonlocalized rather than as a realized vortex tube.
6. The proposal graph now reflects its own dependency ledger: L1–L5 are
   unearned, N2–N7 are active, and terminal success is revoked pending the
   missing microscopic coupling construction.

## Load-bearing bridge finding

The current N3 `alpha` coefficient is obtained by replacing the exact relative
deformation with `h - skew(Phi)` at first order and then retaining its quadratic
norm. Exact Green–Lagrange line stretch cancels a free frame rotation because
`R.T R = I`; straight-tube tension alone therefore cannot supply the uniform
relative-rotation energy. The exact moment matching is reusable, but an
Euler-derived frame-locking or core-polarization interaction is still needed
before `alpha = L_v T/6` is a derived modulus. This agrees with the manifest's
pre-existing `L1: unearned` state and with attempt 0019's note that the
`verify_cst002` velocity forms fail the momentum residuals.

## Route verdicts and next construction

- N5 repair route: **established conditionally**; evidence scope `EXACT` plus
  well-conditioned `NUMERIC_REGRESSION`.
- N6 repair route: **established as scoped** for the coherent signed-response
  closure; it does not claim vanishing microscopic quadratic fluctuations.
- N7 repair route: **established as provenance + equation regression**; it does
  not claim the elementary field is an EPS vortex tube.
- N2→N3 candidate A bridge: **blocked** by the missing Euler-derived
  frame-locking interaction. This is a route verdict, not an obligation no-go.
- Next executable route: rederive N2 with the momentum-residual-correct velocity
  forms from attempt 0019, then compute the candidate-B core-polarization or
  inter-tube Biot–Savart interaction that couples macro rotation to the
  independent director. Success supplies L1 and makes the exact N3/N4/N5
  conditional machinery promotable without changing the user's objective.

## Commands

```text
PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P251-cosserat-from-vortex-euler/verify_cst003.py
PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P251-cosserat-from-vortex-euler/verify_cst005.py
PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P251-cosserat-from-vortex-euler/verify_cst006.py
PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P251-cosserat-from-vortex-euler/verify_cst007.py
PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python scripts/validate_repository.py
PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python scripts/validate_changed.py --base 16dae73077d219cb8d864384fd1e07adae279967
ruff check proposals/P251-cosserat-from-vortex-euler/verify_cst00{3,5,6,7}.py
git diff --check
```

All commands exited 0 at the repaired boundary; concise tallies and numerical
receipts are in `stdout.txt`, and stderr was empty.
