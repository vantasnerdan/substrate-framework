# Validation and exposing-oracle receipt

## Executions

The short exact verifier was first executed with the shell `python`
interpreter. That command, its complete stdout, empty stderr, and exit `0` are
preserved in `finite-window.initial-system.*`. It was then replayed with the
repository interpreter:

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P253-euler-particle-mechanisms/attempts/0048/verify_finite_window.py

The second command and its complete stdout, empty stderr, and exit `0` are
preserved in `finite-window.repository.*`. Both executions printed the same
thirteen passing predicates. After the nonlinear continuation exposed the
internal critical-harmonic scale, two exact predicates were added and the
repository interpreter was run again; that append-only receipt is
`finite-window.repository-v2.*` and contains fifteen passes. This is an exact
symbolic check, not production numerics; no scripts-pane coordination was
needed.

## Checked predicates

The verifier independently expands the exact source-defined rescaled Cao
equation through `q^2` and checks both cell recurrences. It also checks:

1. the exact physical energy normalization `k/sqrt(M lambda)`;
2. the resonance index, coupling, Kelvin spacing, and their `1/log` ratio;
3. the `O(delta^2)` width of one excluded parameter interval and the
   `O(delta/log)` dyadic bad-set measure;
4. the low- and high-index summands in the Schur-tail estimate; and
5. the discriminant sign which distinguishes a same-Krein avoided crossing
   from an opposite-Krein Hamiltonian--Hopf crossing; and
6. `j_crit*Omega_pattern=O(1)` with
   `j_crit=1/(delta^2 log(1/delta))`, while every fixed order remains below
   that scale as `delta->0`.

The functional-analytic statements—DA support closure, Sturm--Liouville Weyl
law, graph-domain Green-kernel expansion, KKS reconstruction, and analytic
Fredholm projection—are proved in `derivation.md`; the script does not replace
them with a finite sampled matrix.

## Exposing failures

- Omitting the `q/(1+qx) partial_x w` term changes both cell equations.
- Dropping the physical energy normalization loses the factor `k/n` and
  falsely makes the accumulated Kelvin modes strongly coupled.
- Using the unearned transfer through `q=1/(4 Omega)` reports `n^-2`; the axis
  form-domain audit in the derivation rejects that shortcut and retains the
  proved `n^-1` decay.
- Comparing the coupling with the bending frequency rather than the local
  Kelvin spacing misses the window-exchange issue.
- Replacing the plus sign in the same-Krein discriminant by a minus sign
  creates a false instability at exact resonance.

## Verdict boundary

The receipts license the exact local Cao second-cell recurrence and the
finite-window scaling algebra. They support, but do not alone establish, the
analytic graph-domain Riesz argument. The attempt establishes the linear
finite-`delta` enlarged spectral construction and every preassigned nonlinear
order `N<=p-2` for a fixed Cao member. Any larger finite `N` requires choosing
a correspondingly larger `p` before the carrier is fixed. It does
not license the exact nonlinear rotating branch: convergence at fixed
`delta` requires the critical-layer transparency/range theorem stated in the
derivation. It licenses no stability, quantum, or particle claim.
