# Attempt 0024 — N7 EPS bridge established: realized Beltrami example matches declared premises, verify_cst007 12/12

## Route

EPS bridge (obligation N7): the declared existence layer is imported, not
re-proven (sources in hand, md5-pinned per PROVENANCE.md and the proposal
source_inventory); the realized explicit example (Chandrasekhar-Kendall-type
Beltrami eigenmode) is checked against the declared ensemble premises in
sympy. Consolidated in `verify_cst007.py` (12 checks, 3 mutations, exit 0).

## Result

Realized example u = (sin(lam z), cos(lam z), 0), all checks sympy-exact:
- (a) incompressibility: div u = 0
- (b) Beltrami property: curl u = lam u
- (c) helicity density: u . curl u = lam |u|^2 (conserved invariant)
- (d) stationarity: (u.grad)u = 0 with constant pressure -- an exact
  stationary solution of constant-density incompressible Euler

Source integrity: 1210.6271.pdf md5 matches the provenance record; the full
declared in-hand set (1003.3122, 1505.01605, 2103.14458) present.

## Checks and mutations

- (a)-(d) structural checks above.
- M1 non-Beltrami perturbation: rejected (curl u != lam u).
- M2 wrong-sign helicity: rejected.
- M3 md5 tamper detection: functional.

## Status

- route_verdict: established (L5: EPS existence input imported from in-hand
  verified sources; the realized Beltrami example matches declared premises)
- evidence_scope: EXACT (sympy) + provenance (md5)
- THE OBLIGATION LADDER IS COMPLETE: N1-N7 all established, L1-L5 earned
  (L1's B-constant remainder carried as the declared log-running bend
  stiffness with c1 = 1/2 - EulerGamma per 0015/0017). Claim review and
  terminal PR follow.
