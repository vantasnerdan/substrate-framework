# P250 terminal campaign: canonical M5 walls, phase slips, and spherical bags

## Objective and disposition

This PR completes issue #195 at the owner-approved boundary and promotes the
corrected C-M5W-001..008 chain. It constructs the canonical invariant wall
sector, orbit-fixed phase-decoupling mechanism, rigorously isolated local
Maxwell crossing, numerical slice kink and tension, exact reduced bag laws,
seven status-zero spherical slice bags, and the reduced radial stability split.

Issue relationship: `Fixes #195`.

The exact global identification `omega_c^2=omega_*^2` is not asserted or used.
The owner directed on 2026-09-03 that its hour-plus exact computation belongs
to a separate successor campaign and approved this PR without it.

## Independent scientific review

The first review pass found and repaired coupled defects that could not be
fixed safely as prose alone:

- C-M5W-001 inferred that every full-theory wall lies in the slice. Exact
  SymPy differentiation proves the slice is invariant, but the two rotating
  shear Hessians differ and contain `-omega_sq`; the universal quantifier was
  removed.
- C-M5W-003 bounded only the scalar phase gradient. The exact transformed
  tensor-plus-scalar gradient excess is `iota(q')^2/2`,
  `iota=f^2+4b^2`. The accepted statement is an H1 energy-bookkeeping theorem,
  not a global stationary two-frequency solution.
- C-M5W-004 promoted a local Krawczyk/Gershgorin certificate to global
  nonnegativity and exact heteroclinic existence. It now states the rigorous
  local crossing, exact global upper witnesses, and separately sourced
  numerical kink.
- C-M5W-005 treated `delta*iota_B/2` as exact pressure and called
  `F_omega+omega^2 I/2` physical energy. Canonically,
  `p=-V_omega(B(omega))`, `dp/d(omega^2)=iota_B/2`,
  `Q=omega I`, `E0=F_omega+omega^2 I/2`, and
  `E=F_omega+omega Q`.
- The original bag computation reached `solve_bvp` status 1 on six of seven
  rungs, used the endpoint pressure approximation and `omega^2 I` charge,
  checked the wrong radial virial, and clamped a nearly singular translation
  mode. Attempt 0006 replaces it with a moving-coordinate, phase-fixed Robin
  BVP and canonical observables.

Eight individual correction reviews accept the resulting claim statements.
The earlier grouped review records are retained as superseded history.

## Strongest promoted results

- C-M5W-001: exact invariant aligned real-`psi` wall sector, canonical metric
  and inertia, Euler–Lagrange system, first integral, and scoped least-wall
  Jacobi functional.
- C-M5W-002: exact absence of a zero-frequency clock-active bulk state and
  Derrick exclusion of nontrivial static radial shells.
- C-M5W-003: exact orbit-fixed locus and zero mismatch-dependent bulk energy
  for the stated H1 cut-and-paste family, with the full twist density.
- C-M5W-004: exact deep-branch Maxwell system, unique certified local crossing
  at the displayed interval, positive fixed-frequency Hessian there, and
  exact global upper witnesses `omega_c^2<5/3<45/16`.
- C-M5W-005: exact canonical thin-wall Routhian, pressure, charge/energy
  decomposition, `R=2sigma/p`, and `dE/dQ=omega` at reduced-family scope.
- C-M5W-006: status-zero numerical slice kink with
  `sigma_0=0.72929841786(58)` and its itemized error budget.
- C-M5W-007: seven status-zero spherical slice bags for
  `delta=0.001..0.007`; maximum collocation residual `1.00e-8`, boundary
  residual below `1e-10`, radial Derrick relative residual `3.25e-10`, crossed
  tolerance/matching sensitivity, `chi=1.000328..1.002286`, physical
  `Q=2.3430e10..6.8628e7`, and envelope error at most `1.847e-4`.
- C-M5W-008: exact fixed-frequency reduced-radial curvature
  `F''(2sigma/p)=-8 pi sigma<0`, plus six negative numerical Q/R secants; no
  constrained fixed-charge minimum is claimed.

## Artifacts and governance

Reusable APIs and tests live in `src/substrate_framework/m5_wall_clock.py` and
`tests/test_m5_wall_clock.py`. The final campaign record is
`campaigns/P250-shell-bubble-clock/`, including attempt 0006, eight individual
reviews, adjudication, and the content-addressed validation receipt. Corrected
claims remain accepted in `governance/claims.yaml`; release v0.174.0 pins the
scientific correction commit and generated docs/memory agree.

The campaign has no in-boundary debt. Full wall classification, a global
stationary two-frequency interface, global-infimum closure, fixed-charge
constrained stability, full spectra, thick-wall continuation, curvature,
dynamics, force laws, particle interpretation, and gravity remain explicit
future work rather than premises of this PR.

## Validation

```text
attempt 0006 exact_stitching.py: ALL 12 CHECKS PASS
attempt 0006 stability_recheck.py: ALL 8 CHECKS PASS
attempt 0006 bag_bvp_repair.py: 7/7 status-zero rungs; all internal checks pass
tests/test_m5_wall_clock.py + tests/test_m5_exterior_clock.py: 27 passed
scripts/validate.sh --full: 2503 passed; all repository workflow checks pass
git diff --check: clean
```

Validation receipt:
`campaigns/P250-shell-bubble-clock/validation-receipt.yaml`.
