# G-a2 unblock request (beacon 0116, for shepherd ruling/routing)

Repro: `ga2_repro.py` (exit 1, run receipt in `tool-receipts.md`).
Finding: numpy 2.5.1 + scipy 1.18.0 + sympy 1.14.0 + **scikit-fem 12.0.2
present**; firedrake/dedalus/dolfin/mpi4py/petsc4py absent; NO in-tree
axisymmetric Grad–Shafranov solver; NO member field data anywhere.

## What fails without it

`ga-status.md` rows 1–5 (Q/λ/s/R/ε/H/G numerics) + 8–9 ((25a) orbit, K
blocks): the G-a1 pipeline is proven working on synthetic data but
input-starved. Unit G completion stays blocked at G-a2; G-b consumes the
curve downstream.

## Exact build requested (either route; shepherd rules)

- Route A (in-venv, recommended): frozen-design numerical task on
  scikit-fem 12.0.2 (present, no install): axisymmetric r-weighted
  Grad–Shafranov discretization of 0077 (21)–(23) + L_c Maxwell solves +
  τ = g² IFT continuation at fixed (κ, I_z) in window 0080 (37), governed
  by `small-ratio-numerics` (error budgets, λ-floor, observed-order).
  New packages: NONE. Command when designed:
  `PYTHONPATH=src .venv/bin/python attempts/<ga2>/build_member.py`.
- Route B (heavy): approve + install a full 2D stack (FEniCSx or
  firedrake + MPI/PETSc system deps, or dedalus + FFTW/MPI). Cost:
  system packages + long build + new-env validation; needs owner
  approval per repo policy. No version pinned until routed (latest
  stable at build time + hash lock).

## Acceptance (turns repro exit 0)

Member field arrays (ω_g, B_g, χ_g + grid/metrics) land in-attempt with
norm certificates → `ga_pipeline.py` feeds → ga-status rows numeric with
error bars → G-a2 DONE, G-b unblocked.

## Minimal repro

`.venv/bin/python proposals/P253-euler-particle-mechanisms/attempts/0116-beacon-ga2unblock/ga2_repro.py`
→ exit 1 (BLOCKED line above). No implementation in this attempt by
shepherd order.
