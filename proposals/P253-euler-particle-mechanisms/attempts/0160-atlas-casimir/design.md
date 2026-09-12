# 0160-atlas-casimir — Casimir functional + discretization SCOPE (atlas)

Charter: shepherd sequencing decision (owner-delegated, recorded) — climb
the 0155 wall (0155-atlas-chart build-report: energy–Casimir functional +
discretization ABSENT from banked code). Serves the R-B coercivity verdict
and the R-C C2 falsifier (0150 design: Hessian indefinite → C2 dead).
Falsifier + stop FROZEN here pre-compute. Banked intermediates, in order:
(1) functional candidate, (2) discretization, (3) Hessian estimate.

## Bracket identification (what the reduced model is)

Banked model (0120 run_a3.py, read-only): 2 rings × N nodes, Saffman-local
self-induction + Rosenhead–Moore mutual, axisymmetric leapfrog orbit.
The continuum parent is 3D Euler (vorticity Lie–Poisson); Casimirs of 3D
Euler include helicity ∫v·ω and generalized enstrophies. The reduced ring
model inherits candidate invariants: ring impulse (momentum map — NOT a
Casimir, do not mislabel), ring volumes/areas, and the axisymmetric
Casimir family ∫F(σ) with σ = ω_φ/r (potential-vorticity-type tag).
Scope step 1 decides which of these survive as Casimirs of the REDUCED
bracket (not the continuum one) — continuum Casimirs that the
discretization breaks are DECORATION, named as such.

## Construction (frozen, three banked intermediates)

1. FUNCTIONAL: write the reduced-bracket Casimir candidate C(state) on the
   (R,Z,m-mode) collective coordinates; prove {C,·}=0 in the reduced
   bracket on paper (exact algebra) or record BRACKET-ABSENT if the reduced
   model has no Hamiltonian form to check against.
2. DISCRETIZATION: C evaluated on banked orbit states (chart-table.txt +
   Mono npz states); consistency check vs continuum limit at N=64/128
   (two banked resolutions — divergence across resolutions = FAIL).
3. HESSIAN: second variation of (energy − Casimir) on the member orbit in
   the comoving chart; definiteness verdict feeds C2 (definite → C2 lives;
   indefinite → C2 DEAD per 0150) and R-B (lower bound or wall).

## Falsifier (frozen pre-compute)

- Hessian indefinite on the member orbit → C2 DEAD + coercivity DEAD;
  functional/discretization stand as banked intermediates (reusable).
- No reduced-bracket Casimir identifiable (BRACKET-ABSENT) → build STOPS
  at intermediate 1; wall stands as priced; C2 stays open-pending (not dead).
- Resolution divergence (N=64 vs 128 disagree on C) → discretization FAILS;
  no smoothing, no tuning — record and stop.

## Stop (frozen pre-compute)

- Any intermediate fails its gate → STOP at that intermediate, bank it,
  no second rounds, no repair loops.
- New solves: NONE licensed beyond banked orbit states + short linearized
  flows with frozen machinery (same metrology class as sep_probe.py).
  Full re-shoots/re-monos need a new charter.

## Cost (frozen)

Paper + banked-number compute. Drift review invited pre-compute; build
runs on shepherd charter (given), findings fold on arrival.

## Build record (below freeze line)

- 2026-09-11: scope frozen, drift review invited, 0160 claimed in ID-LEDGER.
