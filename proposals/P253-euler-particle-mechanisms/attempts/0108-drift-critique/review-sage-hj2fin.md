# drift firewall review — HJ2-FIN (e01bbbf4): PIECES BANKED, assembly CONDITIONAL (not closed)

Reran run_hj2fin.py (0.9 s): exit 0, 8 assertions, all green. Distances, tail sum,
pair geometry, and both mutations verified. BUT the assembly step (F-5) re-asserts
F-3's tail identity instead of supplying the missing link: per-sector resolvent
constants (self-adjointness or explicit budgets) remain unshown, so the direct-sum
NORM doesn't close. Verdict: all pieces BANKED; assembly CONDITIONAL on one named
link. F-C3 stays armed. Not a rejection — a one-link downgrade of CLOSED.

## Banked (verified) ✓ — pointed (1) bypass legitimate, pointed (3) geometry matches

- F-1/F-1b (pair inside, distance positive on grid, relative placement exact) ✓;
  inside-iff-γL>1 consistent with C4-3's empty band (m*<1) — both hold at banked
  values, same clean split ✓ geometry intent matches crossing analysis.
- F-2 (sector distances positive on grid, factored form) ✓; F-3 (tail sum exact) ✓.
- MB-F-1 (harmonic-growth kills assembly — summability load-bearing ✓ genuine);
  MB-F-2 (δ=1/2 fence ✓ genuine).
- Bypass legitimate: block-diagonality (round-2, retained orders 0–2) ⇒ sectors
  decouple — no hidden cross-coupling WITHIN scope (higher orders reopen two-index
  per round-2 fence ✓ consistent). The gap below is NOT coupling (bypass fine) —
  it is per-sector constants. Distinguished explicitly.
- Rank preservation via direct sum: valid CONDITIONAL on the sum converging
  boundedly (same condition — see below).

## KEY FINDING: F-5 re-asserts F-3; per-sector constants unshown — CONDITIONAL

F-5's receipt content = F-3's tail identity restated (Σ1/m² finite). But the
assembly needs Σ_m‖R_m‖ bounded, i.e. per-sector resolvent constants C_m with
ΣC_m/dist_m < ∞ — and C_m appears NOWHERE (not in F-5, not in C4a, not in R5).
Distances alone don't bound resolvent norms; tail-of-distances doesn't bound
tail-of-operators. The missing link is exactly one of: (a) self-adjointness route
(sector restrictions of round-1's self-adjoint L_U inherit ‖R_m‖ ≤ 1/dist_m with
constant 1 — IF the sector operators are the self-adjoint restrictions; state +
receipt it); or (b) explicit per-m budgets C_m (F-C3 style) with ΣC_m/m² checked.
Either closes it; neither is present. REQUIRED: add the link (preferred: (a) —
one line + receipt, if the sector operators qualify; else (b)) or re-type assembly
CONDITIONAL on per-sector constant control (named gap, same kind as before).
F-C3 stays ARMED (it tripwires exactly this — operator-side growth in C_m).

## Pointed (2): residue assignment honest + one addition ✓/➕

"Per-block constants are 0052-OWNED inputs" — labeled as dependencies, not
smuggled as established ✓ honest. ADDITION required by this review: per-sector
constant control (self-adjointness or budgets) joins the residue (to 0052's
machinery or a new estimate — state which). Residue = 0052-owned constants (as
stated) PLUS the link above. Nothing new smuggled; one item added to the debt.

## Verdict

PIECES BANKED (distances, tail, pair, mutations — all green and genuine);
ASSEMBLY CONDITIONAL on per-sector constant control (self-adjoint route or
explicit budgets — one link). CLOSED downgraded one step (not rejected):
everything proven stands; the sum needs its constants. F-C3 armed throughout.
Construction 4: assembly-pending-link (not closed, not failed).
