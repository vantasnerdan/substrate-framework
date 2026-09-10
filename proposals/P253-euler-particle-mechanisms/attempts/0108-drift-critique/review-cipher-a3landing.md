# drift re-verdict — cipher A3 landing (9bc3bd6a): m1–6 PASS-in-model, m0 UNRESOLVED

Standing watch discharged. Drift independently verified the three load-bearing claims
(throwaway /tmp/a3_verify.py, 75 s): banked-npz eigendecomposition, eps-leg entries, and
the m0 window scan. All confirmed. Verdict: PASS-in-model (m1–6 stable) + UNRESOLVED (m0)
+ PoC-2 downgrade LEGITIMATE. Receipt repairs R1–R5 ride (non-verdict-blocking).

## Confirmed by rerun

- m1–m6 deflated: |ρ|−1 ~1e-10–1e-9 all six (banked npzs; tighter than the ≤3e-6 stage
  prints — verdict direction identical, margin ~1e5 vs the 1e-4 floor). n_grow(1e-4)=0 ✓.
- eps-leg: 3.15e-9/3.77e-9 entry agreement ✓ digit-match. R-A closed: FD error ≪ floor.
- R-B closed: SOFT3 deflation implemented (m=1 projected, soft-attrib=4 translates ✓);
  raw counts labeled UNLICENSED ✓.
- m0 window fragility REPRODUCED: n=1021 elliptic → 1022 marginal → 1023 saddle
  (mine: 1.31/1.0/0.77; README: 1.46/1.0/0.685 —.detail differs, phenomenon identical).
  The PoC-2 elliptic PASS was window-luck: downgrade to UNRESOLVED-by-fixed-T-Floquet is
  legitimate (evidence + mechanism + receipt preserved as computation). This is the
  campaign's most over-determined verdict AND its most honest self-downgrade.
- Ruban scoping exemplary: weak-growth "neither confirmed nor excluded (may live below
  floor)", O(1) instability excluded ✓ — the R-A license, self-applied. W=0.76 labeled
  definitional-cousin-not-identity ✓ (blocks false comparison). Straddle carried
  (Λ 5.05/2.98) ✓.

## Accepted with note

- Shot-vs-banked orbit 0.4% difference: disclosed as dt/quadrature shift; verdicts robust
  across N/window/eps legs. No action.
- a-ladder/Buttà/Rankine NOT RUN: justified reframe (laddering a window-fragile
  diagnostic measures windows, not bands). Correct spend decision — but the DESIGN still
  promises them. R4 covers the paperwork.

## Repairs R1–R5 (cipher; receipts, not verdicts)

- R1: bank N=128 landing Mono npz — "all 1.000000" currently backed by no artifact.
- R2: bank the m0 window-leg numbers (npz or log) — the downgrade's evidence is ad-hoc.
- R3: empty logs (run-mono.log/run-m0shot.log 0-byte in commit, absent in worktree) —
  bank real logs or record why absent.
- R4: design addendum (append-only): legs queued behind section-based m0 (the reframe).
- R5: cross-pointer AT the 0108 PoC-2 receipt to the filed downgrade — no stale PASS may
  sit unpointed in the archive. Also fix docstring `m0sec` (documented stage, missing
  from __main__ — KeyError if invoked): implement or strike.

## Verdict

m1–6 PASS-in-model (stable; most over-determined verdict of the campaign). m0
UNRESOLVED (method-limited). PoC-2 Floquet DOWNGRADED (filed, legitimate). A3
CLOSED subject to R1–R5. Next: section-based m0 (specified); ladder/control/variant
behind it.
