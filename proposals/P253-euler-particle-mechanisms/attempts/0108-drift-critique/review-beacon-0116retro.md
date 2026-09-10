# drift firewall review — beacon 0116 unblock spec (RETROACTIVE, closes ledger L-1)

Context: shepherd's Route-A ruling consumed 0116's acceptance; 0117/0120/0121/0122 derive
from it — all reviewed by drift WITHOUT 0116 itself being reviewed. This closes that gap.
Drift verified skfem 12.0.2 present in-venv; read repro + request + receipts at source.

## Spec soundness — SOUND on all three parts

- Repro exact: capability matrix (present numpy/scipy/sympy/skfem; absent heavy solvers;
  grep-probe for in-tree GS solver; find-probe for member data) with correct exit-1 logic
  (blocked iff heavy-missing AND no solver). The "key surprise" (skfem already in-venv,
  verified by drift) is factual and load-bearing — it made Route A a no-install task.
- Route framing honest: Route A scoped to 0077 (21)–(23) + L_c Maxwell + τ=g² IFT at fixed
  (κ,I_z) in 0080 window (37) under small-ratio-numerics — exactly the build 0117 executed
  (staged: GS first, Maxwell deferred per 0121 — staging discipline held downstream). Route B
  costed (system deps + build + validation + owner approval) with no version pinned until
  routed. Recommendation explicit, decision left to shepherd/owner. ✓
- Acceptance ladder well-posed: member arrays + norm certificates → pipeline feeds → rows
  numeric with error bars → G-a2 DONE (exit-0 flip defined); G-b downstream named. The
  ladder's norm-cert/error-bar requirements are precisely what caught every downstream
  difficulty (stall, branch-hop, floor, lemma margin). ✓

## Retroactive verdict + inheritance — PASS, chain stands

No spec defect propagated: free-boundary order (p≈0.4), branch non-uniqueness, narrow basin
are execution DISCOVERIES within the ladder's discipline, each kept BLOCKED with mechanism —
exactly what the acceptance criteria demand. 0117/0120/0121/0122 inherit a sound license;
no retroactive invalidation, no re-verdict of downstream work needed. U3 (0114-repair
verification) independently corroborated by drift's own reproductions. No repairs.

## Verdict

Retroactive PASS. Ledger L-1 CLOSED. The unreviewed-consumption pattern itself: process
note to shepherd — rulings consuming acceptance should route the spec through firewall
 FIRST (as now done retroactively); no harm this time because the spec was sound.
