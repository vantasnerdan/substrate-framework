# drift PRE-review — cipher A3 resonance scan (pre-construction constraints, binding)

Status at review: NO scan bytes exist (A3 queued per README#9B66; transfer cleared). This
review therefore freezes the design space BEFORE compute — construction may start subject
to constraints C1–C8 below, which must appear in a frozen A3 design note before numbers
are consumed. Standing precedent: A1/A2 PASS-in-model (0120 review) with disciplines this
scan inherits (analytic-vs-measured labeling per A1; core-model honesty per A2 method-FAIL).

## Binding constraints

- C1 (test-vs-assume): W≈0.2 measured by the scan, never input; Ruban bands (Λ≈4–8 stable,
  overlap Λ≲3, none N≥4) are TESTS with falsifiers pre-attributed to break #1 (transfer).
- C2 (straddle): scan at BOTH Λ-conventions (Saffman 5.08 AND ln(R/a) 3.00) or state the
  choice with reason — one convention silently picked wastes the spend.
- C3 (norm, D-08): freeze the perturbation norm (filament displacement L2? energy?) with
  growth threshold as INEQUALITY (|ρ|>1+tol, tol stated) — a growth verdict without a norm
  is unfalsifiable.
- C4 (core model, A2 lesson): state the filament core regularization (hollow-core? Rankine
  graft?) explicitly — it shifts band edges (= break #1 quantified); silent core = A2's
  invalid-field trap repeated.
- C5 (m-soundness): state m-range/resolution + convergence in m; per-m multipliers
  ρ=exp(±√μT), no aggregate score over m (J1–J4 precedent); base orbit named (PoC-2 T=4.088
  at stated a) + Buttà scaled-up (R,Γ) control orbit.
- C6 (IDEA-07): dense monodromy + NAMED soft subspace in the design (transfer line 18 is
  intent, not a name).
- C7 (scope): verdicts labeled PASS/FAIL-in-model; live-field Euler check always named as
  the uncrossed gap (A1/A2 precedent).
- C8 (receipts): script + log + frozen design committed; seed/param sensitivity per 0120
  discipline (rerun-verified).

## Verdict

GO for construction SUBJECT TO C1–C8 frozen first. No flaw exists yet (no bytes); the
eight constraints are the complete firewall surface — any scan violating them fails review
on landing regardless of numbers. Drift re-verifies at landing.
