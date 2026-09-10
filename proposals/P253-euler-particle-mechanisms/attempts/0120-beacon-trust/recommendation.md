# Recommendation: (a) fitted mesh vs (b) error-bar acceptance (beacon)

Question shepherd must rule: how does G-a2 close — (a) free-boundary-
fitted mesh rebuild, or (b) accept the approximate trust member with
residual-propagated error bars?

## Costed evidence

- (a) fitted mesh: cost ≈ 3+ obligations (cut-cell/r-adapted mesh fitted
  to src=0 contour, re-solve bordered system, re-verify FD/Jacobian/
  rows, re-feed). Payoff: lifts observed order p≈0.4 toward ~2 →
  potentially 10–100× residual reduction → tol-level member. RISK:
  formulation risk repeats (focusing p=6 + free boundary + border on a
  NEW mesh; the merit-plateau mechanism is mesh-independent per 0118,
  so fitting attacks only the quadrature half of the floor).
- (b) error-bar acceptance: cost ≈ 1–2 obligations (Maxwell O(g) stage
  on the trust member + Q/H perturbation lemma in δω + production feed
  with budget). Structural facts favoring sufficiency:
  (i) Q_F is EXACTLY translation-invariant (P_L is a Fourier multiplier;
  one-line proof) — the solver's known soft direction does not move the
  QOI at all; (ii) current margin λ≈11.2 vs plausible error O(0.1–1)
  from the 6.6e-3 floor (estimate — the lemma must replace it);
  (iii) G-a consumes INEQUALITIES (σ_min ≥ threshold with margin), never
  equalities. Shared prerequisite: Maxwell stage either way.

## Recommendation: (b) first

(b) is cheaper, uses tested code paths, and its first step (the lemma)
DECIDES the question: if margins hold, G-a2 closes without (a); if the
lemma fails (error ≥ margin), (a) is then justified with a quantified
target instead of a hope. (a)-now spends 3+ obligations before learning
whether tol-level fields were ever needed. Awaiting ruling before
building (per order).
