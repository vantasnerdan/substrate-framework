# Perturbation lemma with measured gap (beacon 0121)

Claim to decide: do the trust-member witness constants hold with margin,
i.e. δλ < λ_ω = 11.2, from the measured residual alone?

## Measurement (dense eigendecomposition, free block, 800 dofs)

- Spectrum: one negative (−0.51, separated saddle — expected, harmless
  for the bound), soft cluster 0.016–0.075 (8 modes), bulk to 43.
  ARPACK-SM fails on this spectrum (8000 iters, 0/6); dense succeeds.
- Residual soft-projection share: 50% of ‖R‖² sits in the 8 softest
  modes (amplification ×15–60).
- Linearized error: ‖δu‖_rms = 0.047 (3.7% of umax) via full eigen-
  expansion (all modes incl. negative; no deflation needed — nothing
  is within 1e-2 of zero).
- δω transport: ‖δF‖_2 = 2.08 via r³jf²-weighted mass form (jf ~ O(10³)
  in the core — the dominant amplifier).
- Bound: δλ ≤ 2·s_ω·‖δF‖ = 2(5.92)(2.08) = 24.7 vs λ = 11.2.

## Verdict: LEMMA FAILS (margin ratio 0.45 < 1)

Per the (a)/(b) ruling, (a) ACTIVATES with quantified target: cut
‖δF‖ 2.08 → ≲ 0.1 (≈20×) AND halve the soft-mode residual share.
Brute refinement cannot (p≈0.4 needs 6 halvings); fitted mesh or
soft-aware solver (deflation/smoother targeting the 8-mode cluster)
is required. First-order analysis; reg-smoothing; translation modes
provably harmless (Q exact-invariant, separate one-line proof).

## Maxwell O(g) stage: DEFERRED with reason

G-rows inherit the same δω-driven bound (B_g via χ_g(member) +
witness assembly); building them now spends work behind the same
failed inequality. Build Maxwell AFTER (a) delivers a member meeting
the quantified target. No work lost: L_c-Poisson path + tag
approximant spec stand ready in 0111/G-a2 docs.
## Appendix (drift repair): translation invariance, written not softened
Claim: Q_F is exactly invariant under spatial translations of F.
Proof: let (T_a F)(x) = F(x−a). Cross with a constant vector commutes
with T_a, and the Leray projector P_L is a Fourier multiplier, hence
also commutes: S_{T_a F} y = P_L((T_a F)×y) = T_a P_L(F×y) = T_a S_F y.
T_a is unitary on L², so Q entries ⟨S_i,S_j⟩ are preserved. ∎
Scope: translations only (rotations act nontrivially but unitarily on
the triple — Q conjugates by the rotation, spectrum preserved too).
This removes the translation null direction from the lemma's error
budget exactly; all other soft modes remain in the measured bound.
