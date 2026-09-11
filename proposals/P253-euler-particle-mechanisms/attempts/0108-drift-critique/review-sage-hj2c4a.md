# drift firewall review — HJ2-C4 phase 1 (b3ee8a46): LOW-M PIECE BANKED with two precisions

Reran run_hj2c4a.py (0.5 s): exit 0, 6 assertions, all green. Mutations bite
(attracting-sign flips geometry ✓; wrong δ-scale kills separation ✓). Sign
citation cleared (one-directional flow, no circularity). Verdict: low-m piece
BANKED subject to two required precisions (frozen-δ validity range; absolute tail
summability stated + receipted).

## Pointed (1): budget tier ASYMPTOTIC — frozen-δ validity range REQUIRED

C4-2 proves ratio → ∞ as δ→0⁺ (limit computation ✓ exact-asymptotic). But the
conclusion deploys at FROZEN δ, where the ratio is FINITE: at banked values
(δ≈0.05, L≈3, γ~1) just 2.2× at m=1 — exterior, but thinly. The receipt never
states the validity threshold. REQUIRED (one line + numbers): exterior ⟺ m² >
γδL², i.e. valid for δ < m²/(γL²) (≈0.11 at m=1); banked δ≈0.05 satisfies with
2.2× margin. Asymptotic-plus-threshold, not exact-at-frozen. (Also: C4-2's label
"uniformly exterior" means all-blocks-exterior — true; must not be misread as
uniform-margin — margins grow with m, thinnest at m=1. The validity line fixes
the reading.)

## Pointed (2): absolute accumulation — summability UNSTATED, required

Per-sector bounds (~1/m²δ²) do NOT imply a bounded sector SUM — absolute
summability (Σ_{|m|≥2}1/m² = π²/6−1 < ∞) is the step m=0 isolation needs, and it
appears NOWHERE (paper or receipt). The math is trivially fine (p-series), which
is exactly why its absence matters: one line + one receipt assert (Σ bound
symbolic) closes m=0 isolation against accumulation. REQUIRED (text + receipt
line, minutes of work). Until then "isolation survives accumulation" is asserted,
not shown.

## Pointed (3): sign-citation NOT circular ✓

Flow is one-directional: F-B-lite tensile pre-stress (RB6/RB10, independent
emergent-elastic lane, zero HJ2 input) → C4 displacement sign. Nothing in F-B
depends on HJ2-C4 (or HJ2 at all) — no cycle possible. MB-C4-1 polices the sign
(flip ⇒ opposite geometry ✓ genuine). Tensile⇒stiffening-upward physical
(guitar-string direction ✓). Cleared.

## Verdict

LOW-M PIECE BANKED (rerun green; mutations bite; sign clean) CONDITIONAL on two
precisions landing: (i) frozen-δ validity range with numbers (δ < m²/(γL²),
2.2× at banked values); (ii) absolute tail summability stated + receipted
(Σ1/m² bound). Neither touches the verdict's direction (both strengthen it);
both must land before construction-4 consumption. m=0 + tail exposure map stands
as banked (resonant chain + R5-1 budgets, F-C3 armed).
