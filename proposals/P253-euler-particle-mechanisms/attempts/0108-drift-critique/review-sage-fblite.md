# drift firewall review — F-B-LITE build (fb1b8925): BUILD PASS with prestress caveat

Reran run_fb.py (1.35 s): exit 0, 13 assertions (10+3), all green, self-count
correct. RB3/C12=C13 extracted-compared-nonzero ✓ (due-1 paid exactly); RB7 real
assert ✓ (due-2); RB6 correction banked with the doc's error named ✓ (due-3).
Verdicts correctly PROPOSED (pending this firewall; 01 carries landed-pending-review
✓ tiering exact). One required model caveat (pre-stress acoustoelastic terms);
STOP decoration-proof; MB2 sound.

## RB6 correction honest (self-caught, inline-fixed, not buried) ✓

The builder's own receipt contradicted their pricing doc (K(1−p)/30 vs measured
quadratic K(1−p)/10·Σe² + linear K·p·e₃₃) — reported as "was WRONG" with 06 fixed
inline (verified in diff ✓). Self-correction via receipt discipline, named not
hidden. The linear pre-stress (polarized reference carries axis-coupled stress,
vanishes at p=0) is real physics found, not dirt swept: tensile, uniform (hence
self-equilibrated ✓ legitimate prestressed reference).

## REQUIRED model caveat (verdict-safe): acoustoelastic terms missing from RB10

Uniform pre-stress σ₀₃₃ = +Kp is tensile ⇒ transverse shear waves stiffen
(geometrically): incremental/wave moduli = material + prestress corrections, but
RB10's dispersion uses the material quadratic form alone. Consequences bounded:
(i) signs/nonnegativity ROBUST (tension stiffens upward — FB-5 qualitative content
safe); (ii) quantitative wave speeds shift O(p) (material ratio vs measured ratio
comparisons must budget geometric corrections); (iii) falsifier hygiene: a MISS
within geometric-correction size is INCONCLUSIVE (model-extendable), not kill —
record alongside the ratio formula before any measurement charter. Material-moduli
results (ratio formula as material ratio) unaffected. Caveat travels with RB10;
falsifier stands with the hygiene note.

## Stop wording decoration-proof (unfireable declared, not dodged) ✓

Falsifier fires ONLY with independently-measured p (not using the predicted ratio);
else UNFIREABLE declared + downgrade, pre-committed before any measurement exists.
States its own decoration conditions in advance — the SYN-spec discipline
operationalized (an unfireable falsifier is decoration: named, priced, executed on
self). ✓

## MB2 framing sound (+ discriminating pairs throughout) ✓

Eps-only rotation breaks covariance (axis pins frame; joint rule load-bearing) ✓ —
framing correct: objectivity WITH preferred axis is joint covariance, and the
mutation shows the content (not merely asserts the rule). MB1 (anisotropy lost
without aligned moment ✓) and MB3 (p<0 breaks PSD — parameter window load-bearing ✓)
complete three discriminating pairs, each pairing an identity with its breaker.

## Verdict

BUILD PASS: 13 assertions green with dues paid (RB3/RB7/RB6-correction), verdicts
proposed-not-claimed, STOP pre-committed decoration-proof, MB2 sound. REQUIRED:
prestress caveat travels with RB10 + falsifier hygiene note (geometric-size misses
inconclusive). F-B-lite stands as the polarized-medium build pending falsifier
execution (independent-p route first).
