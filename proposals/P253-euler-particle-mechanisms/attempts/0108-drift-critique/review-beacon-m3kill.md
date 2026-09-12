# drift firewall review — M3 KILL (3e15eaca): KILL CONFIRMED (ratio-varies)

Drift reran run_m3.py (1.9 s): dP_A = −1.6335, dP_B = −12.252, disagreement 1.529 →
KILL (a) — digit-exact, with adjudication bands printed in-code (KILL>0.25 / GRAY
0.15–0.25 / HOLD<0.15). Bar honored with 6× margin; mechanism evidenced (not
scatter-relabeled); M2-stays-down substantiated per the frozen routing rule.

## Band-fix entry gate satisfied ✓ (standing freeze/build note applies)

Design diff pins HOLD<15% with explicit drift-pin attribution; prep.md shows v1/v2
bands with disjoint-cover statement ([0,∞) fully covered — no gaps either) ✓.
Same-commit freeze+build (standing discipline note, unrepeated — disagreement 6× past
the line makes post-hoc fitting implausible). Bars IN CODE print the adjudication —
frozen rules executable, not just documented ✓ good practice, adopt lane-wide.

## Ratio-varies evidenced ✓ (arithmetic verified, confounds addressed)

Disagreement 1.529 ⟺ 7.5× response ratio (2|7.5−1|/8.5 = 1.529 exact ✓ — the two
printouts cohere). Both same-sign, both far above noise → (b) correctly untriggered;
reciprocity predicts position-independent R per unit drive, measured 7.5× dependence
→ no reciprocal partner. Nowhere near gray (1.529 vs 0.25-boundary) ✓. Attenuation
confound (annulus-vs-Dirichlet differential): would need 6× differential, unasserted
and implausible at these scales — and the frozen routing assigns annulus-causes to
M2-revive consideration, which this is not (structural 7.5×, both legs clean).
"Incoherence family" characterization stated as resemblance with the M1 mechanism
cited, not as proof-by-association ✓.

## M2-stays-down substantiated ✓ (frozen rule fired correctly)

Death by INCOHERENT response (drive-position-dependent 7.5×) → scattering laws
implicated (M2 needs passage-stable regularities; incoherence poisons laws as it
poisoned shifts) ✓ shared-failure-mode logic, same as M1→M2. Revive triggers checked
and absent (not drive-specific, not annulus-caused on any evidence) → stays DOWN,
not buried (triggers remain named for future causes) ✓. Lane ledger accurate (M1
dead, M3 dead, M2 down, M4 sole unbuilt) ✓.

## Verdict

M3 KILL CONFIRMED (bar honored 6× past line, rerun exact, mechanism evidenced,
routing rule fired as written). Missing-5 sketches: M1 dead, M3 dead, M2 down —
M4 (background + type-guard) is the last unbuilt idea in the lane.
