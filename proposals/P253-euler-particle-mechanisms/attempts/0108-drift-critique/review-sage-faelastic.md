# drift firewall review — F-A + CONTRAST (38512dbb): FORMULA-HALF PASS with model caveats

Reran run_fa.py (0.65 s): exit 0, all green, self-count CORRECT at 14 (12+2 — the
0147 miscount not repeated ✓). Independently re-derived the |Fn| expansion from
F=RU polar decomposition. Verdict per half as chartered: F-A formula-half PASS
(in-model, μ>0 robust) WITH two required model caveats; window-half PENDING (pricing
honest); CONTRAST model-sound with one REQUIRED receipt repair (R6b vacuous as
computed); LANE-1 PENDING (correctly typed already). F-C charters independently (below).

## Second variation: FORMULA exact, MODEL footnoted (watch-item 1 answered)

Receipt algebra (R1→R2→R3) verified exact; ε_s-coefficients (1/10, −1/30) re-derived
independently and match; R4 tensoriality + MA-1/MA-2 audits sound; R7/R8 follow.
μ_aff = K/10 > 0 STANDS. Two footnotes ride (verdict-safe — μ sign/magnitude-class
unaffected):
- (a) ε-identification: the expansion is exact for BIOT (rotation-free) strain U−I;
  the doc's gloss "ε = sym(displacement gradient)" holds only to first order (at
  second order they differ by rotation terms incl. a |W|²/6 piece the form drops).
  Receipt works purely in symmetric-ε space + R4 conjugation, so receipts are sound —
  fix the gloss to Biot strain, keep the form. μ extraction unaffected (rotation term
  separates additively; deviatoric coefficients verified).
- (b) cutoff-slaving sensitivity (the bigger one): K's ln(ℓ/a) HELD FIXED by stated
  convention — honest labeling ("STATED, not derived," twice ✓). But if ℓ slaves to
  line density (ℓ ~ L₀^{-1/2}), K drifts at O(ε²) WITH the stretch and feeds the
  quadratic form at O(1)-relative (my scratch: factor ~(2−1/2ln) on μ — same order,
  same sign, number moves, possibly ~2×). Stated-convention status ACCEPTED for the
  paper modulus (sign/form robust under either slaving) WITH required priced caveat:
  μ's NUMBER carries O(1) outer-cutoff sensitivity; FB-3's "no unaccounted sink" is
  too strong as written — soften to "sink priced as cutoff sensitivity."
Neither footnote touches μ>0 or the LANE-1 logic; both travel with the formula.

## FB-1 window-half pricing honest (watch-item 2 answered) ✓

τ > 0 by Kelvin-transit scaling with the limit STATED (topology pins Lk not length —
smoothing-shortens-Λ conceded openly); verdict PENDING window; "no static μ(0)"
refuses the overclaim ✓. Pricing, not proof — labeled at every touch. Honest.

## Priced halves not spent (watch-item: verified) ✓

Formula ESTABLISHED / window PRICED / FB-4 trivially-expected-WITH-reverdict-clause /
FB-6 HOLD-conditional / FB-5 conditional-on-window: each priced item labeled, none
promoted. FB-4's ⟨H⟩=0 symmetry-grounded (isotropy ⇒ achiral mean) ✓ plus S4-deviation
tripwire ✓. (Watch-item 3: trivial-expectation + cited-model-property + re-verdict
clause = honest tiering, confirmed.)

## CONTRAST: model SOUND, receipt R6b VACUOUS-as-computed (watch-item 4 answered)

Model argument sound (equilibrated gas at fixed V: W fixed ⇒ μ = 0 — equilibrium
thermodynamics; regime statement honest about fast-deformation limits) ✓. BUT R6b in
code is `W_wave := 0` then assert-zero — asserts its own definition (likewise R6c
checks internal consistency, not discrimination). REQUIRED repair: replace with the
honest receipt — assert dW/dε_dev = 0 given W = W(V)-only (one-line thermodynamic
identity, actually checks content). R6a stands (rules out smuggled K-free terms —
meaningful ✓). Discriminator verdict PENDING receipt repair (model-PASS already).

## LANE-1 verdict per half (as chartered)

F-A formula-half: PASS (in-model + footnotes a/b). Window-half: PENDING (priced).
CONTRAST: model-PASS, receipt-PENDING-repair. LANE-1 overall PENDING — already typed
so in 01 ("neither half alone licenses anything") ✓ no retyping needed.

## F-C charters: CONFIRMED INDEPENDENT (answers routing)

F-C (tilt coupling, own falsifiers) does NOT require LANE-1 (medium story) —
independent construction; charter on its own scope whenever ready. BUT its
missing-5 interpretation inherits LANE-1's pending status (a tilt observable without
a licensed medium is data without a theory home). Precise: build-license YES,
interpretation-license PENDING LANE-1.

## Verdict

FORMULA-HALF PASS with footnotes (a) Biot-strain gloss + (b) cutoff-slaving caveat;
window-half PENDING (honest pricing); CONTRAST model-sound with R6b receipt repair
REQUIRED; LANE-1 PENDING as typed. F-C build-license CONFIRMED independent.
