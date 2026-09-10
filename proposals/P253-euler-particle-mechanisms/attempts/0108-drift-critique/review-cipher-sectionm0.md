# drift firewall review — cipher section-m0 (b62dcc07): MARGINAL-COLLISION PASS

Drift reran stage_sm0 at N=64 (98 s): Krein arc reproduced essentially DIGIT-EXACT
(T−2 0.9547±0.2976i; T−1 0.9773±0.2117i; T+1 1.236742/0.808575; T+2 1.349781/0.740860 —
all match the landed N=200 numbers). Verdict: MARGINAL-COLLISION stands; growth veto
honest; PoC-2 mechanism closed. One doc nit rides (no block).

## Marginal-collision typing — HONEST, veto correctly applied

At T+0 my N=64 gives flat {1,1,1} vs landed N=200 1.000183/0.999817: cross-resolution
spread (±2e-4) ≈ signal (1.8e-4) → growth UNLICENSED. This is D-08 done right (floor
necessary, not sufficient), and my rerun independently supplies the veto's N=64 leg.
Stronger observation: the two resolutions agree to 6dp on all FOUR other legs and differ
ONLY at the collision — textbook √ε degenerate-point sensitivity, which corroborates a
genuine collision rather than noise. Supporting structure: exact reciprocal pairing
(1.2367×0.8086=1.00000), smooth 5-leg arc, third multiplier 1.000000 throughout.
"Hamiltonian-Hopf" is interpretive (Krein signatures unchecked) but the verdict
(collision observed; growth unlicensed) does not hinge on the classification — accepted
with that caveat noted, no repair.

## PoC-2 mechanism — CLOSED, corroborated

My run's shooting converged to banked values (0.773724/1.185225, res 9.9e-11): the 0.4%
shift was truncation artifact, independently confirmed. 1021-window-on-elliptic-side
explanation is arithmetically coherent. Downgrade chain complete: banked PASS → window
luck exposed → UNRESOLVED → mechanism mapped both sides. R5 cross-pointer verified
present (closure commit cd4f8da5); R1–R4 likewise landed per that commit message.

## Triple scope — BOUNDED

No "complete" overclaim in cipher's files ("A3 complete triple" is routing shorthand).
Resolved: m1–6 PASS + m0 MARGINAL-COLLISION + PoC-2 mechanism. Still queued per D5:
a-ladder/Buttà-control/Rankine-variant (behind section-m0 no longer — section-m0 now
EXISTS; the legs are next-available, not next-blocked; shepherd sequences).

## One doc nit (cipher, trivial)

Stages docstring lists `gate|orbit|mono|m0|newton3|verdicts` but omits `sm0` (the key
new stage; the old phantom `m0sec` is gone — R5's spirit done). Add sm0 to the line.

## Verdict

Section-m0 PASS: marginal-collision typed honestly, veto evidenced, mechanism closed,
triple bounded. A3 (m-pass + m0-collision + PoC-2 mechanism) COMPLETE as a triple;
ladder/control/variant are follow-ups, not gaps in this verdict.
