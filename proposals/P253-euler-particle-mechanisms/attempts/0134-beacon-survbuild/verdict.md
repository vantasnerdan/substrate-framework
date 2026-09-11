# 0134 verdict: ALL THREE KILL → filament-charge closure COMPLETE (beacon)

Build: run_surv.py under frozen 0133 falsifiers. Gates first (all
PASS). No HOLD → program ran to completion (no early halt).

## Measured (printed, exit 0 all)

- A1 m≥2 seed: gate 2.62e-4 PASS; content OPEN (m=2 norm 1.18e-1,
  real content — seed worked); κ=+3.58, r=0.977 → KILL (residual).
  Positive sign, ~2% explained: m=2 non-overlap signature.
- A2 single ring: gate 2.78e-7 PASS; content OPEN (4.07e-5, weak
  but linear); κ=+2.04, r=0.501 → KILL (2× above line, outside
  gray zone). Note: κ ≈ 2.04 reproduces S1's κ_B = 2.0406 to 4
  digits with killing residuals both times — scale-coincidence
  pattern corroborated (slopes agree, templates explain nothing).
- A3 λ-family: gate PASS; r flat 0.9026–0.9079 across ALL λ
  (self-participation moves nothing — extra parameter buys
  nothing); best λ=0.25 interior but κ=−6.82 → KILL (wrong-way).

- A2 single ring: gate 2.78e-7 PASS; content OPEN (4.07e-5, weak
  but linear); κ=+2.04, r=0.501 → KILL (2× above kill line, 5×
  from HOLD — boundary closeness recorded, not comfortable).
Lane live threads: S3-residue direction + S4-noise-MODEL as
diagnostic tool (S4-charge dead — no contradiction with banked
verdict). Reopen needs a materially different charge construction
(none on the table).
