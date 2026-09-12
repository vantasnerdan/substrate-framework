# 0110-beacon-s3s9 — S3–S9 primary re-verification on the S2 carrier (beacon)

Bounded obligation (shepherd ack 41d99e99): re-verify each frozen-issue
supplier at its primary source and pin exactly what transfers to the S2
fixed thin Cao ring carrier. Surface: `attempts/0110-beacon-s3s9/` only.

## Method

- Cache-first: Cao/Choi/Gavrilov/HF PDFs re-hashed at the 0110 boundary
  (all 4 match recorded SHA-256; see `tool-receipts.md` V1).
- Fresh primary reads this session (arXiv/NYU fetches V2–V8): full PDFs for
  S3, S4, S5, S6, S7-companion, S9; full CLV PDF (S8 second supplier).
- Memory searched first (D1–D2): no direct Dávila/Slobodeanu coverage beyond
  the frozen candidate universe — no duplication, all rows freshly sourced.

## Verdicts (one line each; detail in per-source notes)

- S3 Dávila (v4): finite-window exact leapfrogging for similar positive rings;
  join with S2 carrier open (profile/separation-regime match, no theorem).
- S4 García (v1): all-time periodic pair in translating frame, no-swirl,
  Cantor-set λ; periodic existence ≠ stability; same S2-compatible regime.
- S5 Slobodeanu-2014 (v3): steady-only duality (converse LOCAL); Dirichlet
  term → forced Euler; axisymmetric-no-swirl finite-energy no-go.
- S6 Slobodeanu-2019 (v4): S³ analytic Q=2 solution; R³ transport FAILS
  (item iv) — sphere ≠ Euclidean particle, by the source itself.
- S7 Faddeev–Niemi: stabilization = E₄ Skyrme term + virial E₂=E₄; existence
  numerical (gradient flow), scale from coupling G; nothing from Euler.
- S8 Gavrilov (hash-verified) + CLV: compact existence only; swirl REQUIRED;
  CLV velocity Hölder/weak; stability+interaction separate.
- S9 Choi–Jeong (v2): L¹-small outward perturbation → linear-in-t tail ∀t,
  coexisting with orbital stability; P1 observables must detect it.

No route verdict changed: R1/R2/R3 inventory stands; R2 time-dependent
transfer still unearned (S5 converse is local; S6 transport fails).

## Files

- `s3-davila.md`, `s4-garcia.md`, `s5-s6-slobodeanu.md`,
  `s7-faddeev-niemi.md`, `s8-gavrilov-clv.md`, `s9-choi-jeong.md`
- `s2-carrier-crosscut.md` — joint transfer table onto the S2 carrier.
- `tool-receipts.md` — hashes, fetches, searches with exits.
