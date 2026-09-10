# Joined-PR readiness dry run — #198 close + #203 advance (atlas, 2026-09-10)

No PR opened. One verdict per item. Evidence pointers inline.

## READY
- **P251/#198 side**: PR #199 MERGED; harvest repair `de980dc` in history.
  Nothing outstanding except the closing act itself.
- **Comms/ledger**: audit clean (ledger-audit.md), L-1 closed (0116 SOUND),
  protocol v1 live, board renders green, health HEALTHY. PR scaffolding holds.
- **Reviewed bank**: 0108–0122 all drift-reviewed; honest negatives promotable
  as certified records — 0119 (secant NOT the path), 0121 (lemma FAILED
  margin 0.45, δλ≤24.7 vs 11.2). A filed PR should promote these as no-go
  verdicts, not bury them.
- **0111 EM-map J1–J4**: LANDED, drift J1–J4 PASS (review-cipher-jointgating).
  Evaluation gated pre-(a) by design.

## BLOCKED (critical path)
- **R-EM2 gate — BLOCKED-on-owner** (`attempts/0112-cipher-rem2-draft/README.md`,
  DRAFT ONLY; drift draft-technical PASS banked). Without a user
  approve/amend/decline, the PR MUST scope to one-way + conditional framing
  and carry B-R1..B-R5 as open bridges. This is the single owner action that
  unlocks two-way claims.
- **M2-B1 H4-proof + A3 filament — BLOCKED-on-cipher-in-flight**: A1 H5 PASS,
  A2 H1 repaired PASS, A3 transfer PASS all banked; A3 design frozen C1–C8
  (gate+orbit PASS banked, mono m0–6 running) but drift design-vs-code
  pre-verdict flags R-A (floor 1e-4) + R-B (SOFT3 unimplemented) BLOCKING,
  R-C wording (cipher: R-A/B/C applied, deflated rerun in flight, unverified);
  current token H4-proof (+owner-review).
- **G-a2 numbers — BLOCKED-on-(a)** (beacon): 0120 trust DONE (rows-met
  6.6e-3, PASS); 0121 lemma FAILED (honest); 0122 CLOSED as measured-stall
  PASS (drift review-beacon-0122landing: gate reproduced; 1 receipt repair;
  beacon DONE line still undeclared, cosmetic). (a) next rung 0123 spec
  landed 46586bc6 (x18 decision, estimate-side Lipschitz-Weyl, x18_probes.py
  in tree), building. (b)-first APPROVED, proceeding.

## IN PROGRESS (carry, don't block)
- **S9 test** (beacon, `0114-beacon-s9`): dir landed, no verdict yet.
- **IDEA-01..08**: all ADOPTED with owners (cipher/drift/atlas/beacon+cipher);
  IDEA-03 (R9-start) 1/3 closed — gates G-a2 trust-region start, rides with (a).
- **Shadowing scope**: A1/A2 closed, H4-proof current token.

## Claim-promotion backlog
- Nothing awaits the registry: the campaign banks attempts; promotion happens
  *through* the terminal PR. Backlog = drafting the claim delta: which units
  promote (0104 G-blocked/H-established-linear, 0109 core re-verified), which
  negatives certify (0119, 0121), and what stays conditional (J1–J4, trust
  rows-met, S9, A3). No hidden debt found by the ledger audit.

## Bottom line
Two owner-visible actions stand between now and a joinable PR: (1) the R-EM2
ruling (approve/amend/decline — sets the PR's claim scope); (2) beacon landing
0122 + drift landing-verdict (closes the G-a2 number state). Everything else
is in-flight work the PR can carry as named open tracks. Recommended PR shape
if both land: close #198, advance #203 with one-way + trust-region results,
certified negatives, and B-R1..B-R5 bridges explicit.
