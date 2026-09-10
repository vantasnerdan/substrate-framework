# Joined-PR readiness dry run — #198 close + #203 advance (atlas, 2026-09-10)

No PR opened. One verdict per item. Evidence pointers inline.

## READY
- **P251/#198 side**: PR #199 MERGED; harvest repair `de980dc` in history.
  Nothing outstanding except the closing act itself.
- **Comms/ledger**: audit clean (ledger-audit.md), L-1 closed (0116 SOUND),
  protocol v1 live, board renders green, health HEALTHY. PR scaffolding holds.
- **Reviewed bank + certified negatives**: through 0125 drift-reviewed (0124 NO-GO CERTIFIED then re-certified stronger with /6 correction; 0125 conditional PASS). Honest negatives promotable as no-go verdicts — 0119 (secant NOT the path), 0121 (lemma FAILED, closed via 0124), 0124 (G2-exact transfer fails with named mechanism; floor ~1e-2 robust). A filed PR promotes these, not buries them.
- **0111 EM-map J1–J4**: LANDED, drift J1–J4 PASS (review-cipher-jointgating).
  Evaluation gated pre-(a) by design.

## BLOCKED (critical path)
- **R-EM2 gate — BLOCKED-on-owner** (`attempts/0112-cipher-rem2-draft/README.md`,
  DRAFT ONLY; drift draft-technical PASS banked). Without a user
  approve/amend/decline, the PR MUST scope to one-way + conditional framing
  and carry B-R1..B-R5 as open bridges. This is the single owner action that
  unlocks two-way claims.
- **M2-B1 H4-proof + A3 filament — BLOCKED-on-cipher-in-flight**: A1 H5 PASS, A2 H1 repaired PASS, A3 transfer PASS all banked; A3 design frozen C1–C8 (gate+orbit PASS banked, mono m0–6 running); drift design-vs-code pre-verdict flagged R-A/R-B BLOCKING + R-C wording; cipher repairs landed 538918bb (1e-4 floor, SOFT3 deflation, D4 addendum) — drift re-verdict pending; current token H4-proof (+owner-review).

## IN PROGRESS (carry, don't block)
- **G-a2 numbers — 0125 PASS CONDITIONAL** (beacon + drift review-beacon-0124correction): conditional lemma rho*~0.03 banked; 0124 RE-CERTIFIED STRONGER with /6 correction; R1–R3 DONE 37625ab1 (stale framing fixed, corrected receipts, sweep banked). G-a2 construction track complete — remaining: R-EM2 scope decision governs what it licenses.
- **IDEA-01..08**: all ADOPTED with owners (cipher/drift/atlas/beacon+cipher);
  IDEA-03 (R9-start) 1/3 closed — gates G-a2 trust-region start, rides with (a).
- **Shadowing scope**: A1/A2 closed, H4-proof current token.

## Claim-promotion backlog
- Drafts EXIST: attempts/0126-beacon-claims/drafts.yaml (C-EUL-001/002, registry schema, `review: draft` — no registry writes, namespace free as of 2026-09-10).
- Path to promotion: drift firewall + individual review first (0123/0124 verdicts + /6 correction cited in-entry), then R-EM2 scope decision sets what the claims may license.
- No hidden debt: ledger audit clean; negatives bank (0119/0121/0124) rides as certified no-go records.

## Bottom line
Two actions stand between now and a joinable PR: (1) the R-EM2 ruling
(approve/amend/decline — sets the PR's claim scope and what C-EUL-001/002 may
license); (2) drift's A3 re-verdict (cipher repairs landed). G-a2 is fully
banked (0125 conditional PASS, 0124 re-certified, R1–R3 done, drafts filed).
Recommended PR shape: close #198, advance #203 with the conditional lemma,
certified negatives (0119/0121/0124), and B-R1..B-R5 bridges explicit.
