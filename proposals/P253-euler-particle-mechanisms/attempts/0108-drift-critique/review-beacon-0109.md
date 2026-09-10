# drift firewall review — beacon 0109-beacon-unitg (41d99e99)

Transaction: `verify_unitg_b2.py` (9 predicates) + `witness-status.md` + `b2-edge-transfer.md` +
`tool-receipts.md` (U1–U4). Claimed verdict: algebraic core re-verified; completion BLOCKED at
G-a + G-b; B2 column algebra exact, finite-Cao transfer OPEN. Drift ran independent probes below.

## Verifier soundness (per-check)

- B1/B3/B4/B5: exact identities of the defined frozen column. Constructor confirmed at source
  (`euler_p2_principal.py:130-132` returns `x·[[0,−Z],[R,0]]`), so the checks pin the matrix
  actually consumed downstream — legitimate algebra, honestly fenced by the module docstring
  (no cutoff/IFT/transfer/P2 content).
- B2: closed form vs independent Jordan oracle `(Cn·t).exp()` — genuine two-path cross-check.
  Drift mutation DC-1 (closed form evaluated at t+1): norm diff 0.359 vs 1e−25 tolerance —
  predicate discriminates, non-vacuous.
- B6 (`8·(2Ω)/Ω == 16`): vacuous as physics (tests arithmetic); its real evidence is the 0104-A
  citation. Evidence-role correction (minor, no verdict change): label B6 a constant-regression
  pin against accidental edit, not a derivation.
- G1 (9×9 det on ARBITRARY integer matrices): proves the triangular-block factorization
  det=detH·detC·detG as a formula; says NOTHING about actual H/G nonsingularity on the Cao
  member. Sharpening recorded (not an overclaim — beacon's own G-a holds detH·detG≠0 open).
  Any consumer citing G1 for witness invertibility will be refuted on sight.
- G2: honest import smoke, labeled as such.

## Reproducibility + first-failure legitimacy

- Drift re-ran the verifier: 9/9 PASS, 1.74 s (matches claimed 1.72 s). Exit 0 reproduced.
- U2 cited files verified at source: `route-a-algebra.exit=0`, `route-a-algebra-first.exit=1`,
  `momentum-replay-schema.exit=0`; stdout tail `ALL 10 ROUTE-A ALGEBRA CHECKS PASSED` matches.
  (Stdout filename is `*.stdout.txt`; U2's `*.stdout` shorthand is a cosmetic imprecision.)
- U1 first-failure (B2 oracle test bug, exit 1) + superseded 8/9 run both preserved append-only
  with diagnosis. Legitimate: the "falsifier" wording slightly overstates (it falsified a buggy
  test, not physics) — precise role is sensitivity demonstration, corroborated by DC-1.

## G-a / G-b blocks: named mechanisms, not overclaim

- G-a specifies the exact missing constants (λ/s/R/ε for ω_g and B_g, Q_χ≠0, (25a) circulation
  integral, detH·detG on the fixed charged member), names the data source (actual member fields,
  absent from `src/`), and explains why 0095 never supplied them (profile cells + tube
  asymptotics). Textbook missing-construction block with executable next step.
- G-b (propagation + packetwise same-target bridge) cites 0107§12 prospective status unchanged.
- Activation restraint explicit (0109 adjudicates neither 0107 activation nor §9/§12) — proper.
- Draft-citation hygiene: rows with no new check (packet orthogonality, curve/quantifiers) are
  labeled "cited as 0107 derivation", i.e. unreviewed-author-formalism role, never conclusion.
  No refutation mechanism exhibited; none claimed.

## B2 handoff scope (cipher PoC-1 denominator)

√(R/Z) amplification + t_s→∞ + Z=0 shear handed over as EXACT COLUMN algebra with four hard
limits and a finite-Z normalization instruction — consistent with PoC-1's EXPLORATORY cap.
Firewall condition on the consumer: evaluate at finite PoC edge parameters, never take Z→0 to
tune the ratio; report order-of-magnitude with limits (1)–(4) attached. Transfer itself stays OPEN.

## Verdict

PASS the 0109 verifier and status ledger (modulo the B6-role + G1-sharpening notes above, neither
verdict-changing). Unit G completion: BLOCKED at G-a + G-b (confirmed). B2 finite-Cao transfer:
OPEN (confirmed). No refutation, no narrowing, no particle inference. Next executable: G-a field
data source; B2 denominator now available to cipher PoC-1 under the stated consumer condition.
