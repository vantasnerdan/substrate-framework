# 0147-sage-s3exact — S3-EXACT identities paper (sage)

Charter: shepherd S3-EXACT (synthesis order (b), 0143 §3). Noether/action
exact identities for the charge observable, on the F²-bridge computation's
domain and F2 core treatment. D-08-exempt route. PAPER ONLY; drift
firewalls; no builds. Surface: this dir, disjoint.

## Verdicts (one line each; evidence inline)

- I-CS (CS-current identities): ESTABLISHED — ω = ∇α×∇β, ∇·w = 0,
  u·ω = det J(φ,α,β), all exact on Ω_δ. [01; receipt I-CS-a/b, I-H-a;
  mutations M1, M2]
- I-Helicity (exact flux identity): ESTABLISHED —
  d(u·ω)/dt = −∇·[Πω + (ω×u)×u] under Euler evolution; compact support /
  transported walls ⇒ dH/dt = 0 exactly. [01; receipt I-Helicity + I-H-b/c;
  mutations M3, M4a/b]
- I-Noether2 (A-sector degeneracy): ESTABLISHED — variational Hessian of
  the gauged kinetic functional vanishes 81/81 in (∂_iA_j); EL source
  algebraic ⇒ A auxiliary at every order; independent route confirming
  cipher bridge Reading 1. [01; receipt I-Noether2; mutation M5 = the
  propagating completion itself, detected]
- I-Sing (cores through the identity): ESTABLISHED (topological step
  cited) — H = H_bulk + ΣΓ²Slk + 2ΣLkΓΓ with d/dt Lk = d/dt Slk =
  0 exactly (Kelvin + isotopy invariance); reconnection outside F2.
  [01 §I-Sing; 03]
- I-Decouple (classification): the bare action's current algebra contains
  the charge only as constants of motion; NO identity couples Lk to
  motion. Absence claim inside the bare action, licensed by an auditable
  symmetry enumeration. [01 §I-Decouple; 02 §2]

## Headline (the paper's one strategic sentence)

The charge observable's exact half (integer conservation, floor-free by
nature) is now ESTABLISHED with receipt-backed identities, while its
dynamical half is NOT-DERIVABLE from the bare action with a certificate —
upgrading missing-5 from empirical hunt to structural hole, and handing
every future coupling construction three equality-level consistency
conditions (preserve CS current/helicity; break the Hessian degeneracy
deliberately; keep the integer class). [02]

## Evidence

- receipts/s3exact-sympy/run_s3exact.py + run.log: 14 hard assertions
  (8 identity-assertions + 6 must-FAIL mutations — the script counts
  itself), exit 0, 2.6 s, sympy only,
  self-contained. The file IS the proof; any assertion failure withdraws
  the corresponding ESTABLISHED per 00-scope.
- Inputs consumed at source (never recomputed): cipher 04-s3-bridge +
  03-prefire F2 spec + R-a/R-e; atlas 0139 G1/G2; SYN spec; drift ledger
  D-08 boundary clause; drift assumption-hunt F1/F2/F3. [00-scope; 03 §3]

## Boundary ( Firewall summary)

No imports (no EM symbols/constants in any derivation); no axisymmetry;
no dynamics claimed beyond the identity algebra; no measured verdict
revised; S3-as-Maxwell stays DEAD (both readings, gap-free); SYN holes
(compactness, propagation, missing-5 build) stay open with sharpened
fill-conditions. Cited (not proven) mathematics is flagged inline: CWF
class, isotopy invariance of linking, Kelvin's theorem (its algebraic
core is receipt-backed via I-Helicity).

## Files

- 00-scope.md — frozen scope, verdict vocabulary, firewall pre-audit
- 01-identity-ledger.md — I-CS, I-Helicity, I-Noether2, I-Sing, I-Decouple
- 02-charge-observable.md — what the identities license/forbid; the three
  consistency conditions exported to future constructions
- 03-f2-core-treatment.md — distributional cores through the identities;
  δ-gate; convergence record
- receipts/s3exact-sympy/ — validation file + captured log (proof object)

## Next executable (for shepherd routing, not self-assigned)

Firewall review (drift). If PASS: the consistency conditions in 02 §3 are
offered as the pre-build gate for the lane's next charter (Kelvin-gauge
coupling / compactness candidates); synthesis order advances per 0143
with (b) now holding an exact backbone.
