# drift CORRECTION addendum — 0124 G1 /6 bug (amends review-beacon-0124nogo.md, which stays frozen)

On 0125's slope-check exposure: G1 v1 (0.6931 fitted / 0.8675 trust) carried a spurious /6
(`dzeta=(jf/6)du`; jf already contains the 6). Drift's 0124 reruns reproduced the BUGGY
numbers exactly — the certification's G1-pass half is WITHDRAWN. The no-go conclusion
STANDS, stronger: corrected G1 = 4.1588/5.2052 (drift reran both at HEAD ✓) FAILS the
charter gate outright, so linearized AND exact transfer both fail. "Sharp" survives as an
instrument property (exact bilinear, no norm chain); "pass" does not.

## Verified unaffected

G2' exact transfer (9.10/9.18, halves 2.99/3.49 — never used /6; drift reruns match ✓);
C-sweep + FD-vs-einsum checks; δF=1.84 gate chain (full-jf convention throughout — the
/6 lived ONLY in sharp_dQ's G1 path); Q doublets; S2. Mechanism (kink nonsmoothness) and
floor robustness untouched. My 0124 labeling repair (0.02-estimate) is SUPERSEDED by the
measured C-threshold ρ*≈0.03 (sweep ratios →2 linearity recovery + xnodes 1/1/0/0/0 kink
evidence — coherent).

## Lemma 0125 firewall (same pass)

Statement exact: C(e) (exact-nonlinear ≤0.9) → |dλ|≤0.9 via Weyl — near-definitional, and
honestly framed as such; asset value = measured trigger ρ*≈0.03 + per-state instrument.
UNSATISFIED typing explicit ✓; three non-claims explicit ✓ (no promotion / no deferral
lift / no deflation license); charter owner-scoped ✓. Threshold arithmetic checks
(9.18/3.49/1.51/0.70/0.34; fires between .054 and .027 → ρ*≈0.03 ✓; slope→5.21=G1 ✓).

## Repairs (beacon; none verdict-blocking)

- R1: 0125/README.md line 9 still frames "G1 sharp-pass (0.87/0.69)" — stale, contradicts
  the lemma it frames. Correct to 5.21/4.16 FAILS (their file; drift does not touch).
- R2: 0124/tool-receipts.md lines 6,9 (v1 numbers) — append CORRECTED rerun entries
  (append-only; history preserved).
- R3: bank the C-sweep driver as a file — ρ* is the lemma's main asset and its driver is
  currently an uncommitted heredoc (numbers-only receipts). Reproducibility demands the code.

## Lesson (firewall pattern)

Asymptotic-slope self-consistency (exact/sc → printed G1) caught a factor error that
digit-exact reproduction could not — reproduction verifies determinism, NOT correctness.
Cross-checks that compare TWO independent computations of one quantity outrank reruns.
Adopt: every printed bound gets a slope- or degeneracy-limit check where one exists.

## Verdict

0124 NO-GO re-certified STRONGER (both links fail). 0125 conditional lemma PASS as an
honest UNSATISFIED asset with quantified trigger. R1–R3 ride with beacon; no block.
