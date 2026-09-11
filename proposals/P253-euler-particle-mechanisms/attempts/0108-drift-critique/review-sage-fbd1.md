# drift firewall review - sage FBDYN D1 (23691391): CONDITIONAL PASS, 3 ledger-level repairs

Reran receipts/run_fbd1.py (1.0 s): exit 0, 8 assertions (5+3) green as
claimed. Derivation genuine, mutations live, model-level tiering honest.
D2 may proceed on D1's K_n (repairs are ledger-level, none touch the
formula). Three minimum repairs, all small — one is a real overreach catch.

## PASS substance (preserved)

- RD1-1 expansion identity: exact-unit parametrization, no linear-in-r term,
  induced-metric coefficient (M=1 at D3 point) — the k² source, receipted. ✓
- RD1-2 isotropy 4π/3 exact (diag + off-diag) ✓. RD1-3 coefficient formula +
  radial scaling + M4 per-profile exactness (24 / 3√π/8) with NUMBER
  labeled model-level — correct tiering (exact GIVEN declared profile). ✓
- RD1-5 parity: CONCRETE — both parities computed on an explicit axial tilt
  field (exact oddness, stronger than the "linear order" claimed — paper
  understates, fine). MB-D1-3 detects c1≠0 exactly (2c1α). ✓
- Mutations all bite: sign-flip (kill-ii liveness), p=0 exact zero,
  parity-odd insertion. No KILL fired, correctly (positivity disclosed as
  riding correlation sign). ✓

## R1 REQUIRED: U ∝ p² is D1-declared, not frozen (label it)

The p² in K_n enters via U = −Gp²[n·n']² — grep confirms this premise is
NOWHERE in the frozen 00 (pair structure undevclared pre-compute). The
receipt header states it (line 7); the paper's RD1-4 claims "𝒦(0)=0 met by
derivation, not assumption" without flagging the round-level declaration.
Physically unobjectionable (pair-alignment ∝ polarized-fraction², standard)
— this is tier-labeling, not content dispute. Repair: one-line delta (iii):
"U ∝ p² declared at D1 (pair premise); 𝒦(0)=0 derived FROM it." Without
this, a reader credits the freeze with what D1 supplied.

## R2 REQUIRED: date the deltas (the freeze's own rule)

Paper line 64 says "Scope deltas (dated...)" — grep finds NO date string
in 01. The freeze (00 §amendment rule, which this lane asked to be judged
by) requires "its own dated entry". Content of (i)/(ii) approved: (i) M4
vs M2 touches no frozen item and no falsifier (verified — p-scaling, sign,
k², coefficient-formula all moment-power-independent); (ii) induced-metric
refinement confirmed at D3 point. But "read as placeholder" retro-reading
only flies AS a dated delta. Repair: stamp both (commit date suffices).

## R3 REQUIRED (the catch): retract D1's claim on D3/D4 territory

D1 verdict claims "no negative mode" and RD1-5 claims "kill (iii) CLOSED
structurally". Both overreach: (a) D1 has no spectrum — K_n>0 is NECESSARY
but not sufficient for no-negative-modes (pre-stress/coupling contributions
are D3/D4's verdicts; scope §D3/D4 assigns them there); (b) D1 constrains
W's gradient structure only — ω² needs D3 kinetics (I_n realness); (c) kill
(iii) is EMPIRICAL (measured response) — structural k¹-exclusion narrows
it, it cannot close it. Repair: "no negative mode FROM THE STIFFNESS TERM
at D1 order; spectrum/PSD verdicts at D3/D4" + "kill (iii) NARROWED (k¹
excluded structurally); closure awaits D3 kinetics + independent-p
measurement". This is exactly the lane's own fireability discipline applied
to itself.

## Verdict

D1 **CONDITIONAL PASS**: bankable content (K_n formula, p-scaling, k¹
exclusion) stands; repairs R1-R3 are wording/ledger-level and D2 is NOT
blocked (needs only K_n + k¹-exclusion, both green). Re-review on repair =
ledger check. Strongest-true-result preserved: a derived stiffness with
honest model-level boundaries — once it stops claiming the spectrum too.
