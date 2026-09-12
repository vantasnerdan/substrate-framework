# drift fence-consistency audit — pre-PR surfaces (goal-restated proactive)

Audited SYNTHESIS.md (A7/D-sections + F discipline) vs
JOINED-PR-ASSESSMENT-DRAFT.md vs PR-READINESS.md vs underlying
verdicts (drift reviews + attempt bytes). No fence contradiction
found. Findings, severest first:

## F3 (medium): BOARD aims at a stale readiness surface

BOARD.md:20 points atlas's COMMS row at PR-READINESS.md, frozen
2026-09-10. That doc states the FB-director lane "OPENED... P3-B-dyn
owed-with-D4" — CLOSED since #106 — and reviews "through 0125".
Posture ("NO PR until native solution") still matches the draft
assessment ("No PR action authorized"), so nothing false is
currently licensed — but the pointer aims a future PR assembler at
dead state. Repair: BOARD row → JOINED-PR-ASSESSMENT-DRAFT.md;
PR-READINESS.md marked SUPERSEDED (or refreshed) before any PR
assembly. Owner/shepherd call, pre-PR required.

## F4 (medium-low): tilt + pair lanes unfolded in SYN

The 0157 tilt lane (T1 screw-silence #101, edge-T1′ J-firing #102
with ALIVE-transfer fence, T2 back-reaction bound #103) and the
pair lane (D-D1 #107; D-D2 c1a34784 banked) appear NOWHERE in SYN —
fences live only in drift reviews + draft-assessment L9. The joined
PR inherits SYN's buckets; unfolded lanes either ride silently or
get dropped. Repair: SYN tilt/pair fold (or explicit deferral with
owner sign) before PR. D-D2's drift-review status is additionally
unclear from my queue — shepherd to confirm routing.

## F1 (low): A7(ii) verb overstates the gap

"dipole moment analytic mu.b²/2pi reproduced numerically (dE/dd
0.01881 vs 0.01592)" — the gap is ~18%; the underlying verdict
shows the pair with NO verb. Numbers visible (not hidden), but
"reproduced" in a PR-facing synthesis implies closer. Repair:
"agrees to ~18% (sign + order exact; lattice-vs-continuum)".

## F2 (low): A7(i) shorthand is misreadable

"perturbed defect restores E to 5dp" — receipted meaning is
Efin = E_clean to 5dp (NOT restoration of perturbed E0). Repair:
append "to E_clean".

## F5 (info, no action): 0062-R1 repair in-tree, pending routing

2569867b exists with re-type + MB1-3 + wording fix (goes beyond the
ruling: new biting mutation). SYN 11:10Z citation accurate
(fcaa6543 = the ruling, not the repair). Awaiting drift re-review
on shepherd routing — SYN correctly does not claim acceptance.

## Verified consistent (no action)

- LANE-1 gate: "consumed" (SYN D2) + "governance flag, drift
  neutral" (assessment L24) — recorded both sides, unresolved both
  sides, no drift.
- "Family complete" (D4) vs "no lane-level conjunction fires" (F):
  compatible — built-everywhere ≠ claimed-together.
- Assessment posture (success NOT MET; ALIVE scoped survival-tier;
  no electron/carrier/measurement anywhere) matches every fence I
  banked. The draft assessment is the accurate live surface.
