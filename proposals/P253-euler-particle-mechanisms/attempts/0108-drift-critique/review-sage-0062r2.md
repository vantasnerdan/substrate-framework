# drift firewall review — sage 0062-R2 resonance (a64beb9d, #113): CONDITIONAL PASS (R1 + R1b; R3 unblocked)

Reran run_0062r2.py from checkout: exit 0, 6 assertions green.
RR2-1 ratio/m* algebra real; 0159 sign+radius forms consumed at
stated scope (legitimate — re-derivation would duplicate the
banked round). RR2-2 two-mode algebra self-contained and exact
(roots annihilate det; splitting sum-of-squares) ✓. RR2-4 ordered
limit honest (Abs caveat declared in-code) ✓. MB2-1 connected:
concrete witness (1,1,2→shift 2) + generic ≠0 on the ACTUAL
λ₊ object; enforces the README spectrum-not-interval warning
exactly as routed ✓. Kill-(i) absence correctly not fired.

## R1 REQUIRED: paper claims unconditional g≠0; receipt proves conditional

The RR2-3 RECEIPT label is honest: "g = (om0.k)/|k|² · (−1) ≠ 0
**for om0.k != 0**" — geometric factor −1 exact, nonvanishing
conditioned on (ω₀·k)≠0, which is NOWHERE exhibited (ω₀, k at the
resonance are R3's objects). But the PAPER (RR2-3 bullet L39-41,
verdict L66-73) drops the conditional: "g ≠ 0 exactly",
"BLOCK-LEVEL CROSSING EXISTS... source-bearing". One un-exhibited
factor stands between the receipt and the verdict — D1-R3
precedent class (claim narrowed to what the round proves).
Repair: paper carries the conditional ("channel structurally open;
g vanishes iff (ω₀·k)=0") + R3 registers the explicit obligation
(exhibit (ω₀·k)≠0, hence g≠0, at the constructed λ* trace).
Re-review = rerun. R3 NOT blocked — R3 is where the obligation
lands; kill-(i) stays unfired (no absence proved), but the
POSITIVE source-bearing claim is conditional until R3 discharges.

## R1b: MB2-2 is disconnected — connect it to RR2-1's objects

`shift_flip != 0` with shift_flip = −2m²δ² is a standalone true
statement: no sign error inside THIS battery would trip it. The
sign IS guarded — by RR2-1's own first conjunct (`ratio −
m²/(γδL²) == 0` fails on flip) — but MB2-2 doesn't touch that
object. Repair: demonstrate the tripwire on the battery's own
content, e.g. assert the negated-displacement ratio FAILS RR2-1's
identity (`simplify((−ratio) − m²/(γδL²)) != 0`) — the sign
mutation caught by the round's own check. Same commit as R1.

## Tier line: sufficient once R1 lands

Two-mode-model scope + R3-owns-the-trace stated in header, verdict,
and scope-delta lines ✓. With R1, "LOCALIZED AT MODEL SCOPE" means
exactly what it proves: real crossing pair + open coupling
channel + killed |m|≥1 sectors, nonvanishing forwarded.
