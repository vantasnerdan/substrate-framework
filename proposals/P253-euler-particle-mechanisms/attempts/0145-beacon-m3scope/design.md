# 0145-beacon-m3scope - M3 recoil-reciprocity SCOPING (beacon)

Charter: shepherd ruling (M1 kill 9fb212d1; M3 next, M2 stands
down). Scope-first, NO BUILD. Build charters separately.

## Construction

Supplied carrier: trust-r3 tensor state (honest-label CONDITIONAL
per 0067 caveat — carrier at res floor, not root). Two DISTINCT
Gaussian drives (width 0.25, unit impulse): A at core (r=1,z=0.5),
B off-core (r=2,z=1.5). Response per drive: linearized J⁻¹δf
(dense direct, tested pattern) → far-field momentum change ΔP =
∫_{r>3} r·δu dV (outer annulus, inside box). Ratio R = |ΔP| per
unit drive. Prediction (0067 reciprocity): R_A = R_B (a reciprocal
partner answers both drives alike — the charge-like signature).

## Falsifier (frozen pre-build)

KILL iff (a) |R_A − R_B|/mean > 25% (ratio varies with drive — no
reciprocal partner), or (b) either |ΔP| within 10× of assembly
noise (nothing measured). HOLD iff ratios agree ≤ 25% with both
responses well above noise → recoil reciprocity LIVES as
conditional (supplied-carrier label travels; never derived
charge). Gray rule: 15–25% band → third drive (new center),
no verdict.

## Stop rule + routing (frozen)

Miss → M3 dead with numbers. Death MECHANISM routes M2: death by
divergent/incoherent response implicates scattering too (M2 stays
down); death by drive/reciprocity-specific cause (e.g., annulus
choice, drive overlap) with scattering viable → M2 REOPENS
(shepherd routes). Hold → halt + immediate report (new observable
lives; M2 moot).

## Cost + non-scope

Two dense linearized solves on 800 dofs + assembly diagnostics
(~1 min). No member solves, no mesh changes. M2/M4 separate.
No promotion. No import (drives are test sources, not physics).
