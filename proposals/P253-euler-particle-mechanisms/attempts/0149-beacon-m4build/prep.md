# 0149 M4 build prep: ENTRY GATE shown (beacon)

Charter: shepherd M4 BUILD (scope PASS eb77ba2a). Entry gate:
drift's <15% pin satisfied pre-build — shown here, verified at
landing review.

## Pin shown (0148 design.md v2, committed pre-build)

v1: HOLD ≤ 15% OVERLAPPED gray 15–25%. v2 (current): HOLD < 15% /
GRAY 15–25% / KILL > 25%, disjoint. Type-guard (2% tilt) ordered
FIRST as gate (fires before scaling). Supplied-background flag
travels.

## Build (frozen remainder)

Tilted single ring (0.1 rad, parked companion R = 50) under BS +
U_bg = Ω_bg ẑ×x, Ω_bg ∈ {0, 0.1, 0.2, 0.4}. Normal tracking via
per-frame SVD. Guard → rates → linear fit → frozen adjudication.
