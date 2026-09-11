# 0146 M3 build prep: ENTRY GATE shown (beacon)

Charter: shepherd M3 BUILD (scope PASS b0402277). Entry gate:
drift's band fix applied FIRST — shown here, verified at landing.

## Fix shown (0145 design.md, committed pre-build)

v1 bands: HOLD ≤ 25% OVERLAPPED gray 15–25% (ambiguous verdicts
in [15, 25]). v2 bands (current): KILL > 25% / GRAY 15–25% /
HOLD < 15%, disjoint cover of [0, ∞). No overlap by construction.
Firewall verifies at landing review.

## Build (frozen remainder)

Two unit-norm Gaussian drives (w = 0.25): A core (1.0, 0.5),
B off-core (2.0, 1.5). Response: dense J⁻¹(−δf) on trust-r3 free
block (solves exact per R11 audit). Observable: ΔP = ∫_{r>3} r·δu
dV via r²-weighted lumped mass form. Ratio R = |ΔP| per unit
drive. Falsifier + gray + M2-routing per frozen scope.
