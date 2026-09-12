# B1 F-bar freeze (pre-compute; D-08 form) — phase-charge lock

Construction: source = ribbon EDGE curve of framed unknot (core circle R=1 +
Frenet-frame helix, L turns ⇒ linking(edge,core) = L by construction);
A via Biot–Savart quadrature of unit current on edge; test circuit = core
circle (links edge L times); phase Φ = (1/2π)∮A·dl. Control: far circuit
(center offset 5R, phase ≈ 0). L ∈ {1,2,3}. Honesty scope: KINEMATIC lock
(geometry + BS kernel) — supports P2 structure, not dynamical charge;
missing-5 conditionality stands.
- B1-PASS: |Φ − L| ≤ 0.05·max(L,1) all L AND |Φ_control| ≤ 0.05.
- B1-KILL: deviation > 0.25 any leg, or |Φ_control| > 0.25 (spurious phase).
- Gray between → UNRESOLVED. Tolerance 5% = 5× quadrature floor (~1% at
  N=256, verified by N-leg 128 vs 256 inside the run).
Run reports (Φ_1, Φ_2, Φ_3, Φ_ctrl, N-leg) + verdict, nothing else.
