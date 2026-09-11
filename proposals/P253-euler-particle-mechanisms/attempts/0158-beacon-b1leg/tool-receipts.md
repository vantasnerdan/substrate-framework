# Tool receipts - 0158 (beacon B1-leg supply)

- Design frozen + committed BEFORE compute (design.md);
  ledger claims 0158 (HIGH-WATER → 0158).
- `b1_leg.py`: exit 0 (19.9 s). Calibration replicates cipher
  B1: 1.0000/2.0003/3.0012, control −0.0000, N-leg spread 0.011
  (within B1 5% bars → supply PASS per frozen bars).
- Supply: advected loops read ≈0 post-crossing (min-traj-dist
  0.0009 → link class 1→0 changed mid-flight, endpoint guards
  insufficient); N-ladder 128/256/512 agrees to ~0.005
  (quadrature sound — readings correct, not floor failure).
- A3 run_a3.py (ring_state/rhs3/rk4_3/constants): READ ONLY.
  Cipher run_ab.py geometry replicated read-only (local copy).
- Ownership: Φ TABLE supplied; B1 leg adjudication is cipher's.
