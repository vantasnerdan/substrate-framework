# Tool receipts - 0132 (beacon)

- `run_3d.py gate`: exit 0. eps-linearity rel diff halves 2.63e-4
  (PASS < 1e-3). (~43 s.)
- `run_3d.py fit`: exit 0. kappa* = -6.6596; m=1 r=0.9047 KILL;
  m=2 r=1.1158 (inconclusive leg). (~44 s.)
- Per-m content probe (inline, exit 0): response norms m=1..8 =
  2.05e-2 / 3.91e-5 / 1.23e-7 / ... (m=1 real, m=2 noise-level).
- A3 code (run_a3.py) + cipher run_kappa.py: READ ONLY.
- Bug trail: tangents() axis convention for (N,3) vs (2,N,3)
  inputs (fixed to axis=-2/-1).
