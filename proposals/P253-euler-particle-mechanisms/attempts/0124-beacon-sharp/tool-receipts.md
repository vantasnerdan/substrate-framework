# Tool receipts - 0124 (beacon)

- S0 runs: throwaway `python -c` probes (commands in receipts of 0123
  pattern); no new module — permanent code waits for the build it verifies.
- `sharp_dQ.py --mesh fitted` (default state member-fitted.npz):
  exit 0. G1 = 0.6931, margin x15.5; G2' central = 9.1023.
- `sharp_dQ.py --mesh tensor --state .../0117-beacon-member/member-trust-r3.npz`:
  exit 0. Q reproduces feed-trust-r3.log digit-exact (11.15715871 x2,
  35.03562717). G1 = 0.8675, margin x12.9; G2' central = 9.1750.
- S2 P2 probe (inline `-c`, exit 0): P1 rms 0.010329 max 0.100745;
  P2 rms 0.008902 max 0.143701.
- P2/ElementTriP2 + Basis intorder/quadrature API verified present
  (skfem 12.0.2) via import probe, exit 0.
