# Tool receipts — defect-dipole build (beacon, 0157/dipole-run, FIRED-ON-#87)

- Envelope: post-f2-build-design.md + #87 fence header (amended,
  committed with this build).
- `dipole-run/defect_dipole.py`: exit 0. Printed: B1 0.36731 /
  0.074% / 0.079%; analytic 0.36647; E(d) monotone 0.0929–0.5067;
  attraction numeric/analytic; E_ann/L 0.0929; tail -2.000; B2
  economics 0.3665 vs 0.1832 + license audit.
- Two banked-honest mid-build failures: (1) B1 fired at coarse
  h (core-edge sampling) — retuned to matched-h pairs, holds;
  (2) B3 probed a nodal line (exz=0 on-axis by symmetry) —
  diagonal probe, exponent exact. Both repairs in-tree, history
  stands per ruling.
- Banked inputs read-only: 0154 02-fa-build (W, mu, lam), 0151
  skirt/core prior, 0117 L~1.9 scale map, 0147 I-Sing license.
  No banked-file edits.
