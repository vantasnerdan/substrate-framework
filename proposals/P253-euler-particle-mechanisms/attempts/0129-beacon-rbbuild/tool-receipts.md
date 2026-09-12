# Tool receipts - 0129 (beacon)

- `rb_experiment.py` (default scales 1.0/0.5/0.25): exit 0.
  Printed: extension max|D|/unit-du = 0.0888, cpts = 43;
  1.0x shape-err 3.1292 / 0.5x 1.6974 / 0.25x 1.0373; zero detJ drops.
- Bug trail (all in-worktree, no verdict impact): skfem
  interpolator takes (2,M) not (M,2); tensor mesh z in [0,3]
  (upper-half contour only — half-field revolve convention shared
  with all banked Q numbers); tensor mesh ships mixed-wound
  800/1600 (detJ gate made RELATIVE); mesh.p read-only property
  (fresh MeshTri(doflocs(2,N), t) per constructor).
- Debug probes: inline `-c` mesh-shape + contour-range checks.
-Cipher 0128 sketches + 0071 review/verdicts + 0065/0067
  derivations: READ ONLY, nothing touched.
