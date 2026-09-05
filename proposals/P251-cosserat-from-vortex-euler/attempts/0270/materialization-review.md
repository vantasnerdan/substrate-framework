# 0270 materialization result

The existing `euler_compact.py` is a jet/coefficient-section algebra and
`euler_joint.py` is a finite observation algebra; neither provides the
requested cylindrical streamfunction/residual or tag normalization, so the
additive module was justified.

Added `src/substrate_framework/euler_compact_ring.py` (SHA-256
`e7d759f7611a501f25f5f2f14f12cc91bb115a54a04733f46f2b48cbea9979e8`) with
the explicit convention `u_r=-psi_z/r`, `u_z=psi_r/r`, `u_phi=I/r`, literal
cylindrical Euler residual including `-u_phi**2/r` and `u_r*u_phi/r`, and
`(j,a,ell_tag_sq)` geometric normalizations.  The API documents that these
are conditional algebra/normalization data, not a finite-R inverse,
existence certificate, or constitutive law.

Added focused tests `tests/test_euler_compact_ring.py` (SHA-256
`741da3bbfa4d4a7b1854773ef282bd76af915637f086b0fd09abc5b4834434c8`).  The
solid-rotation analytic solution has zero residual; omitting pressure
exposes the radial centrifugal residual; `r=0` is rejected; and positive
literal normalization values are checked.  Focused result: `2 passed`, exit
0.  No full suite, central edits, or commits were made.

`route_verdict: established as conditional additive materialization`.
