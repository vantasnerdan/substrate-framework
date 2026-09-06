# Validation and limits

The exact oracle enumerates all 24 proper signed-permutation rotations rather
than hard-coding the averaged tensors. It checks determinant and orthogonality,
derives vector and generic rank-two orbit sums, constructs a generic
symmetric-first-pair rank-three tensor and verifies its full orbit sum is zero,
and evaluates both three-dimensional lattice comparison integrals.

The public API tests repeat the group/cardinality and moment calculations and
exercise invalid domains. These checks expose an omitted rotation, use of
improper signed weights, a surviving low anisotropic moment, or confusion of
absolute with first-moment summability. The analytic text supplies the
stationary Euler virial and the distinction between compact static pressure
and a global Hodge tangent. The scripts do not prove an invariant Euler mode,
Riesz band, nonlinear persistence, or particle interpretation.

Independent scientific review is required before consuming the established
moment theorem. No production numerics are used.

The first repository-interpreter execution completed with four exact passes,
empty stderr, and exit `0`. The focused public-API run passed four tests in
`1.92s`, also with empty stderr and exit `0`. Captured commands and outputs
are stored beside this document. Principal SHA-256 values before this
validation append are:

- construction: `4169530d5c5e52d600e04ac6ab5fd787d82eefa1795fa62bcf0baeee23372d7a`;
- result: `6d382dcff6b79c6264f5e48e0505096f25701ca4f7e0d9d1ae55bed2df869f3b`;
- verifier: `e97c31498b324148e2476cd6db277630f25d4ca46ab52c23b55370379a3241fa`;
- exact stdout: `eeb32162de01dbd947deb12e83c177b7e1a2764e174e9dee0aa7db0e2168bfbe`;
- API: `8f8be82e3d5f956851b04b348ce7d0eb423a3d3b50970d607ed29490321ed77a`;
- tests: `288a4cd3f9c51b68e97d3bd720f6b6ec9f521aa5947feee677c048a0db521280`;
- focused stdout: `b4a5d3759e611cd9dde20a8804ee3e75fd0949886be9c311379bcfc57e88e536`.

The tensor-to-Hodge decay interpretation in that first run was then exposed
as one order too strong. `multipole-order-correction.md` preserves the false
green and derives the corrected `r^-5` velocity order. The changed verifier
and strengthened antisymmetric-rank-two API test were rerun with the repository
interpreter; four exact assertions and four focused tests pass, both exit `0`.
