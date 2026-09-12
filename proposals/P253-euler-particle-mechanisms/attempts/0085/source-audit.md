# Source and authority audit

## Reviewed framework inputs

- P253/0068 as independently reviewed by P253/0071:
  `attempts/0071/review.md`, SHA-256
  `50e1cb63d6c115a92c3fe6720eaeae0b707c5b81871db663a7a804a9a3127155`,
  and `verdicts.yaml`, SHA-256
  `fd59cc3209460f85aa97cacf66ea5a88d6d5cd9ce8a1630b8e37cf1f31d0f1c0`.
  Consumed scope: the declared constrained `U(1)` action, transported current,
  field energy/momentum balances, constraint kernel, and local flow.  It does
  not supply a charged carrier or coercivity.
- P253/0077 as independently reviewed by P253/0081:
  `attempts/0081/review.md`, SHA-256
  `09f8580b2ff4120ef85046b22d85875ea1371e6b2a028913fe888a28cbcb00ce`,
  and `verdicts.yaml`, SHA-256
  `27ae74404344242b78b1757d7ffb00af307248dfbbc87b6464289cd81f95e099`.
  Consumed scope: subluminal traveling Maxwell elimination, tag locking on a
  regular nonisochronous band, affine Coulomb domains, and `g` parity.  Its
  branch remained conditional there.
- P253/0074 as independently reviewed by P253/0078:
  `attempts/0078/review.md`, SHA-256
  `676fd487d310690cfe89246a5f0a8ba71c676c39c02aa69b4dd2a393a4b2c95e`,
  and `verdicts.yaml`, SHA-256
  `33725bbee5e149ca152146583614a8b8cc02bcad6ec49ed9e24ea77a06c2fe3f`.
  Consumed scope: two distinct simple positive-Krein uncharged Cao modes.
  Their fluid-fiber Riesz isolation is not joint Euler--Maxwell isolation.

P253/0080 is now consumed at its independently reviewed P253/0084 scope:
`0084/review.md` SHA-256
`3fb7d2a69489668f6168d7f57eb9b4e86e35621064dd5f55172ca24270279aac`
and `0084/verdicts.yaml` SHA-256
`7793688a23156d1b015ed82505f5e48b88223e19e8a36563c2f57ef63f0a74c7`.
Its use is the exact finite-window charged branch with fluid/tag/speed change
`O(g^2)` and field size `O(g)`.  The saddle proof is unconditional at `g=0`;
the uniform small-`g` transfer uses only this reviewed branch regularity.

## Primary source

Cao, Lai, Qin, Zhan, and Zou, *Uniqueness and stability of vortex rings for
3D incompressible Euler equation*, arXiv:2206.10165v2.  Cached PDF SHA-256
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`
and text-extraction SHA-256
`aa808fcf435a307c40815d4fdf81d4711c4a2750275e9810fbfc9b501f6085fd`.

The source definitions (1.8)--(1.9) and Theorems 1.2 and 1.8 concern
axisymmetric vorticity rearrangements and the maximization of `E-WP`.  They
support the fluid sign and the axisymmetric-only boundary used in Section 3.
They do not include Maxwell radiation, a material tag, arbitrary full-3D
perturbations, a joint Hessian, or joint spectral stability.

## Direct derivations

The field completion (5), Fourier Gauss-fiber minimization (8)--(12),
regular-band tag identities (21)--(22), compact-curl graph Weyl construction
(22a)--(24), and
radiation-shell criterion (25)--(27) are direct calculations from the
reviewed joint equations.  The exact helper tests rederive their local
algebra.  No oracle is used as evidence for functional closedness, the Cao
high-frequency compactness argument, a limiting-absorption theorem, or a
Fermi-golden-rule sign.

No result from P253/0083, no particle comparison value, and no quantum or
neutrino model is imported.
