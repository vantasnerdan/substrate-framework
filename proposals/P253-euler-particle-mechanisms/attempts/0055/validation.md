# Validation

The exact oracle is `verify_scale_causality.py`.  It derives the similarity
weights, pressure Hessian contraction, axisymmetric quadrupole and shell-sum
boundary from symbolic expressions.  Its checks would fail under a reversed
Euler pressure sign, a Frobenius rather than crossed index contraction, an
incorrect angular-momentum power, a topology rescaling, or promotion of the
borderline `r^-3` lattice tail to an absolutely summable one.

The focused API tests exercise the same public definitions plus invalid-scale
domains.  Neither oracle tests a carrier spectrum or a foundational
extension; the result types both as continuations rather than inferred
physics.  The pressure theorem is a compact-source far-field and initial-time
statement.  It is sufficient to expose the absence of a strict full-field
domain of dependence and does not claim uniform long-time asymptotics.

Independent scientific review is still required before any claim promotion.

The first repository-interpreter execution completed with exit `0`.  The
initial exact oracle printed five passes, and the focused public-API run passed seven
tests.  Captured commands, stdout, stderr, and exit files are stored beside
this document.  At this checkpoint the principal SHA-256 hashes are:

- initial `verify_scale_causality.py`: `7a0917ca88468b254dc8b591708a31384598592c43fb898a243b05f20d814d2f`;
- exact stdout: `c0b917e062f32d5e9bbaf92204c83b9c3ca4ae730c566266bfa33403d187ef4f`;
- `euler_scale_causality.py`: `b15f59edfc3d951d031d4927e50e210e321d6470ccb419c21158c6ec02e280c8`;
- `test_euler_scale_causality.py`: `afa4a16fbaa9f6da505efd126519e220170e985b2813ec08a43e1ff9369597f8`;
- focused pytest stdout: `e6eae79f2ecf62b9d87a811f2f0ee6d9af22657a9582b085e2318955b72f2e3c`.

The extension-ledger continuation added the exact classical-action prefactor
test and the first-spatial-moment lattice boundary.  Its updated verifier is
`51c0d4ae72c49d4eba525f9d55ac1a535d998dac7d34878aaabf0d976212f5b0`;
the repository-interpreter replay printed six passes with stdout
`ea56f577a3df62732bcb84c6e363a91d3f890b61851084313290a218916ccbe2`,
empty stderr, and exit `0`.  The initial receipt remains preserved rather than
being relabeled as evidence for the added assertion.
