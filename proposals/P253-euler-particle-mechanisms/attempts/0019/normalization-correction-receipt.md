# Bounded packet-normalization correction receipt

## Trigger and scope

Independent review `P253/0026`, Unit A, identified one bounded omission in the
written WKB-to-PDE bridge: the raw orbit generator
`a_N=N^-1 A exp(iN vartheta)` has order-one vorticity but order-`N^-1`
Biot--Savart velocity.  This receipt records the minimum correction only.  It
does not alter the carrier, cocycle, gain ratio, accessibility argument,
theorem scope, or symbolic oracle.

The corrected construction defines

    qhat_N=q_N/||v[q_N]||_2

and states the geometric-optics residual and exact-semigroup Duhamel error as
relative `O(N^-1)||vhat_N(0)||_2`.  Equivalently, one may rescale the raw
generator by a factor comparable to `N`.  Scalar rescaling preserves the
coadjoint-orbit tangent, axisymmetry, finite energy, translation quotient, and
energy-gain ratio.

## Review-boundary hashes

| Artifact | Before SHA-256 | After SHA-256 |
| --- | --- | --- |
| `derivation.md` | `b7bc7926e75d1eaae1b681b509b7508063403757673b8d14a14c29914591bfba` | `1f56e289e85c0ffd463067d82da2ae05a0b5d9cfca275f0a1956529018f6e3ae` |
| `result.yaml` | `e6fcc43e552d48d07376cfed47e4cb65bfaef2c39c7cd91bbcaed78d3450d757` | `b320c5216095f32f765e99fa1cde4bd61fd71be726b58fe7457926aa1d04acca` |
| `validation.md` | `295237f75d31791fd1fdbabe85c14209729543e5cea2d8aed0538269bf91e79e` | `107423185d6750e7f85482ed7d3961d3f5a136757c335a0e13b928ea1aa433a6` |

## Unchanged oracle provenance

`verify_return.py` does not test packet normalization and no predicate changed.
It was therefore not rerun.  The original receipts remain:

| Artifact | SHA-256 |
| --- | --- |
| `verify_return.py` | `f47f7554f495956b1382da2754b55ee341482871e92e099aeb6fec51ad73c98d` |
| `first-run.command.txt` | `e09341ec8a88f9f48d8bd2021559c8c24dad6fa7b88fe1e651aaa6f19af7da65` |
| `first-run.stdout.txt` | `f83812fc7fad60ea23b417eaff17c1630bfc8e8c0e0532eccc1789b4559c31b5` |
| `first-run.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

No other `0019` artifact was changed in this repair.

## Post-acceptance formatting cleanup

The `0026` bounded correction check accepted the scientific repair and noted
one duplicated phrase only.  The accepted-boundary `derivation.md` hash was
`1f56e289e85c0ffd463067d82da2ae05a0b5d9cfca275f0a1956529018f6e3ae`.
The phrase `squared-norm relative squared-norm error` was changed to
`relative squared-norm error`; no equation, estimate, claim, or oracle
predicate changed.  The post-cleanup `derivation.md` SHA-256 is
`f6e92640ab7440200c1ffa52dc1fc7f366bb9c5fa37e9b634010270116952fd7`.
