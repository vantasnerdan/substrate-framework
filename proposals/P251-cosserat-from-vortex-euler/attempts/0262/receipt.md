# 0262 evidence receipt

The actual local Kelvin blocks and exact-in-K simultaneous Gram construction
are independently supported by0264/review.md and its bounded correction
0264/periodic-supplier-review.md, at the explicit periodic-background and
invariant-gap scope. The current-transfer consumer is separate and requires
the actual physical history supplier; the parent remains active.

| Artifact | SHA256 |
| --- | --- |
| local-normalizer.md | 13a00561e486f7d305bd392bd69eb91a3d4d23d8d6471c07f2f9d1b2bef19729 |
| exact-gram.md | dcde8a0fe47aabd145909b0154fdebeff248a7706278531ac6132fb60152517d |
| periodic-supplier.md | 6cb6175df0ea5bc2c4ae0ca7ccb5d224c1f516b4cbe3efa34bd287703ecc69dd |
| current-transfer.md | f33ba31b0dcfbd9e6c97ad0f9c98ad804055cb59a560fc5b84c004abe8d9d45b |

The direct API tests pass: tests.stdout/tests.exit, eight tests, exit0.
The exact finite Gram oracle passes both its initial algebra run and its
API-consuming run (exact-gram.stdout/exit and api-gram.stdout/exit).
The Bloch supplier oracle passes in bloch-supplier.stdout/exit, testing the
nonconstant Beltrami product rule, ordinary-bracket mutation, full Fourier
projector curl and affine force constraints. Analytic uniform estimates
are proved in the source and reviewed; no symbolic tally substitutes for
them. Targeted Ruff passes. Fixed repository validation passes in
fixed-validation.stdout/exit:270 accepted claims,1063 valid memory records,
43 pre-existing memory warnings. No accepted claim or release was changed.

All scripts ran in separate Herdr shell w3:p2 with the working Python at
/home/dan/substrate-framework/.venv/bin/python and PYTHONPATH=src.
New code is an additive conditional algebra API; prior APIs are unchanged.
