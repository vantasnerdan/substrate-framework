# Sol-High tag-dimension and coefficient-normalization correction

Date: 2026-09-06

The complete source state in the comparison is `(u,chi)`. The corrected
principal ledger therefore retains the transported material factor:

    chi_(Euler+tag)=(lambda+i U dot k)^3,
    chi_(Maxwell+tag)=(lambda+i U dot k)
                         (lambda^2+c_EM^2 |k|^2)^2.

Equation (9) remains explicitly the velocity-only versus Maxwell-transverse
comparison. The stronger physical-dimension mismatch in (9a) includes the
tag/Gauss row and does not remove the unresolved fluid velocity from the 0068
extended state.

The coefficient ledger also distinguishes displayed inputs from physically
identifiable combinations. Under `(A,phi)->a(A,phi)`, the equivalent
parametrization is `epsilon->epsilon/a^2`, `mu->mu*a^2`, `g->g/a`; hence
`epsilon*mu` and `g^2/epsilon` are invariant. A tag rescaling adds another
normalization convention. No coefficient value is selected by Euler.

Pre/post SHA-256:

| Artifact | Pre | Post |
|---|---|---|
| `derivation.md` | `152a5760c9ea7ee45c05c1b7ff2a9f6cd06a7c6ac9e0661dd40bf660ddac7c63` | `44efb1909bd0a9475c185a7297cd60b5809c50b236a4256ad132b3f75224b525` |
| `result.yaml` | `33fd846aa1f7335ea20e7c82b922c4d745d7d40cbe8a8ef242ddd330fe94080f` | `343d6feaae6c4cae42d1eaba805d72f89c138539089c85f8d2356c52836b79e8` |
| `verify_gauge_emergence.py` | `7173000df7133b575d3d2db842f6c298a2815a8a097b9732754f4f0b3bbabc3b` | `50696a3f0401ae99a9872746f68532af843a9111cb4f1b8e7cd87f20bb19293a` |
| `src/substrate_framework/euler_gauge_emergence.py` | `8f281970d078bbe7888f46c3ef416d96a92d32ee20a70d0e787ec755b84419af` | `3a21e065edf2749ae655700b90fa97d0dbf7ef949e44c4b964d5451e10afc6fd` |
| `tests/test_euler_gauge_emergence.py` | `ded183c43e71d4303bfb61dc7c05cd023286383f273e1d8d31f01e25e0db63f4` | `8f5ec1768397cdcb1c6600195d5de0a813d5d77f2048c9c9b4fc2fa7a506ddaf` |

The repository interpreter rerun is captured in
`tag-dimension-correction.{command.txt,stdout,stderr,exit}` and passes ten
exact assertions with exit zero. The focused API rerun is captured in
`tag-dimension-pytest.{command.txt,stdout,stderr,exit}` and passes three tests.
The roots, no-local-conjugacy verdict, and active structured-background routes
are unchanged.
