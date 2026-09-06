# Hodge multipole-order exposing correction

The first oracle correctly computed the octahedral vector, rank-two, and
rank-three orbit projectors but the prose mapped tensor rank to Hodge decay
one order too aggressively. For compact divergence-free vorticity the first
spatial moment is a trace-free antisymmetric rank-two tensor and produces the
`r^-3` velocity dipole. The second spatial moment is rank three, symmetric in
its two spatial indices, and produces the `r^-4` quadrupole. Equations (4)
cancel both. The third spatial moment has rank four; the proper octahedral
group admits a cubic rank-four invariant, so an `r^-5` octupole can survive.

The corrected robust result is therefore `v=O(r^-5)`, `grad v=O(r^-6)`. This
still gives finite zeroth and first lattice moments in three dimensions. The
API rank-two projector was strengthened to accept the physical antisymmetric
vorticity moment, and its test now exposes that trace-free orbit explicitly.

The initial `exact-check.*` and `focused-pytest.*` receipts are retained as an
exposing false green. The corrected repository-interpreter receipts are
`multipole-order-correction.*` and `focused-pytest-repaired.*`; they report
four exact passes and four focused tests, empty stderr, and exit `0`.

Post-correction SHA-256 values before this receipt and validation update:

- construction: `1e5d2c769daf81ea688d1a029fa2f2c46ed5f98b5ad57559b10d1a1062cda372`;
- result: `46f8a0da69e30a01bc46cee0d946f0fe3ad3411deba00eefd84d41a57bf459f7`;
- verifier: `895cbd35955209b8b56bab3560b59809affec613190ff26b8ca3cf4c49369c0f`;
- corrected stdout: `f3a908db95d4fd54c845b0f58223111d596199df74b0f2f568a173f0379e0b07`;
- API: `b666df746bf7bb96a04fa7abfa3e08b230622d2c8b89cc6bf5d59fb843f4fe63`;
- tests: `a3b59c32618fa2e5847a3df0aa7bd81ccec258fed0cf232fbe66c24c880da56a`;
- repaired focused stdout: `50b027a51efd2571ad9c171d98a5ccef20afaedbdcce9b2ed22bed18d22b292e`.
