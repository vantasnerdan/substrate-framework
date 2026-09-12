# Validation and evidence boundary

The exact scientific inputs are the individually reviewed `0001/0004`
material balances and pressure-quadrupole counterexample, and the separately
reviewed `0009/0011` Euler split, energy exchange and solved shear memory.
`construction.md` rederives their common state map and the nonlinear
conditional unresolved solution rather than treating review prose as proof.

The additive API exposes the same-field pressure Poisson source and the exact
bounded linear block-elimination signs.  Its focused run used the repository
interpreter on:

- `tests/test_euler_retained_memory.py`;
- `tests/test_euler_material_balance.py`;
- `tests/test_euler_shear_memory.py`.

Result: `12 passed in 5.06s`, exit zero.  The new pressure mutation distinguishes
`trace(Du*Du)` from the Frobenius contraction, and the two-block test derives
the Markov, unresolved-initial and memory signs.  The existing material and
shear tests are affected-consumer checks.  They do not prove an
infinite-dimensional orthogonal-dynamics semigroup.

Pinned pre-review hashes before final result files:

- corrected README: `9c37f046edde581edd3e03e0c69fadaffb5d5ef91c7b1001c87bc1797c013a87`;
- construction: `6e9993547fea6617150e4a7f8f72cefa75dc443026b4215f1e173e2767d1d277`;
- API: `75ab59d079de41481f1bb85d6f9a97a7d230ab00d3458236e43456f4f9362f0d`;
- tests: `46159071e1c7b90952643c2f24efb626ec3e2e26e63b435abea1f0dddf422e92`;
- focused stdout: `4095c3a2c9393f1fc49c7273a9083197179e91f6bb23de146808ec4b1bdcf57c`;
- focused exit: `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

No production numerics or small-ratio quantity is involved.  Full repository
validation and independent synthesis review remain separate checkpoint tasks.

The single bounded `0045` review correction subsequently tightened only the
analytic domain: finite-moment material tags, the finite-rank smooth-range
Sobolev projection, forced-Euler energy-method construction and `Q`-range
preservation, a distinct observable projection `P_obs`, and conditional
integrated ambient momentum. The API and its exact finite-dimensional algebra
did not change, so the focused oracle was not rerun. Current hashes and the
full before/after map are pinned in `review-correction-0045.md`.
