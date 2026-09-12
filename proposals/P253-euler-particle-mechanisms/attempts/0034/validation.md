# 0034 validation and claim boundary

The strongest oracle is the direct Euler variational calculation in
`canonical-construction.md`. It derives the regular-strip Clebsch equations,
the finite-domain impulse surface term, the actual solitary-wave end-jump,
the full nonlocal linearized operator, and the material Routh reduction. Six
focused API tests corroborate the signs and factors; they do not establish the
still-open threshold propagation or nonlinear persistence theorem.

The final focused invocation in the dedicated Herdr scripts pane was

    /home/dan/substrate-framework/.venv/bin/python -m pytest -q \
      tests/test_euler_axisymmetric_action.py

and passed six tests in 3.44 seconds with exit zero. Earlier three- and
five-test receipts are retained as append-only construction history. No
production numerics, soft eigenvalue, force, or small energy splitting was
computed, so no numerical floor is used to support the result.

The Varholm--Wahlén--Walsh transfer was checked against
arXiv:1811.08024v2. Its one-negative-direction Hessian hypothesis fails for
the actual accessible traveling form, whose regular-interior high-frequency
directions have both signs. The source is therefore not used to claim
stability. The next scientific construction is the fixed-Casimir
zero-frequency threshold problem for the full operator, followed by the
nonaxisymmetric sectors.

Independent `0036` establishes the result after the bounded whole-space
relative-Routhian correction recorded in `review-correction-0036.md`. The
unchanged API and source predicates were not rerun. No accepted registry
entry, release, generated record, electron/neutrino identification, or parent
status changes through this result.

## Pinned author artifacts

- `README.md`: `552a165559dea55526b3fb49f67cf3bd4a752d3eda24eaa218ae343d5983b500`
- `canonical-construction.md`: pre-review `5a72e4e5b7bf7a83d8d5add8b0f2e87f2d6bf961f4c9a57764f40137b259ce12`; post-correction `6f233ee844fb3d61c63efb4b7281ee04eaf2536f9753d3958ceba544fe8b3512`
- `review-correction-0036.md`: `16a65d783d01a84121b0303cac27b75d45f3d7f12d97adc7978485196a1fc99d`
- `source-transfer.md`: `dc20c11e0aaceed08f1939a1c210ee348a624944cba6b654315ed8416d54b74d`
- `coordination-model-repair.md`: `18062a99050d919897398b690e43e2893441a037e51a46fe4512027c5e2d9511`
- `result.yaml`: `9f1b1df7fca625a7a5698d1f971fda477a97719aa4050670be7dd1d1df5235e6`
- `jump-api.stdout`: `b3d09481c8a1a1f5f642953336e9b973edf6c09ba57f8a583edcde90501117b4`
- `jump-api.exit`: `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`
- `src/substrate_framework/euler_axisymmetric_action.py`: `e062c38657bfbc911b1174eea2bce77f27d6801e560aa6814cb3aac0b6db181b`
- `tests/test_euler_axisymmetric_action.py`: `93ca6f21eea7554b84414a51b533ef3470c272879934a964194488e4cb97fc11`
