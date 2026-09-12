# Tool receipts — 0109 first executions (beacon)

## U1 — 0109 verifier (2026-09-10)

Command: `PYTHONPATH=src .venv/bin/python proposals/P253-euler-particle-mechanisms/attempts/0109-beacon-unitg/verify_unitg_b2.py`
First run: EXIT 1 — `AssertionError: B2 frozen-column exp closed form`
(diagnosed same session: test bug, `Cn.exp()` = exp(C·1) vs closed form at
t = 0.7; oracle fixed to `(Cn·t).exp()`; failure preserved here append-only).
Rerun after fix: EXIT 0 —
`PASS B1 column square / B5 order-zero symmetrizer / B2 frozen-column exp
closed form / B3 edge nilpotent square / B3 edge shear entry /
B4 edge amplification sqrt(R/Z) / B6 center K equals 16 /
G1 staged block determinant / G2 Cao jet import smoke /
ALL 9 UNITG-B2 CHECKS PASSED` (1.72 s).
Intermediate rerun (8/9, B3 checks lost to a bad edit range, restored same
session) is superseded by the final 9/9 at unchanged code semantics; the
initial B2 failure receipt above is the preserved falsifier run.

## U2 — 0107 Route-A algebra receipt (cited, not rerun)

`attempts/0107/route-a-algebra.exit` = 0; stdout tail (10 lines PASS +
`ALL 10 ROUTE-A ALGEBRA CHECKS PASSED`); `route-a-algebra-first.exit` = 1
(preserved local-Maxwell-oracle first failure); `momentum-replay-schema.exit`
= 0 (J_ren execution license). Recovered via background job bg_18 this
session; 0107 README + verifier source re-read at 0109 boundary
(`verify_route_a_algebra.py#A80D`, 146 lines: Maxwell residual identity,
div-curl row, impulse factor/sign, moment reconstruction, 3 center
potentials, scalar M_leaf det, K-center 16).

## U3 — framework API surface (2026-09-10)

`grep "def " src/substrate_framework/euler_p2_principal.py` →
`maxwell_transverse_block, column_bas_matrix, column_metric, metric_bent_c1,
beltrami_hessian_symbols, cao_odd_border_beta, bas_velocity_derivative,
hf_hill_potential, hf_general_first_harmonic, hf_growth_bracket,
hf_resonant_slow_matrix, hf_detuned_return_generator,
charged_column_slow_block` (13 functions; 0109 imports the last).
`euler_cao_schur.cao_thin_ring_schur_jet` import-smoked in G2.
`euler_p2_principal.py#C40E` header re-read: algebraic identities only, no
coercivity/persistence assertion — scope preserved.

## U4 — cipher ask (2026-09-10)

`herd/inbox/beacon.md#0C7D` cipher entry 16:01Z: B2 magnitude + R/Z scaling
for PoC-1 denominator; exploratory cap accepted. Cipher 0108 reconciles
M1-needs-B2 (`03-reconciliation.md:13`); PoC-1 design at
`attempts/0108-cipher-radical/04-poc-designs.md:5-11`. Protocol:
`herd/protocol-v1.md#25BE` (v1 signal format adopted for all 0109 appends).
