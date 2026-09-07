# P253/0090 author completion receipt

The prescribed-current calculation is complete at its declared author scope.
The strongest exact result is

    P_out=pi/(2 epsilon_EM)
      integral |Omega| delta(c_EM^2|k|^2-Omega^2)
                    |P_T J_hat|^2 dk,

with the exact coarea sphere weight, conserved switching completion, and
late-time free-field observable derived in `derivation.md`.  For a reciprocal
mode, the amplitude width is `P_out/(2 nu_phys A_mode)` and the fractional
energy-decay rate is `P_out/(nu_phys A_mode)`.  No source-specific Cao current,
resonance, P4/P5, electron, neutrino, or action/statistics selection is claimed.

Pre-review/pre-correction core SHA-256 values:

```text
073fe1a19f64ec1590a6dc2920795c9b5f865bb018519214025eb6da18587d79  README.md
20c882b24bc788ea4bb6c962d89b789efbec1e19300b0a35e76cfd5c4a11e802  derivation.md
dc1450eb9ec25b14797a3f04c30c1ecad7650e42b4edbdef29c03f6a2314f60b  source-audit.md
f64defb0400dc348ff68ced7356630cbec4717677a307728262d47230e66bd25  result.yaml
869722796ad79c8d25437d0b85f1765878a539a5ad6ce32a933da47ca53e403b  validation.md
b536d374ab54c9a7f54931a3644334162bc6627ea294aa0e8514faa4a74522d1  verify_shell_flux.py
74881415a94e53350d2db5edb328a0d8611d3eace9b1c48fb680f0c0db14fb4d  src/substrate_framework/euler_maxwell_radiation_flux.py
1bc141cc20da1f5cb8c628a3f714da51e6c0b1428fbb5c33fd50304907d618d9  tests/test_euler_maxwell_radiation_flux.py
```

The exact verifier emitted nine `PASS` lines, exit `0`, empty stderr.  The
focused public API run reported six passed tests, exit `0`, empty stderr.  The
repository validator exited `0` with the standard `WORKFLOW VALID` line and
empty stderr.  The final supervisor normalization audit independently matched
equations (3)--(18) and requested only the explicit full-real-current
definition of `F_T`, now present in equation (22); no oracle rerun was needed
for that prose precision.

Independent review0091 subsequently requested one bounded analytic repair to
the plateau theorem.  The final derivation now defines the exact pushforward
density and full-real endpoint/conjugate correlations, proves the rectangle
first-moment remainder, types the shifted-anisotropic Helmholtz LAP bridge and
literal neutrality domain, and scopes equation (25) to its finite moving-
volume identity.  The correction is recorded in
`0091-bounded-correction-receipt.md`; the unchanged exact and focused oracles
were not rerun.

Final post-correction claim-bearing SHA-256 values:

```text
c2c39e2ef0585647dfe3d3238c9c2999cf6bbe970389df8528706697985b0360  derivation.md
43140e6eecfe2411fd7b31f550eaba6416ac040d6babdbebe1077c74305a789e  source-audit.md
4a5ff2e67ea6d3713dede6fd269eb4402ac8d5e9c1973c338376546578f87b3c  result.yaml
a5cd1523bc6c653d5902577a68159de8196645125f173517de1927323436a2ea  validation.md
```
