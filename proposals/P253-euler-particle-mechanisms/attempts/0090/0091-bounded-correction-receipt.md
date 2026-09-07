# P253/0091 bounded correction receipt for P253/0090

This is the single correction transaction requested by the independent 0091
review.  It changes no shell, coarea, Gaussian-source, switching-current, or
reciprocal factor-two formula and does not rerun the unchanged oracle/API.

The repair supplies the previously compressed analytic bridge for
`E_late(T)=T P_out+O(1)`:

- the resonant pushforward density is defined with its exact `1/(8 epsilon)`
  real-phasor/late-energy coefficient;
- the rectangle autocorrelation identity gives a remainder bounded by twice
  the first correlation moment;
- fixed endpoint profiles include every ramp and `a_prime*K` term, with
  endpoint energy and rectangle--endpoint correlations separately bounded;
- the conjugate temporal branch has detuning at least `|omega|`, and its cross
  correlation is included explicitly;
- neutrality now requires a literal zero mean or enough integrability to
  evaluate continuity at zero;
- the comoving symbol is completed to a shifted anisotropic Helmholtz symbol,
  typing the `L2_1 -> L2_-1` outgoing bridge; and
- equation (25) is scoped to the exact finite moving-volume balance.  Its
  infinite-world-tube equivalence with late free-field energy remains open.

Pre/post SHA-256 values:

```text
derivation.md
  pre  20c882b24bc788ea4bb6c962d89b789efbec1e19300b0a35e76cfd5c4a11e802
  post c2c39e2ef0585647dfe3d3238c9c2999cf6bbe970389df8528706697985b0360
source-audit.md
  pre  dc1450eb9ec25b14797a3f04c30c1ecad7650e42b4edbdef29c03f6a2314f60b
  post 43140e6eecfe2411fd7b31f550eaba6416ac040d6babdbebe1077c74305a789e
result.yaml
  pre  f64defb0400dc348ff68ced7356630cbec4717677a307728262d47230e66bd25
  post 4a5ff2e67ea6d3713dede6fd269eb4402ac8d5e9c1973c338376546578f87b3c
validation.md
  pre  869722796ad79c8d25437d0b85f1765878a539a5ad6ce32a933da47ca53e403b
  post a5cd1523bc6c653d5902577a68159de8196645125f173517de1927323436a2ea
```

Unchanged evidence hashes:

```text
b536d374ab54c9a7f54931a3644334162bc6627ea294aa0e8514faa4a74522d1  verify_shell_flux.py
74881415a94e53350d2db5edb328a0d8611d3eace9b1c48fb680f0c0db14fb4d  src/substrate_framework/euler_maxwell_radiation_flux.py
1bc141cc20da1f5cb8c628a3f714da51e6c0b1428fbb5c33fd50304907d618d9  tests/test_euler_maxwell_radiation_flux.py
```

The completion receipt and refreshed artifact manifest are outside this
receipt's self-reference set and are updated after this file is hashed.
