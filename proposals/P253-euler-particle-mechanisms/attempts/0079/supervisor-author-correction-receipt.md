# Supervisor author-correction receipt

This append-only receipt preserves the activated README at SHA-256
`9a61ba2d5df5677540924392453149eadab760b02eb75584563c24fe6d210942`.
The frozen candidate routes did not change.

## Exposed author-stage defects and repairs

1. The preregistered and first author-stage leakage ledger contained only the
   outbound `QVP` kernel. It could bound `Q` generated from `P` but not the
   `PVQ` feedback into `P`. The final derivation uses distinct `eps_out` and
   `kappa_back`; the isolated propagator and its inverse plus return closeness
   yield the two-sided absolute action estimate (28c).
2. Harmonic conjugation was incorrectly allowed to classify a counter-rotating
   term as `V_A`. The final derivation requires the explicit positive/negative
   frequency blocks before applying `V_C=(V-JVJ)/2` and
   `V_A=(V+JVJ)/2`.
3. The first crossing candidate used fixed `n` as `delta->0`, outside 0074's
   massive regime. The final calculation uses rational rays
   `n_1=Np,n_2=Nq`, fixed positive column wave numbers, and uniform `C^0`
   graph-Riesz bracketing. A continuous exact-carrier path covering geometric
   `delta` remains explicit; it is not attributed to 0074/0078.
4. The optional `C^1` `O(N^-2)` location formula is retained only as a
   symbolic refinement. Exact existence needs only endpoint error `o(h_N)` and
   IVT, giving `delta_N=delta_N^0+o(N^-1)`.
5. The `m=0`/`m=1` polarization route was repaired to use
   `tau_1=-b_1>0` and a sign-changing boundary root. No derivative is taken
   through the value-level `O(k^2)` remainder.
6. The high mixer order is measured in the physical covariant graph:
   `delta_N ell_N=O(1)`. The conditional Fourier-volume dilution and gate-time
   scale are separated from the still-open absolute KKS-dual normalization,
   finite-row bounds and nonlinear long-time construction.
7. Completed independent 0078 review is now the authority boundary for all
   five established 0074 units. Review SHA-256 is
   `676fd487d310690cfe89246a5f0a8ba71c676c39c02aa69b4dd2a393a4b2c95e`;
   verdict SHA-256 is
   `33725bbee5e149ca152146583614a8b8cc02bcad6ec49ed9e24ea77a06c2fe3f`.

## Preserved and corrected executable provenance

Before the leakage repair:

- `schwinger_hopf_blocks.py` SHA-256
  `5fad718c67e29fe2fb33b5445596b4045665fb2bdc16ac2f459b5073f7467d5b`;
- `verify_schwinger_hopf_blocks.py` SHA-256
  `cdc037aba87214214256c0ed559b70094f621d6fed7fd9f068378b3eeacf6e40`;
- first stdout SHA-256
  `b9c783b618fdc516e120fff9b3bbc1b204673cf03f2ffedbec9dc54fec27dab7`.

Final corrected artifacts before this receipt:

- `derivation.md`
  `a2020660662dc6ce4964bd97e0a7eb3e00a33d740b1f97369ee624ef55e7c7d5`;
- `source-audit.md`
  `8a2419ecd0dfdfb605e8e99db3e8ebf2d8fb7f1e2172efa512b77bdd5feccb6d`;
- `result.yaml`
  `d401ffe21845e27716ad25d59e7f2af9a1887d6049b973d3754305d9947e19ed`;
- `validation.md`
  `f4858179ee6d57af394de7f69fb9183367fc1d7efed631094937f3060515d913`;
- `schwinger_hopf_blocks.py`
  `fc72be3b906cd6c7fae956915c674bbfb27b3bdf200df9b27415722250cae8a1`;
- `verify_schwinger_hopf_blocks.py`
  `6366d0ec7698edf8d620eb626483e3af216ca6a16a3079c0d875efbc5868debd`;
- final symbolic stdout
  `6c0de54d9bca2039472ed1c53b4b3e786da0aebb6db4aebb59b5ab3c265e3844`;
- final symbolic exit
  `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

The first artifact-validator invocation is preserved in
`validation-first.*`; it exposed only a literal `O(N^-2)` versus `O(1/N^2)`
notation assertion. The corrected final repository-interpreter execution is
in `validation-v2.*` with exit exactly `0`. No production numerics were run.
