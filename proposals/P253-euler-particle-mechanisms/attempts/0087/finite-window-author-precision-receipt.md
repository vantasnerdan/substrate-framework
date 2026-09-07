# P253/0087 bounded author precision receipt

This single claim-scope repair records the finite-ceiling consequence of the
Debye estimate without changing any equation, API, test, verifier, or route
verdict.

The pre-correction hashes were:

- `derivation.md`: `247079dd059c574ee3b96b98305d7504403106e5e5f84ea08acbc2c0d8a17f82`
- `result.yaml`: `356ebe4a5ed0b66dabb26bb6c4294896f502aa08808d586ad2d8abaa0457aabf`

The post-correction hashes are:

- `derivation.md`: `44a6d6ad277b51f7557fb0047e24ed6c0ca2cdce4ab443aaeaeaf99db2b5a8a2`
- `result.yaml`: `7c14ccf452d4d72e974855c4a4701dbb17cfb64390cc8fd18c419f6cc7e4ed5e`

The corrected statement introduces a constant-dependent lower threshold
`H_loss`.  A finite-gate suppression conclusion requires an actual integer in

    max(H_graph,H_Riesz,H_KKS,H_gate,H_asym,H_loss,...)
      < H < min(H_EM,H_rad).

Because the upper endpoint is finite, exponential domination of polynomial
losses as `H` tends to infinity does not itself prove this intersection is
nonempty.  A later bounded correction replaces the former fixed-`c_EM`
high-index-limit wording by a finite-`J` predictor plus a quantified
`C_asym/J` remainder and `J_asym` threshold.  The exact shell maximum, Debye
bound, fixed-mode ceiling, and all route verdicts are unchanged.

The unchanged focused and exact checks were not rerun because their predicates
do not encode this claim wording.
