# P253/0087 finite-subluminal-range author correction receipt

This single bounded scientific-scope correction was applied after author
completion and before independent0089 froze its content-blind target hashes.
It replaces a nonexistent fixed-`c_EM` subluminal `J`-to-infinity limit by the
finite-`J` predictor and quantified joint-transfer remainder required by the
already frozen finite-window route.

## Pre-correction hashes

- `README.md`: `8e4fd5886a6bd7ff28e51786e12639ac875df55caf3530b4980da6bc9dcea7f2`
- `derivation.md`: `44a6d6ad277b51f7557fb0047e24ed6c0ca2cdce4ab443aaeaeaf99db2b5a8a2`
- `source-audit.md`: `e29f6812f432efc995f0e3c6dd205905b1e5a75722f5994d9bd5b16dc6f2206f`
- `result.yaml`: `7c14ccf452d4d72e974855c4a4701dbb17cfb64390cc8fd18c419f6cc7e4ed5e`
- `validation.md`: `484e5ee5885ad9dd2f5f684aed4bc8c28e602edbcfc2363b6a2d6f6179647780`
- `author-completion-receipt.md`: `81bcfedf97d3e1c3432b64273ef7ddc719094d554f5b1922e0c0825cc8ac2ca1`
- `finite-window-author-precision-receipt.md`: `fe0eb811204e12cee5169473b733758b0157ba3445cdd4cfc7f7408552f7e306`
- `postactivation-precision-receipt.md`: `1c1fcbdc12e323c9841bc96e8e793491c1e0c09cad5c965ebc0e20f7987db9c6`

## Post-correction hashes

- `README.md`: `498f0592474b6d5a9297fa1a26755436e830bda8d2768708a1d48a4e7334a4ee`
- `derivation.md`: `ff2ac36934271382bee99a00cc759e0449b3d9e6fb111c18ab2cb248200ff79f`
- `source-audit.md`: `e3a40e1e970c4f78754d0db2b47e605d849d552e25294fa95cb98dca78f1e8bd`
- `result.yaml`: `bb0845624141562d62ecffd0485a1015ec8c49d62793cc77bdbbd49b80002968`
- `validation.md`: `7e5bc72716a242d8633c1af12c224d1b369bf8cb50658387727264ec3ebe2741`
- `author-completion-receipt.md`: `921d15d00888a70558aac8d1bfb24076d50660f9accf0674f1769f9181ad9e9c`
- `finite-window-author-precision-receipt.md`: `c4d4bd234005336533da6cba53290ae23558e63e32034ffc48fb1ae73cbcf41d`
- `postactivation-precision-receipt.md`: `53b813a8cc42d0db47bc159cb3ccfcdae2c674a9819da86feb6090749718d3bf`

For each actual adjacent transferred pair, the repaired statement defines

    q_0(J)=C_0/sqrt(c_EM^2-c_g(J)^2),
    Q_i(J)=q_0(J)+R_i(J),
    |R_i(J)|<=C_asym/J.

An exact Debye margin requires an actual integer with
`q_0(J)<=1-2*eta_rad` and `J>=J_asym>=C_asym/eta_rad`, as well as every other
lower threshold, below `min(J_EM,J_rad)`.  The current fixed-`k` transfer does
not supply `C_asym` or `J_asym`, so Route B remains a conditional finite-window
criterion.  The exact physical clock, shell maximum, vector Bessel orders,
Debye bound, fixed-index finite ceiling, and route verdicts are unchanged.

The repository schema replay on the corrected README is recorded by
`finite-j-activation-schema.*`; it printed `WORKFLOW VALID`, left stderr
empty, and exited exactly `0`.  The unchanged executable predicates were not
rerun.
