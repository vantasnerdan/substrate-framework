# P253/0088 Sol-High recovery and rederivation receipt

## Runtime audit

The coordinator verified the newest rollout record at
`/home/dan/.codex/sessions/2026/09/06/rollout-2026-09-06T12-49-28-01a07656-64cf-74e2-ac94-9ba86652daef.jsonl`,
line 13774, ordinal 13772, timestamp `2026-09-07T01:38:09.080Z`:

- model: `gpt-5.6-sol`;
- reasoning effort: `high`;
- approval policy: `never`;
- sandbox policy: `danger-full-access`;
- working directory: `/tmp/substrate-particles`.

All post-fallback text was treated as quarantined scratch until independently
rederived in that verified context.

## Quarantined starting hashes

- `derivation.md`: `27f18bcba5da45f8c5d589304e33ce8bead11ea45e57ee11d5a212d4002ed6a0`
- `result.yaml`: `32bc7085bda9a40c6bd9636d2d8497af792ecbebb0943d32b9302414656ea8f0`
- `source-audit.md`: `8bcd8fb3dcb9a8c92e164e0c1b2b2bd048cb822e4928f299c6c284a943c8c5bb`
- `validation.md`: `03ae3947dc79e04cbd2189cc735a8c453c4ac97bec5df3198b5fc117f4dc5a9c`
- `verify_hessian_response.py`: `44a6068760891513594b96a2bc3cc5d700ca825db6fd8b5ef2ce5836bf745955`

These hashes are provenance only and carry no scientific authority.

## Re-derived result

The local formula (35a.1)--(35a.3) is retained but is not used as a global
Leray representative.  The recovered proof instead:

1. cancels the column `W` factor exactly in the physical KKS row, producing
   an explicit continuous `L2` vorticity covector;
2. proves exterior-extension independence by
   `integral curl(delta a) dot Bq=integral delta a dot q=0`;
3. constructs the regular-axis/interface/decaying-`K1` adjoint Sturm field;
4. uses the algebraically simple rank-one dual Riesz range and one compact
   positive-core DA test to identify that field with the intrinsic Hessian
   row, without assuming compact-core density in the edge-degenerate closure;
5. derives the singular-axis recurrence and
   `A_1^sharp=rho_0 A_1!=0`;
6. proves `gamma_12^col>0` with the explicit punctured-annulus DA seed; and
7. combines that result with independently reviewed 0083/0089 to establish
   source-specific unnormalized nonzero response on sufficiently thin exact
   same-family Cao crossings.

The uniform physical `Y4` high-index response scale, autonomous histories,
two-sided gate control, P2, and P4 remain open.

## Final author hashes

- frozen `README.md`: `6b01cbe08fc0321610b64315725e46c6c582be7592939d58b9795a44c5a65393`
- `derivation.md`: `f9852decae5879980f2d93d3ba3f85213ef61a5956c4308302fd95e68c89edd7`
- `result.yaml`: `7f7c5ac3e90c0feb6f6cb2fb9c64cc08fa59341cf3344f7bed106c224c9edf3a`
- `source-audit.md`: `dd578cfd4b482cb282fb57fb6dd0074604363e828e15395008a3e3fe018dc797`
- `validation.md`: `ce0e7e216eb86f3444b29a8a9660747bb53af05edac8df1435a7d25725d3e610`
- `verify_hessian_response.py`: `89f8d3d085a14b651fd96c9157ab22b8c62bbf6f2566e47365cb650df25f1021`
- `validate_attempt.py`: `dcb252e00bc0871228942c04a2b7404454d4e1c2109ecf0368c8492a45adaf6f`
- `independent-review-target.md`: `0af89036bb4c2e87fbf3d6c77795036438c86944a3af9631f231cb1d6730d34c`
- final exact-oracle stdout: `2afe1c43a0ce04c24e68a54a28f2e414587b406cf3790f5c3152136ce3e38525`
- final static-validation stdout: `534b8afc0b191b944b40858cb0c6e2b9254fcd717f3114409ce9ec4b0d5bef7e`

The materially changed exact oracle used the repository interpreter, exited
zero, and has empty stderr.  Its preceding structural-equality implementation
failure is preserved in `hessian-oracle-dual-riesz-first.*`.
