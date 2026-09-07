# P253/0080 author completion receipt

The activated README remains SHA-256
`27f71a514462a5eb9ea01cf4b88e919a55b6ea3917858b878da781cf00c4280b`.
The candidate routes and parent objective did not change.

## Interrupted-file recovery

The workstation interruption left `derivation.md` with literal leading `+`
characters, a literal heredoc terminator, and an embedded shell command plus
`result.yaml` body. That transient malformed file was never frozen, reviewed,
or used as evidence. Root rewrote the derivation and separate result file,
then checked that no leading patch markers, shell commands, or heredoc
terminators remained. Existing activation and failed/successful oracle receipts
were preserved.

## Author-stage proof repairs

The final proof establishes the projected-complement HSE as the index-zero
`Q A Q` isomorphism in derivation (17), and replaces the unnecessary full
fixed-`(mu,c)` inverse by two distinct bordered derivatives:

- `B_I=D_(zeta,mu,c)(F_0,kappa,I_z)` at the uncharged base, used for the
  charged branch at fixed circulation and physical impulse;
- `B_R=D_(zeta,mu,c)(F_0,kappa,M_1-R_*kappa)`, where
  `R_*=M_1(zeta_0)/kappa_*` is a frozen target number, used for the exact
  fixed-circulation/fixed-mean-radius external-core path.

The normalized contradiction uses the threshold scale `H_epsilon`, raw density
scale `D_epsilon=epsilon_core^(-2)H_epsilon^p`, and product norm
`||h||_Xepsilon/H_epsilon^p+|eta_0|/H_epsilon+|eta_1|/H_epsilon`.
The limiting equation exhausts every allowed-growth Lane--Emden cell: the
`m=0` scaling cell and even `m=1` radial translation remain, axial translation
is removed by parity, and higher harmonics are removed by the source kernel
classification. Exact integration by parts gives

- `integral q_a=2 Lambda_p/(p-1)`;
- `integral q_t=0`;
- `integral y_1 q_t=-Lambda_p`.

After retaining the physical Jacobian, the impulse border has diagonal
`(M_a,-2*pi*rho_m*R_*s*Lambda_p)` to leading order, while the mean-radius
border has `(M_a,-s*Lambda_p)`. The positive constant-linear affine matching
matrix excludes parameter-dominated sequences. Both bordered derivatives are
index-zero isomorphisms at the base.

The charged IFT is stated only in an explicit nonempty logarithmic window
`L_border<=L<=L_EM`, with a uniform subluminal Maxwell margin. Full
fixed-parameter `A_epsilon` invertibility, the sharp finite-core two-by-two Schur sign, quantitative
`epsilon_core`-to-`delta` coverage, the finite charged harmonic hierarchy,
full three-dimensional coercivity, and action selection remain explicit
continuations. The branch is a declared Euler--Maxwell extension with new
`U(1)`, `g`, `epsilon_EM`, and `mu_EM` inputs; it is not adopted as a bare-Euler
replacement.

## Executable provenance

The first six-predicate exact oracle exposed an implementation-only structural
comparison and exited `1`; `green-schur.*` and
`first-green-schur-failure.md` preserve it. Replacing the raw structural
comparison by exact SymPy simplification produced the first successful
six-predicate receipt in `green-schur-corrected.*`. The proof-completion oracle
then added the exact scaling/translation density moments, physical two-cell
determinant, and affine matching determinant. Its final captured receipt is
`green-bordered.*`, with stdout SHA-256
`601816ad3de2218e31d33f9fe9fcca191e8a65be06b5a25a6cd3b0a8893f3c83`
and exit exactly `0`. The focused direct-consumer receipt reports two tests
passed, exit `0`; no production numerics were run.

## Final core hashes before this receipt

- `derivation.md`: `f36d577639f07fbd649e9b97f450dafe4d67c667a82d554d62f5668dba4958f2`
- `source-audit.md`: `0f5aff991e941de0a4bc7a9e9b4638f9f57dbec71084b50d55a8bf0c4c88d4d7`
- `result.yaml`: `b4e9f4d2dd69bc322182ed3f564ce5190b85b1da9ac6019d810cc1ee114ef8e3`
- `validation.md`: `e4fbdfa1fd59bb66656ad1ffb061f0325812ee419dd60c4899e33e9828e72352`
- `verify_green_schur.py`: `575554224922341826c08138d7a4bcc1b4a91d5a8baad5dca9adee67ac19f3fa`
- `src/substrate_framework/euler_cao_schur.py`: `92a621bb6537fbb1178062a26b91309fbc1583c82aefa6486ce6cb045d9d8f44`
- `tests/test_euler_cao_schur.py`: `feaaae535265d5b8456bbb52e379f2f012f32b3252d2bc72739481be65ab4495`
