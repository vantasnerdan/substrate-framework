# Post-activation supervisor precision correction

Date: 2026-09-06

Central activation occurred against README SHA-256
`80e39cc34576876da61c3d713ec8f0e2b6d78fd8fdedcb905455d537eb159834`.
Before any 0066 claim body was written, the supervisor identified four frozen-
contract precision issues. The README was corrected without changing its P2
objective, three routes, route selection, or success boundary. Its corrected
SHA-256 after those four edits was
`773674c9360af3a1d14d7871ec2c47c0ef6837241d32fa857e01497c3b785d7f`.
One stale topology bullet still described the edge as
`zeta~kappa d^p`. It was then replaced by the reviewed
`zeta=O(d^p)` upper vanishing, with a two-sided `d^p` weight expressly
conditional on proving the nonzero `kappa`/comparability output. The final
corrected README SHA-256 is
`d9b3ffcd049cfea6812ccb3558a958251ada118bea1d4577a02334b2107311fc`.

The correction (1) distinguishes independent-review 0057 from author-scope
0058/0062 inputs; (2) uses only the reviewed `zeta=O(d^p)` collar bound and
makes nonzero leading coefficient/two-sided comparability an output; (3)
freezes `p>=6`, `s=4`, hence `H^5` displacement and `H^3` ambient vorticity;
and (4) fixes the translating nonrotating frame `Omega_0=0`, records exact
`n`-block rotation-shift covariance, and defines `E_2` as the relative Hessian
of `H-c_0 P_z` on the constrained coadjoint leaf.

Post-correction schema replay used
`/home/dan/substrate-framework/.venv/bin/python scripts/validate_repository.py`.
It returned exit `0`, empty stderr, and stdout:

    WORKFLOW VALID: 271 claims, 271 accepted, 13 proposals, 4 skills; MIGRATION QUEUE: 218 units, 0 pending, 0 partial

This is a transparent post-activation contract correction receipt, not a
scientific result or a second activation claim.

After `activation-schema.exit=0` was materialized, one stale status sentence
still described the attempt as awaiting activation. The sentence alone was
updated to state the existing central activation; no route, target, premise,
or success condition changed. The resulting README SHA-256 is
`472119ddb209687c16bc6e3af9a226e6e94e5ca69a8739e55f2a548470db5c0f`.
