# P253/0077 author completion receipt

The activated README boundary remains SHA-256
`0ee7c78850839ab77030a00cc03b399f91dcfeb6a6e0668b10ff54d1938bba35`.

The exact steady-map oracle first exposed an implementation-only stabilizer
assumption error and exited `1`; the command, stdout, stderr, exit, and diagnosis
are preserved under `first-exact.*` and `first-exact-failure.md`.  Removing the
incorrect `nonzero=True` assumption from the unknown radial displacement gave
the first successful run.  The later substantive Unit-A correction replaced
the false local-in-time inference by exact period integration and strengthened
the oracle accordingly.  The final repository-interpreter execution is
`periodic-corrected.*`, exits `0`, and passes all eight exact checks.

Final core SHA-256 hashes before this receipt and manifest refresh:

- `derivation.md`: `28a137fdc6b157895d89508aa6378b344a005e4255f7ed010d80f0efe6df0ae2`
- `source-audit.md`: `75844ab749d45575a322f7f5a4a24e6e275203f7343e065886a25226bc6c4383`
- `source-receipt.yaml`: `cef1ecd0260396e1eb5722fafcbca0d456ee0ea1b08130cd977801509dc3a18f`
- `result.yaml`: `9c6130d3730de0f2034aeb563414268b129952607ec4fa6af46517201feb8033`
- `validation.md`: see the append-only P253/0081 evidence-ledger correction
  receipt for the final post-review hash;
- `verify_charged_steady_map.py`: `47feb24f2ea8bad37f1533bce07d99313b105c8ba6f59675f5b0cd12d0316c17`
- `periodic-corrected.stdout.txt`: `38249ca4c0d3bdbeca674a5dad8601c293459b82324f5afda32272c259af11d7`

The strongest established results are the exact-period fixed-tail
compatibility theorem and the full subluminal axisymmetric no-swirl
Euler--Maxwell steady reduction with action-tag stabilizer locking.  The exact
charged Cao branch remains conditional on two named achievements: the exact
steady complement theorem `HSE` on fully declared axis-weighted spaces and the
finite-core two-by-two circulation--impulse Schur determinant `S_epsilon`.
They are the next analytic attempt, not hidden claims here.  No production
numerics were run.

P253/0081 requested one bounded evidence-role correction after opening this
package.  Only `validation.md` changed: it now lists the eight predicates the
oracle actually checks and leaves the traveling-potential identity and full
Ampere primitive to the analytic derivation.  The unchanged oracle was not
rerun.  Exact pre/post hashes are recorded in
`review-correction-0081.md`.
