# P253/0100 bounded correction requested by P253/0101

The independent `0101` substantive pass preserved every scientific equation
and route verdict. It found one evidence-ledger overstatement in
`failure-derived-continuation-receipt.md`: the attempt derives the exact Cao
magnetization superpotential formula but does not prove that
`grad Phi cross grad P` is nonzero on the actual charged branch.

The single correction changes `nonzero` to `potentially nonzero`. No source,
API, test, verifier, predicate, or captured output changed, and no unchanged
oracle was rerun. A stale-language scan found no other unconditional
magnetization-nonvanishing claim in the final claim-bearing package.

Pre-correction SHA-256 values:

```text
c32d1298554a36074fb38fdde67269d67dffb0fd22c8906784bd2aef9029a3ed  failure-derived-continuation-receipt.md
53fde48bb1495daee2a4d055f504cd501a1d7d742a2b80a0ccb1f5e81a825397  author-completion-receipt.md
dff6ddccd1166019720f6a545344d246195075f9cf462833ea68886c6af8d54f  artifact-hashes.sha256
```

Post-correction SHA-256 values before manifest refresh:

```text
e44fcde0af9cb45217045cd0a94ccfcb0ca5c2f99f5c834f7cbebe123154e978  failure-derived-continuation-receipt.md
49711c3a0c6746f10a98c3af8eebaf9d839690a8a4be90c006033c53e51b4619  author-completion-receipt.md
```

The manifest intentionally does not hash itself.
