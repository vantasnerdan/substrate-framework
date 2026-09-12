# P253/0095 checkpoint whitespace receipt

`git diff --cached --check` exposed one surplus blank line at EOF in each of
three files after staging.  This checkpoint correction removes only that
extra newline, leaving exactly one terminal newline and changing no source
token, scientific statement, API behavior, or verifier predicate.

```text
3d2eefa07163fbc808ee2dc820242243af470ee953be168259c353434222ff5c  ->  e950c7a2684b9813e33f6d35855bbc9dd1bbea41fe0ba3caebc903e52c2c2ead  packaging-test-failures.md
03cd19b896aeb46591fa7bcfc7b59805f279da39359807d3780a544565092334  ->  aee01bd519897b09191b72e2e9f50e7204180d17b3fbabd259f3489d9f70274f  recovery-audit-6-receipt.md
c22eb7847075c4410f37c69cea4d336c342ed7610c83c0adcd3239ee624b97f9  ->  7087cb2edf4ea8d9b6ea71ca027f34f117ade2cf269430b3057aee799fdc2a84  ../../../../src/substrate_framework/euler_p2_principal.py
```

The pre-cleanup manifest SHA-256 was
`dc15211dd0a2b4799a72c89446592202cf3256898f896118ad31915067bfaba4`.
No oracle or test was rerun because executable tokens and predicates are
unchanged.  The refreshed manifest records this receipt and all post-cleanup
hashes.
