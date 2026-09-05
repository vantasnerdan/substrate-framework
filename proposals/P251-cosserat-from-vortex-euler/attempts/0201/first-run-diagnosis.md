# First execution and bounded implementation repair

The first execution evaluated and passed all 30 scientific predicates, then
exited 1 because the verifier called nonexistent `CheckLedger.summary()`.
The shared API actually exposes `finish()`, as the frozen 0195 consumer
already demonstrates. The repair changes that terminal method call only.
It changes no equation, gate, source profile, operator or physical observable.
`first.stdout` is retained; the corrected execution is `repaired.stdout`.

This is a verifier implementation failure, not a scientific route verdict.
