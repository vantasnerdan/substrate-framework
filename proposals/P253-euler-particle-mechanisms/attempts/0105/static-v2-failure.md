# Static-v2 assertion failure

`static-validation-v2` exited 1 because its prose scan searched for the
contiguous phrase `only nontrivial regular full-core cancellation family`,
while the derivation places a line break between `only` and `nontrivial`.
Every scientific predicate, exact-v2 check, focused-v2 test, YAML parse and
repository validation passed.  The repaired static-v3 check searches the
invariant equation and verdict fields rather than a line-wrapping-sensitive
phrase.  No scientific claim or implementation changed in response.
