# P253/0086 postactivation wording receipt

The first activation ran against README SHA-256
`221f14d72046066a0ae3612a6617040a1fb7c4c819737a7e3c5fe41974eba82e`
and exited exactly `0` with empty stderr.  It raced one content-blind wording
repair requested before execution.  No P253/0085 body, API, test, verifier,
receipt, manifest, or captured output had been opened.

The corrected README SHA-256 is
`ecc0e615d1e9e424c276b35d402398e451353f5e115948f3d16f8989181d0d7e`.
Unit A now says `field-eliminated Schur remainder`, leaving its
finite-rank/compact/bounded classification to the independent review.  No
unit, target hash, criterion, or verdict boundary changed.

The original activation quartet is preserved under
`superseded-preprecision-activation-schema.*`.  Both the original and replay
quartets have command, exit, stderr, and stdout SHA-256 values respectively

    5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f
    9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa
    e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
    d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f

The replay on the corrected README printed `WORKFLOW VALID`, left stderr
empty, and exited exactly `0`.
