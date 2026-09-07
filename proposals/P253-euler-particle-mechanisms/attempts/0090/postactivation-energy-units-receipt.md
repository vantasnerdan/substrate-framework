# P253/0090 postactivation energy-units correction

The first central schema activation used README SHA-256
`518ef68ca3e3f7955098ad98a57c823716825b66bee40bb85302092705f29f2f`.
Its command/stdout/stderr/exit files are preserved as
`preenergy-activation-schema.*`; their SHA-256 values are respectively
`5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

Immediately after activation, dimensional inspection of Route C found that
the finite-gate expression divided radiated energy by the physical mode action.
The README now divides by the mode energy `nu_phys A_mode`.  This changes no
route, supplier, source convention, or exclusion.  The corrected README
SHA-256 is
`073fe1a19f64ec1590a6dc2920795c9b5f865bb018519214025eb6da18587d79`.

The central schema was replayed against those corrected bytes with the same
command.  It exited exactly `0`, wrote
`WORKFLOW VALID: 271 claims, 271 accepted, 13 proposals, 4 skills; MIGRATION QUEUE: 218 units, 0 pending, 0 partial`
to stdout, and left stderr empty.  Body work performed after the first
activation consists only of the normalization API and its exposing tests; it
does not consume the corrected Route-C criterion as prior authority.
